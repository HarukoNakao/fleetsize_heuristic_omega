######################################
# code for solving the fleet size problem problems 
# MILP model formulation : see the pdf
# Tai-Yu MA,  11-15-2022
# Please do not diffuse without permission
######################################

using DelimitedFiles, DataFrames, CSV, XLSX
using Random, Distributions, Distances
using LinearAlgebra, Combinatorics 
using JuMP, StatsBase, Plots , Gurobi
using Profile #see https://www.julia-vscode.org/docs/dev/userguide/profiler/

const DIGITS = 3 #  precision digits after the decimal point for round()
const GRB_ENV = Gurobi.Env()  # to silence the solver

include("data_structure.jl")
include("read_data_fleet.jl") 
include("feasibility.jl") 
include("utils.jl")
include("generate_init_solution.jl")
include("utility.jl")
include("local_search.jl") 
include("darp_solver.jl")
include("utils_chg_scheduling.jl")

function generate_demand_scenario(scenarios_demand, demand_data, global_parameter, instance)

    n_c, n_bs, n_departure_train= instance.n_c, instance.n_bs, instance.n_departure_train
    n_scenario, max_passenger_request = global_parameter.n_scenario, global_parameter.max_passenger_request
    for s in 2:n_scenario
        demand_tmp = deepcopy(demand_data)
        demand_tmp.id_bus_stop = [rand(1:n_bs) for _ in 1:n_c]
        demand_tmp.id_train = sort!([rand(1:n_departure_train) for _ in 1:n_c])
        demand_tmp.first_last_mile = [rand(0:1) for _ in 1:n_c]
        n_station_used= demand_tmp.n_station_used
        if  n_station_used > 1
            demand_tmp.id_train_station = [rand(1:n_station_used) for _ in 1:n_c]
        end 

        demand_tmp.no_person=[rand(1:max_passenger_request) for _ in 1:n_c]       
        scenarios_demand[s] = demand_tmp
    end
end

function init_demand_scenario(n_c, n_ts, input_request::String, nodes_info, scenarios_demand)

    requestInfo = DataFrame(CSV.File(input_request))
    id_request =collect(1:n_c)
    id_bus_stop = requestInfo.id_bus_stop 
    no_person= requestInfo.No_person 
    id_train = requestInfo.id_train 
    first_last_mile = floor.(Int32, requestInfo[:, "is_last"]) # 0 first mile and 1 last mile
    id_train_station = requestInfo.id_train_station 
    demand_data = Demand_data(n_c, n_ts, id_request,id_bus_stop,no_person,id_train,first_last_mile,id_train_station)
    scenarios_demand[1]= demand_data
    return demand_data
end

function setup_chargers(total_chargers, fast_chg) 

    # setup one charging station with two DC fast
    n_chargers = sum(total_chargers)
    for cs in 1:2
        for _ in 1:total_chargers[cs]
            add_charger(fast_chg, cs) #DC fast 
        end
    end
    idx_chg=  fast_chg.set_id_chg_installed
    @show(fast_chg.n_fast_chg_installed, fast_chg.num_chg_installed_by_site, idx_chg, 
          fast_chg.pw[1:n_chargers], fast_chg.cost_fast_chg)

end

# get the solution (cs configuration) with the minimum obj_l1 value and its fleet composition
function output_fleet_solution(idx_best, vec_obj_l1, vec_obj_l2, vec_sol_scenario, instance, comnbis, global_parameter, rho)

    n_request= instance.n_c
    n_scenario=global_parameter.n_scenario
    co2_reduc_target = global_parameter.co2_reduc_target
    obj_l1 = vec_obj_l1[idx_best]
    obj_l2 = vec_obj_l2[idx_best]
    vec_tmp = vec_sol_scenario[idx_best] 
    vec_obj = vec_tmp[:,1] 
    threshold = max(1, floor(Int32, rho * n_scenario))
    idx_threshold =   sortperm(vec_obj)[threshold]
    fleetsize_gv = vec_sol_scenario[idx_best][idx_threshold, 2]
    fleetsize_ev = vec_sol_scenario[idx_best][idx_threshold, 3]
    cs_config    = comnbis[idx_best]   
    print("Solution of the fleet size problem :\n")
    @show(instance.instance_name, n_request, n_scenario, co2_reduc_target, obj_l1, obj_l2, fleetsize_gv, fleetsize_ev, cs_config, rho)
    # return Fleet_size_solution(n_request, n_scenario,co2_reduc_target, obj_l1, obj_l2, fleetsize_gv, fleetsize_ev, cs_config)

end

function fleet_size_solver(path, silent_fleet_size_solver, silent_darp_solver, global_parameter, darp, fast_chg, instance, scenarios_demand, demand_data, nodes_info, timeTable, output_solver)
     
    n_scenario, rho, co2_reduc_target= global_parameter.n_scenario, global_parameter.rho, global_parameter.co2_reduc_target
    parameter_energy = instance.parameter_energy
    global_parameter.flag_initialization = false 
    n_chg_station = nodes_info.n_chg_station
    # generate n_scenario and store them in a vector
    generate_demand_scenario(scenarios_demand, demand_data, global_parameter, instance) # generate new demand scenario
    # suppose there are max two charging station sites 
    n_chg_station > 2 && error(" nodes_info.n_chg_station !, please revise fleet_size_solver() ")
    n_cs_i, n_cs_j= fast_chg.max_n_chg_by_site[1], fast_chg.max_n_chg_by_site[2]
    comnbis =vec([(i,j) for i in collect(0:n_cs_i), j in collect(0:n_cs_j)])
    sort!(comnbis, by = x -> x[1]+x[2] ) #sorted accoriding the num of installed chargers
    n_combis = size(comnbis)[1]
    @show(n_combis, comnbis)
    # bilevel optimization 
    vec_obj_l1 = zeros(Float32,n_combis);vec_obj_l2 = zeros(Float32,n_combis)
    vec_cost_chg_infra = zeros(Float32,n_combis)
    vec_sol_scenario = Vector{Matrix{Float32}}(undef, n_combis)
    vec_sol_l2= zeros(Float32, n_scenario, 3) # total_cost, fleetsize_gv, fleetsize_ev
    count_cs_setting_tested = 0
    vec_obj_l1_n_chg =zeros(Float32, n_cs_i+ n_cs_j + 1) # shift 1 
    # for idx_cs_setting in 1:n_combis  
    for idx_cs_setting in 1:3  

        t_start_cs_setting=time()
        count_cs_setting_tested += 1
        vec_sol_l2 .*= 0     
        for cs =1:n_chg_station
            for _ =1:comnbis[idx_cs_setting][cs]
                add_charger(fast_chg, cs)
            end           
        end
        @show(idx_cs_setting, fast_chg.set_id_chg_installed)
        sum_obj_l2=0
        vec_cost_chg_infra[idx_cs_setting] = get_cost_chg_infra(fast_chg, global_parameter, n_chg_station)
        ! silent_fleet_size_solver && @show(idx_cs_setting, fast_chg.n_fast_chg_installed, fast_chg.num_chg_installed_by_site, fast_chg.set_id_chg_installed)
     
        for scenario in 1:n_scenario 
            ! silent_fleet_size_solver && @show(scenario)
            demand_data = scenarios_demand[scenario]             
            darp, _ = generate_darp_instance(demand_data, instance, global_parameter, nodes_info, timeTable, fast_chg, darp)
            global_parameter.is_gasoline_fleet = true # get the co2 emission of the gasoline fleet
            global_parameter.co2_threshold = global_parameter.co2_threshold_init  # big_value
            t_start_gaso=time()
            best_solu, global_parameter.co2_gasoline_fleet = darp_solver(path, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver)  
            global_parameter.co2_threshold = global_parameter.co2_gasoline_fleet * (1-co2_reduc_target) 
            t_elapse_gaso = time()- t_start_gaso 
            @show(t_elapse_gaso, global_parameter.co2_gasoline_fleet, global_parameter.co2_threshold )
            # solve the co2 emission constrained mixed fleet DARP with capacitated chg stations
            global_parameter.is_gasoline_fleet = false
            t_start_ev=time()
            best_solu, _ = darp_solver(path, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver) 
            t_elapse_ev = time()- t_start_ev  
            sum_obj_l2 += best_solu.total_cost
            n_route = best_solu.n_route
            # store the relevant info of the solution: obj_l2, fleet of gv, fleet of ev
            fleetsize_ev = sum(parameter_energy.is_electric[best_solu.RI[1:n_route, 5]]); fleetsize_gv =  n_route-fleetsize_ev          
            vec_sol_l2[scenario,:]=[best_solu.total_cost, n_route-fleetsize_ev, fleetsize_ev]
            ! silent_fleet_size_solver && @show("sol_co2_constraint", best_solu.total_cost, best_solu.co2_emission, best_solu.total_chg_time, n_route, fleetsize_gv, fleetsize_ev, 
            best_solu.co2_emission, global_parameter.co2_threshold, t_elapse_ev)
        end 
        vec_sol_scenario[idx_cs_setting] = vec_sol_l2
        vec_obj_l2[idx_cs_setting] = sum_obj_l2/n_scenario
        vec_obj_l1[idx_cs_setting] = vec_cost_chg_infra[idx_cs_setting] + vec_obj_l2[idx_cs_setting]
        reset_cs(fast_chg, n_chg_station) # empty all cs 
        ! silent_fleet_size_solver &&  @show(idx_cs_setting, vec_obj_l1[idx_cs_setting], vec_obj_l2[idx_cs_setting] ) 
        # if min_obj_l1[n-1 chargers installed] > min_obj_l1[n-2 chargers installed], stop
        n_chg_installed = sum(comnbis[idx_cs_setting]) + 1 # need to shift one as 0 chargers at the begining
        if vec_obj_l1_n_chg[n_chg_installed] != 0 #update with the best solution
            if vec_obj_l1[idx_cs_setting] < vec_obj_l1_n_chg[n_chg_installed]
                vec_obj_l1_n_chg[n_chg_installed] = vec_obj_l1[idx_cs_setting] 
            end
        else
            if n_chg_installed > 2 && vec_obj_l1_n_chg[n_chg_installed-1] > vec_obj_l1_n_chg[n_chg_installed-2]
                break
            else
                vec_obj_l1_n_chg[n_chg_installed] = vec_obj_l1[idx_cs_setting] 
            end           
        end 
        # @show(fast_chg.n_fast_chg_installed,fast_chg.num_chg_installed_by_site, fast_chg.set_id_chg_installed)    
        t_elapse_cs_setting = time()- t_start_cs_setting  
        @show(t_elapse_cs_setting, idx_cs_setting)   
    end
    @show(count_cs_setting_tested, comnbis[1:count_cs_setting_tested],vec_obj_l1[1:count_cs_setting_tested], vec_obj_l1_n_chg)
    idx_best = sortperm(vec_obj_l1[1:count_cs_setting_tested])[1] 
    return output_fleet_solution(idx_best, vec_obj_l1, vec_obj_l2,vec_sol_scenario, instance, comnbis, global_parameter, rho)

end


# # function to test mixed fleet darp with co2 and charging station capacity constraints
# function darp_with_co2_constraint(path, silent_darp_solver, total_chargers, vec_co2_threshold, global_parameter, darp, fast_chg, instance, demand_data, nodes_info, output_solver)
    
#     if darp.n_cus > 100
#         idx= 9 + floor(Int32, darp.n_cus/100)
#     else
#         idx= floor(Int32, darp.n_cus/10)
#     end
#     idx > length(vec_co2_threshold) && error("idx > length(vec_co2_threshold) error!!, instance and vec_co2_threshold not match!!")
#     co2_threshold = vec_co2_threshold[idx]
#     global_parameter.co2_threshold = co2_threshold
#     global_parameter.is_gasoline_fleet = false
#     best_solu, _ =  darp_solver(path, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver) 
#     return best_solu
# end


# # function to test mixed fleet darp with co2 and charging station capacity constraints
# function darp_with_co2_constraint_batch(path, silent_darp_solver, total_chargers, vec_co2_threshold, global_parameter, darp, fast_chg, instance, demand_data, nodes_info, output_solver)
    
#     if darp.n_cus > 100
#         idx= 9 + floor(Int32, darp.n_cus/100)
#     else
#         idx= floor(Int32, darp.n_cus/10)
#     end
#     idx > length(vec_co2_threshold) && error("idx > length(vec_co2_threshold) error!!, instance and vec_co2_threshold not match!!")
#     co2_threshold = vec_co2_threshold[idx]
#     global_parameter.co2_threshold = co2_threshold
#     # setup one charging station with two DC fast
#     n_chargers = sum(total_chargers)
#     for cs in 1:2
#         for _ in 1:total_chargers[cs]
#             add_charger(fast_chg, cs) #DC fast 
#         end
#     end
#     idx_chg=  fast_chg.set_id_chg_installed
#     @show(co2_threshold, fast_chg.set_id_chg_installed)
#     global_parameter.is_gasoline_fleet = false 
#     best_solu, _ =  darp_solver(path, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver) 
#     return best_solu, output_solver
# end


# # function to test mixed fleet darp with co2 and charging station capacity constraints
# function darp_with_co2_constraint(path, silent_darp_solver, total_chargers, vec_co2_threshold, global_parameter, darp, fast_chg, instance, demand_data, nodes_info, output_solver)
    
#     if darp.n_cus > 100
#         idx= 9 + floor(Int32, darp.n_cus/100)
#     else
#         idx= floor(Int32, darp.n_cus/10)
#     end
#     idx > length(vec_co2_threshold) && error("idx > length(vec_co2_threshold) error!!, instance and vec_co2_threshold not match!!")
#     co2_threshold = vec_co2_threshold[idx]
#     global_parameter.co2_threshold = co2_threshold
#     # setup one charging station with two DC fast
#     n_chargers = sum(total_chargers)
#     for cs in 1:2
#         for _ in 1:total_chargers[cs]
#             add_charger(fast_chg, cs) #DC fast 
#         end
#     end
#     idx_chg=  fast_chg.set_id_chg_installed
#     @show(co2_threshold, fast_chg.n_fast_chg_installed, fast_chg.num_chg_installed_by_site, idx_chg, 
#           fast_chg.pw[1:n_chargers], fast_chg.cost_fast_chg)
#     global_parameter.is_gasoline_fleet = false
#     best_solu, _ =  darp_solver(path, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver) 
#     return best_solu
# end

