
########################################################
# solve the lower level problems
# 12.2.2024
########################################################
include("main.jl")
 
step_size_rm_route = 250 
is_gasoline_fleet  = false
percentage_full    = 0.5 # 0.5 means the e_init is 50% of e_max, set 1 if it's gasoline fleet 
co2_reduc_target   = 0.9# 0.9 means reduce 50% of the target,   set 0 if it's a gasoline fleet  
n_run              = 5
flag_dev           = true   # to export V_DF the detail for verification
set_seed_number    = false # set fixed (or not) seed number  
silent_darp_solver = false
path               = pwd()*"\\output\\"
path_result        = pwd()*"\\result\\"
vec_res= zeros(45,15); count_0 = count_1 = count_2 = 0
timeTable =0; tmp_res= zeros(7); 
vec_n_cus = vcat(collect(10:10:100),[200,300,400,500])

for id_set_data in 1:1
    count_1 += 1; count_2=0
    str_dataset     = "\\data_fleet_09022023_" * string(id_set_data)  # used dataset 
    path_data       = pwd()* str_dataset * "\\"
    input_tt        = path_data * "Timetable_6h_23h.csv"
    input_parameter = path_data * "Parameters.csv"   
    input_xy        = path_data * "Coordinate_6h_23h.csv" # coords of bus stops, transit station and the depot is one the last line
    input_vehicle   = path_data * "VehicleInfo.csv"
    input_cs        = path_data * "CSInfo.csv"   # please note that this data set contain additional information about the chargers costs
    for id_data in 5:5
        count_0 += 1
        count_2 += 1
        n_cus= vec_n_cus[id_data]
        instance_name   = "c_" * string(n_cus) * ".csv"
        input_request   = path_data * instance_name
        @show(input_request)
        global_parameter = parameter_settings(set_seed_number, percentage_full, step_size_rm_route)  
        set_co2_threshold(id_set_data, n_cus, co2_reduc_target, is_gasoline_fleet, n_run, global_parameter) 
        instance, darp, scenarios_demand, nodes_info, timeTable, fast_chg, output_solver, demand_data = initialization(global_parameter, input_tt, input_parameter, input_request, input_xy,input_vehicle, input_cs, instance_name)
        # @show(input_request, parameter_energy.max_ec[1:30])
        # @show(instance.parameter_energy.E_init) 
        best_solu, output_solver = darp_solver_batch(path_result, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver)  
        n_route = best_solu.n_route ; veh_ids = best_solu.RI[1:n_route,5] 
        vec_veh_type = instance.parameter_energy.veh_type[veh_ids] 
        n_type_veh = global_parameter.n_type_veh
        # @show(global_parameter.n_type_veh, instance.parameter_energy.veh_type)
        vec_n_used_veh_by_type = zeros(Int32, n_type_veh) 
        for r in 1:best_solu.n_route
            vec_n_used_veh_by_type[vec_veh_type[r]] += 1
        end
        tmp_res[1:2], tmp_res[3:3+n_type_veh-1] = [output_solver.avg_obj, output_solver.obj_best_sol], vec_n_used_veh_by_type
        tmp_res[3+n_type_veh:end] = [output_solver.co2_best_sol, output_solver.charging_time_best_sol, output_solver.t_cpu_run]
        vec_res[count_0,1:2] = [count_1,count_2]
        vec_res[count_0, 3:2+length(tmp_res)] = tmp_res
        writedlm( path_result * "res_all_0.5_0.9_dataset1.csv",  vec_res, ',')
    end
end 
  
# @profview darp_solver_batch(path_result, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver)  
        
