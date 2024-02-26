################################################################
# use data set 3 to calibrate the algorith's parameters
# code to calibrate penalty, penalty_veh, n_iter_max, n_max_stagnant,step_size_remove_route
# 12.2.2024
# ################################################################

include("main.jl")
step_size_remove_route = 250
silent_darp_solver     = false; 
co2_reduc_target       = 0.5#  reduce XX% CO2 
percentage_full        = 0.5 # percentage full of initial battery level
flag_dev               = true# to export the detail for verification
set_seed_number        = false # set fixed (or not) seed number  
path                   = pwd()*"\\output\\"
path_result            = pwd()*"\\result\\cali\\"
global_parameter = parameter_settings(set_seed_number, percentage_full, step_size_remove_route)  
global_parameter.n_run             = 3
global_parameter.is_gasoline_fleet = false
global_parameter.co2_reduc_target = co2_reduc_target # % of co2 reduction, user-defined parameter 
vec_co2_max_milp = zeros(Float32,14) 
vec_penalty_passenger      = [1, 5, 10, 15, 20,40] #penalty
vec_penalty_veh            = [50,100,150,200,250,300] #penalty_per_veh_used
vec_iter_max               = ([25,50,100,150,200,300] .* 1000) #N_ITER 
vec_n_stag                 = [25,50,100,150,200,250] #max_stagnant_multiplier
vec_step_size_remove_route = [50,100,150,200,250,300] #step_size_remove_route
# vec_n_imp                  = [100,200,300,400,500,600]
# vec_T_red                  = [100,200,300,400,500,1000]
# vec_t_max                  = [0.3,0.6,0.9,1.2,1.5,1.8]

# co2 for all gasoline vehs of dataset 2 using the reference parameters
# need to update for using dataset (to do)
vec_co2_max_milp[1:5]=[21.89301737,	36.59451329, 62.94336657, 75.64283257, 95.88865796] # max co2 emission for c10 to c50 instance based on the milp solutions, used for testing the darp solver with co2 constraints
vec_co2_max_milp[6:14]= [107.1183484, 122.9277477,	139.096208,	163.1255187, 186.6766963, 309.4690702, 458.3253716, 570.2280747, 716.4381396] # obtained from SA using all gv
vec_co2_threshold = vec_co2_max_milp.*(1-co2_reduc_target); # as above
 
timeTable =0; total_chargers = [2, 0]
 
str_dataset     = "\\data_fleet_09022023_2"  # used dataset 
path_data       = pwd()* str_dataset * "\\"
input_tt        = path_data * "Timetable_6h_23h.csv"
input_parameter = path_data * "Parameters.csv"   
input_xy        = path_data * "Coordinate_6h_23h.csv" # coords of bus stops, transit station and the depot is one the last line
input_vehicle   = path_data * "VehicleInfo.csv"
input_cs        = path_data * "CSInfo.csv"   

n_instance = 10; count_1 = count_2 = 0
# select the parameters to be tested here
vec_para_test = vec_penalty_veh

length_vec_para_test = length(vec_para_test) 
res_all = zeros(Float32, n_instance* 6, 14) 
res_avg = zeros(Float32, length_vec_para_test, 4) 
res_avg[:,1]= vec_para_test
@show(vec_para_test) 
# global_parameter.parameter_sa.N_ITER=  100*1000
# global_parameter.parameter_sa.max_stagnant_multiplier=  200

for ite = 8:10
    count_1 = count_2 = 0
for id_data in 1:n_instance
    tmp = id_data*10
    instance_name   = "c_" * string(tmp) * ".csv"
    input_request   = path_data * instance_name
    @show(input_request)
    count_2 = 0   
    global_parameter.co2_gasoline_fleet = vec_co2_max_milp[id_data]
    global_parameter.co2_threshold_init = vec_co2_max_milp[id_data] 
    global_parameter.co2_threshold      = vec_co2_threshold[id_data] 
    @show(id_data, global_parameter.is_gasoline_fleet, global_parameter.co2_reduc_target, global_parameter.co2_gasoline_fleet,
    global_parameter.co2_threshold_init , global_parameter.co2_threshold)
    for val_tested in vec_para_test
        count_1 += 1 ; count_2 += 1
        instance, darp, scenarios_demand, nodes_info, timeTable, fast_chg, output_solver, demand_data = initialization(global_parameter, input_tt, input_parameter, input_request, input_xy,input_vehicle, input_cs, instance_name)
        # best_solu, output_solver = darp_solver_batch(path, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver)  
        # penalty, penalty_per_veh_used  
        # N_ITER, step_size_remove_route, max_stagnant_multiplier 
        global_parameter.parameter_sa.penalty_per_veh_used=  val_tested
        # global_parameter.parameter_sa.step_size_remove_route=  val_tested   
        @show(val_tested)   
        best_solu, output_solver =  darp_solver_batch(path_result, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver)  
        n_route = best_solu.n_route ; veh_ids = best_solu.RI[1:n_route,5] 
        n_ev = sum(instance.parameter_energy.is_electric[veh_ids]) 
        tmp_res= [output_solver.avg_obj, output_solver.obj_best_sol, n_route-n_ev, n_ev, output_solver.co2_best_sol, output_solver.charging_time_best_sol, output_solver.t_cpu_run]
        
        # set_id, data_id, avg_obj, best_obj, co2, chg_time, n_gaso_used, n_ev_used, t_cpu_run
        res_all[count_1,1:2] = [id_data, count_2]
        res_all[count_1,3:9] = tmp_res 
        writedlm( path_result * "res_cali_n_penalty_veh_$ite"  * ".csv",  res_all, ',')
    end  
end  
 
for i in 1:length_vec_para_test
    seq_tmp = [i+ length_vec_para_test * (j-1) for j in 1:n_instance]
    res_avg[i,2:4] = [sum(res_all[seq_tmp, 3]), sum(res_all[seq_tmp, 4]), sum(res_all[seq_tmp, 9])] # avg_opt_value, best_opt_val, cpu_time
end
res_avg[:,2:4] ./= n_instance
writedlm( path_result * "res_cali_n_penalty_veh_avg_$ite"*".csv",  res_avg, ',')

end