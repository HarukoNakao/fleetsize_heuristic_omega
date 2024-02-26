################################################################
# not ready yet
# 12.2.2024 
#  
################################################################
 
include("main.jl")

########################################
# fleet_size solver: run this section (Bettembourg)
########################################
flag_dev = true # to export the detail for verification
set_seed_number = false # set fixed (or not) seed number  
#  percentage_full =  1 # 1 100% full,  percentage_full = p, 0<p<1, p% full for all EVs,  
#  percentage_full = -1, set as 30%,40%,...80% 
percentage_full = 1 # percentage full of initial battery level
step_size_remove_route = 250 # 250 for Bettembourg

instance_name   = "c_500.csv"
str_dataset     = "\\data_bettembourg" # used dataset 
path            = pwd()*"\\output\\"
path_data       = pwd()* str_dataset * "\\"

input_tt        = path_data * "Timetable_bettembourg.csv"
input_parameter = path_data * "Parameters.csv" 
input_request   = path_data * instance_name
input_xy        = path_data * "Coordinate_6h_23h.csv" # coords of bus stops, transit station and the depot is one the last line
input_vehicle   = path_data * "VehicleInfo.csv"
input_cs        = path_data * "CSInfo.csv" 

global_parameter = parameter_settings(set_seed_number, percentage_full, step_size_remove_route)  
instance, darp, scenarios_demand, nodes_info, timeTable, fast_chg, output_solver, demand_data = initialization(global_parameter, input_tt, input_parameter, input_request, input_xy,input_vehicle, input_cs, instance_name)

silent_fleet_size_solver= false; silent_darp_solver = false
global_parameter.n_run =3 # num of darp runs for each demand scenario
global_parameter.parameter_sa.step_size_remove_route = 250 #250
global_parameter.co2_reduc_target = 0.9

# t_start = time()
fleet_size_solver(path, silent_fleet_size_solver, silent_darp_solver, global_parameter, darp, fast_chg, instance, scenarios_demand, 
                  demand_data, nodes_info, timeTable, output_solver)
# t_elapse = time()- t_start  
# @show(t_elapse)

########################################################
# solve the lower level problems
########################################################
########################################################
# using gasoline vehs: run this section
########################################################  
global_parameter.n_run =3; global_parameter.co2_reduc_target = 0
global_parameter.co2_threshold= global_parameter.co2_threshold_init
silent_darp_solver = false; global_parameter.is_gasoline_fleet =true
best_solu, _ = darp_solver(path, silent_darp_solver,global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver)  

########################################################
# solve the co2 emission constraints with limited num of chargers: run this section
########################################################
# user-defined parameters
percentage_full = 0.5 # 1 100% full,  percentage_full = p, 0<p<1, p% full for all EVs,  
silent_darp_solver = false; 
global_parameter.n_run = 5
co2_reduc_target = 0.9#  reduce XX% CO2 
global_parameter.is_gasoline_fleet = false
total_chargers = [2, 0] # 2 DC chargers at location 1

global_parameter.co2_reduc_target = co2_reduc_target # % of co2 reduction, user-defined parameter 
vec_co2_max_milp = zeros(Float32,14)
if str_dataset == "data_fleet_09022023_1"
    vec_co2_max_milp[1:5]=[26.3277762667357, 39.139795180004, 58.2477553979516 , 72.9762773663471, 85.571741114915] # max co2 emission for c10 to c50 instance based on the milp solutions, used for testing the darp solver with co2 constraints
elseif str_dataset == "data_fleet_09022023_2"
    vec_co2_max_milp[1:5]=[21.893017365244, 36.558900831143, 62.943366326296, 71.973729632130, 93.974812326737] # max co2 emission for c10 to c50 instance based on the milp solutions, used for testing the darp solver with co2 constraints
else
    vec_co2_max_milp[1:5]=[22.682061480231,	41.172339316297, 53.817788060028, 77.071818292654, 93.814572611713] # max co2 emission for c10 to c50 instance based on the milp solutions, used for testing the darp solver with co2 constraints
end
@show(vec_co2_max_milp[1:5])
# below: co2 emissions for c60,70,80,90,100,200,300,400,500 of Dataset 1
vec_co2_max_milp[6:14]= [104.209, 127.889, 142.713, 156.182 , 165.025,  293.256, 449.850, 525.444, 679.207] # obtained from SA using all gv
vec_co2_threshold = vec_co2_max_milp.*(1-co2_reduc_target); # as above
setup_chargers(total_chargers, fast_chg) 
best_solu = darp_with_co2_constraint(path, silent_darp_solver, total_chargers, vec_co2_threshold, global_parameter, darp, fast_chg, instance, demand_data, nodes_info, output_solver)

# profiling to see the computational bottleneck 
# @profview  darp_with_co2_constraint(path, silent_darp_solver, total_chargers, vec_co2_threshold, global_parameter, darp, fast_chg, instance, demand_data, nodes_info, output_solver)
