
# parameters for the simulated annealing algorithm, set as mutable for the calibration and sensitivity analysis
mutable struct Parameter_SA #set as mutable for sensitivity analysis and parameter calibration
 
    N_ITER::Int32
    step_size_remove_route::Int32 #to be removed
    max_stagnant_multiplier::Int32
    T_red::Float32
    t_max::Float32
    n_imp::Int32
    penalty::Float32
    penalty_per_veh_used::Float32
end

# container to store the segment and changed route, added on 1.9.2023
mutable struct Container_route
    
    insert_cost::Float32
    regret::Vector{Float32}
    r::Int32 #route_id
    best_pos::Vector{Int32}
    new_route::Vector{Int32}
end

mutable struct Parameter_chg_infra
    
    unit_charger_cost::Float32
    vec_cs_site_cost::Vector{Float32}
    vec_max_elec_supply::Vector{Float32}
end


mutable struct Global_parameter # set as mutable for sensitivity analysis and parameter calibration
 
    is_gasoline_fleet::Bool
    n_run::Int32 # num of runs to get the average and best result of n runs  
    set_seed_number::Bool # set fixed seed number for reproduction of the same results
    
    penalty::Float32 # penalty of one unserved customer
    # instance related
    n_type_veh::Int32 
    detour_factor::Float32 
    duplcate_charger::Int32

    #algorithm/solver related
    bigM::Float32
    n_init_sol::Int32
    flag_output_to_file::Bool
    tw_width::Float32
    VoWT::Float32
    parameter_sa::Parameter_SA 
    n_scenario::Int32
    flag_initialization::Bool
    parameter_chg_infra::Parameter_chg_infra
    co2_reduc_target::Float32 # pi in eq. (9), between 0 and 1    
    co2_gasoline_fleet::Float32
    co2_threshold::Float32
    co2_threshold_init::Float32
    rho::Float32
    big_penalty_co2::Float32
    max_passenger_request::Int32
    percentage_full::Float32 
    penalty_per_veh_used::Float32
    max_num_chg::Int32 # maximum number of charging operations per vehicle, same as that in parameter_energy
    degree_rand_chg_policy::Float32 
    solve_veh_exchange_t_limit::Float32
    step_size_create_route::Int32
    # max_n_request_removed::Int32 
end


# DARP instance, unmutable object
mutable struct DARP  

    n_cus::Int32
    TH::Float32
    start_depot::Int32
    end_depot::Int32
    ei::Vector{Float32}
    li::Vector{Float32}
    s_t::Vector{Float32}
    Li::Vector{Float32}
    qi::Vector{Int32} 
    set_physical::Vector{Int32}
    is_last_mile::Vector{Int32}
    bus_stop_id::Vector{Int32}
    train_id::Vector{Int32}
    train_station_id::Vector{Int32}
    dist::Matrix{Float32}      # dist matrix with bigM (infeasible arcs) for darp instance (including the depot nodes and all pickup and dropoff nodes only)
    dist_orig::Matrix{Float32} # dist matrix without bigM for the darp instance
    dist_all::Matrix{Float32} # a copy of the distance matrix of all nodes (bus, train, chg_station, depot)
end


mutable struct Veh_info

    n_type_veh::Int32
    type_veh::Vector{Int32}
    cap_veh_type::Vector{Int32}
    cap_fuel_type::Vector{Float32}
    min_fuel_type::Vector{Float32}
    energy_consumption_km_type::Vector{Float32}
    cost_energy_type::Vector{Float32}
    daily_purchase_cost_type::Vector{Float32}
    CO2_type::Vector{Float32}
    CO2_km_type::Vector{Float32}
    fleet_size_type::Vector{Int32}
    is_electric_type::Vector{Bool} 
end

struct Parameter_energy
    pw_charger::Vector{Float32}
    veh_type::Vector{Int32} # vec of veh type for the fleet, 
    is_electric::Vector{Bool} # vec of veh type for the fleet, 1 electric, 0 disel
    beta::Vector{Float32}  # energy consumption rate of each vehicle per minute travelled
    beta_km::Vector{Float32}  # energy consumption rate of each vehicle per km travelled
    E_max::Vector{Float32} # max SOC of each vehicle
    E_min::Vector{Float32}  # min SOC of each vehicle
    E_init::Vector{Float32}  # init SOC of each vehicle
    B_cap::Vector{Float32}  # battery capacity of each vehicle
    max_ec::Vector{Float32} # max energy consumption up to SOC = E_min
    cap_passenger::Vector{Int32} # passenger capacity of each vehicle (this is not related to energy, but this information is storekd here)
    e_price::Vector{Float32}  # euro/kWh
    fleet_cost::Vector{Float32} 
    co2_emission::Vector{Float32} 
    co2_emission_km::Vector{Float32} 
    veh_info::Veh_info  
    max_num_chg::Int32 # maximum number of charging operations per vehicle 
    degree_rand_chg_policy::Float32 
    veh_ids_types::Vector{Vector{Int32}} # veh_ids of different types of vehs
end 

 
# n_c: num of customers, n_bs: num of bus stops, n_k: num of veh, n_ts: num of transit stations, Q_max: cap of veh, v_k : speed of veh
# t_s: customer service time, set as 0.5 minutes at pickup (bus stop) location, 0 at drop-off transit stations on the layered graph
mutable struct Instance # unmutable object
    T_start::Float32; n_depot::Int32; n_c::Int32; n_bs::Int32; n_ts::Int32; n_departure_train::Int32; K::Int32; K_MAX::Int32; 
    vec_n_veh::Vector{Int32}; Q_max::Vector{Int32};v_k::Float32; t_s::Float32; T_max::Float32; parameter_energy::Parameter_energy;
    instance_name::String
end

mutable struct Fast_chargers 

    pw::Vector{Float32}                     # power
    max_n_fast_chg::Int32                   # max number ofchargers can be installed for all chg sites
    idxs_rapid_chg::Vector{Int32}           # index of fast chargers, 1,2,...
    idxs_chg_station::Vector{Int32}         # idx of the charging station (1,2,3,..)
    v_physical_chg::Vector{Int32}           # physical node ids of chargers
    v_chg::Vector{Int32}                    # node ids of chargers, set as 10000+1,+2,...
    max_n_chg_by_site::Vector{Int32}          # max allowed number of chargers can be installed by site
    list_chgs_idx_by_site::Matrix{Int32}    # chg ids of the chargers at each site
    # variables below to be updated when new chargers are created
    n_fast_chg_installed::Int32             # total number of chargers installed  
    num_chg_installed_by_site::Vector{Int32}          # num of installed chargers at each site 
    set_id_chg_installed::Vector{Int32}     # set of installed chg ids 
    cost_fast_chg::Vector{Float32}
    occ_state_chg::Matrix{Int32}
end

# RI: route information
mutable struct Solution 
    succ::Vector{Int32}               # 1:2*n_cus+1 (1 is the depot)
    pre::Vector{Int32}                # 1:2*n_cus+1 (1 is the depot)
    RI::Matrix{Int32}                 # col: 1:first node, 2:last node, 3:num of nodes of the route including the pickup or dropoff nodes only, 4:capacity of veh, 5: veh_id
    unserved_users::Vector{Int32}        # pool of unserved users, users are the pickup nodes of DARP instances
    penalty_unserved::Float32
    n_route::Int32         
    total_cost::Float32               # obj. function value
    total_cost_with_penalty::Float32  # obj. function value with the penaties #added on 2.8.2023
    total_chg_time::Float32      
    vec_chg_event::Matrix{Float32} 
    co2_emission::Float32 # get the co2 emission of the solution
end

# coord of physical nodes etc
struct Nodes_info
    n_bs::Int32
    n_ts::Int32
    n_chg_station::Int32
    coord_bs::Matrix{Float32} 
    coord_transit::Matrix{Float32} 
    coord_chg_station::Matrix{Float32}
    coord_depot::Matrix{Float32} 
    coord_all::Matrix{Float32}  
    v_ts_physical::Vector{Int32}
    v_cs_all_physical::Vector{Int32}
    v_depot_physical::Int32 
    dist_all::Matrix{Float32}  #dist matrix of all physical nodes (bus, transit, charging station and the depot)
end

mutable struct Demand_data

    n_c::Int32 # num of customers
    n_station_used::Int32
    id_request::Vector{Int32}
    id_bus_stop::Vector{Int32}	
    no_person::Vector{Int32}
    id_train::Vector{Int32}
    first_last_mile::Vector{Bool} #0 if customer uses first mile service/ 1 if last mile
    id_train_station::Vector{Int32} 
end

mutable struct Output_solver

    n_run::Int32
    avg_obj::Float32
    obj_values::Vector{Float32}
    avg_co2::Float32
    std_obj::Float32
    std_co2::Float32
    cv_obj::Float32 #coefficient of variation : std/mean
    cv_co2::Float32
    obj_best_sol::Float32 
    co2_best_sol::Float32  
    charging_time_best_sol::Float32 
    t_cpu_run::Float32 
end

mutable struct Fleet_size_solution

    n_request::Int32 
    n_scenario::Int32
    co2_reduc_target::Float32
    obj_l1::Float32
    obj_l2::Float32
    fleetsize_gv::Int32
    fleetsize_ev::Int32
    cs_config::Tuple{Int32, Int32} # number of rapid chargers installed on each cs site
end

mutable struct Relatedness

    dist::Matrix{Float32}
    tw::Matrix{Float32}
    shaw::Matrix{Float32}

end

# # layered graph, unmutable object
# struct LGraph  

#     n_c_lyr::Int32
#     vec_used_lyr::Vector{Int32}
#     nodes_bus::Vector{Int32}
#     G′::Matrix{Int32}
#     D′::Matrix{Int32}
#     G::Vector{Int32}
#     D::Vector{Int32}
#     S::Vector{Int32}
#     v_depot::Int32
#     S′::Vector{Int32}
#     set_d_l::Vector{Int32}
#     set_s_l::Vector{Int32}
#     set_l_r::Vector{Int32}     # customers' original layer indexes
#     set_l_r_new::Vector{Int32} # customers' new layer indexes by deleting unused layers
#     set_d_r::Vector{Int32}
#     set_d_r_physical::Vector{Int32}
#     set_physical::Vector{Int32}
#     lyr_bus::Vector{Int32}
#     lyrs_nodes::Vector{Int32}
#     lyrs_compatible::Matrix{Bool}
#     lyrs_compatible_extend::Matrix{Bool}
#     dist_GD::Matrix{Float32}    #travel time matrix for the nodes of G and D on the same layer 
#     dist_all::Matrix{Float32}   #travel time matrix for the nodes of all physical nodes (including the depot and chargers)
#     coord_V::Matrix{Float32}
# end