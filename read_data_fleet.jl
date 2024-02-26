# adapted from Haruko's code
function initialization(global_parameter::Global_parameter, input_tt::String, input_parameter::String, input_request::String,
                        input_xy::String, input_vehicle::String, input_cs::String, instance_name::String)
    
    n_scenario =  global_parameter.n_scenario 
    max_num_chg = global_parameter.max_num_chg
    timeTable = DataFrame(CSV.File(input_tt))
    data_file = DataFrame(CSV.File(input_parameter)) # parameter
    coordinates = DataFrame(CSV.File(input_xy))
    vehicleInfo = DataFrame(CSV.File(input_vehicle))
    chargerInfo = DataFrame(CSV.File(input_cs))
    requestInfo = DataFrame(CSV.File(input_request))

    K_MAX = max( sum(vehicleInfo[:, "fleet_size"]), 100) # make it bigger as we now have c500 instances
 
    T_start =0; n_depot=1
    T_max = convert(Float32, timeTable.arrival_time[end]) + 40 # time horizon7
    n_c = size(requestInfo, 1) #number of trip requests
    
    n_bs = data_file[1, "no_bus_stop"] #number of bus stops
    n_ts = data_file[1, "no_train_station"] #number of transit stations
    n_type_veh = data_file[1, "no_bus_type"] #number of vehiclce type
    t_s = data_file[1, "service_time"] #service time(getting on/off the bus)
    v_k = data_file[1, "vehicle_speed"] # bus speed
    n_chg_station = data_file[1, "n_char_station"] 
    n_type_veh ==1 && error("there is only one type of vehicle error !",n_type_veh)

    cs_site_cost_daily = chargerInfo[:, "cs_site_cost_daily"]  
    max_elec_supply = chargerInfo[:, "max_elec_supply"]
    cost_rapid_chg	= chargerInfo[:, "cost_rapid_chg"] 

    # setup 
    global_parameter.parameter_chg_infra.vec_cs_site_cost = cs_site_cost_daily
    global_parameter.parameter_chg_infra.vec_max_elec_supply = max_elec_supply
    global_parameter.n_type_veh = n_type_veh
    # @show(global_parameter.n_type_veh)
    percentage_full = global_parameter.percentage_full
    """
    Vehicle info, set up the fleet
    """
    vec_n_veh = Int32[]
    cap_veh = zeros(Int32, n_type_veh)
    B = zeros(Float32, K_MAX)
    E_max = zeros(Float32, K_MAX)
    E_min = zeros(Float32, K_MAX)
    E_init = zeros(Float32, K_MAX)
    beta = zeros(Float32, K_MAX)
    beta_km = zeros(Float32, K_MAX)
    max_ec = zeros(Float32, K_MAX)
    fleet_cost = zeros(Float32, K_MAX)
    e_price = zeros(Float32, K_MAX)
    co2_emission = zeros(Float32, K_MAX)    
    co2_emission_km = zeros(Float32, K_MAX)
    veh_type= zeros(Int32,K_MAX) 
    is_electric= zeros(Bool,K_MAX) # 0 disel, 1 electric
    beta_km_tmp = zeros(Float32, n_type_veh)
    B_tmp = zeros(Float32, n_type_veh)
    e_max = zeros(Float32, n_type_veh)
    e_min = zeros(Float32, n_type_veh)
    veh_type_tmp = zeros(Int32, n_type_veh)
    is_electric_tmp = zeros(Bool, n_type_veh)
    e_price_type = zeros(Float32, n_type_veh)
    fleet_cost_type = zeros(Float32, n_type_veh)
    co2_emission_type = zeros(Float32, n_type_veh)
    co2_emission_type_km = zeros(Float32, n_type_veh)
    ##### Vehicle type 1: Gasoline vehicle
    for i in 1:n_type_veh
        n_k_i = vehicleInfo[i, "fleet_size"]#fleet size of vehilce type 1
        push!(vec_n_veh, n_k_i) 
        cap_veh[i] = vehicleInfo[i, "cap_veh"]#vehicle capactiy 
        beta_km_tmp[i] = vehicleInfo[i, "energy_consumption"] #energy consumption rate for bus type 1
        e_max[i] = vehicleInfo[i, "cap_fuel"] #the max charging state. Big enough nubmer for gasoline vehilce 
        e_min[i] = vehicleInfo[i, "min_fuel"] #the min charging state. 0 for gasoline vehicle 
        B_tmp[i] = e_max[i]  #the initial buttery status 
        veh_type_tmp[i] = vehicleInfo[i, "vehicle_type"] 
        is_electric_tmp[i] = vehicleInfo[i, "is_electric"] 
        ### new parameters ###
        e_price_type[i] = vehicleInfo[i, "cost_energy"] #energy price for vehicle type 1 
        fleet_cost_type[i] = vehicleInfo[i, "daily_purchase_cost"] #purchasing and maintenace cost of vhicle type 1 /day 
        co2_emission_type_km[i] = vehicleInfo[i, "CO2"] #CO2 emission rate 
        co2_emission_type[i] =  co2_emission_type_km[i] * v_k # CO2 per miniute travelled
    end
    veh_info =  Veh_info(n_type_veh, veh_type_tmp, cap_veh, e_max, e_min, beta_km_tmp, e_price_type, fleet_cost_type, 
                         co2_emission_type, co2_emission_type_km, vec_n_veh, is_electric_tmp)
    #  @show(n_type_veh, veh_type_tmp, cap_veh, e_max, e_min, beta_km_tmp, co2_emission_type, co2_emission_type_km, vec_n_veh, is_electric_tmp)
    #total fleet size 
    n_k = sum(vec_n_veh)
    # @show(n_k, vec_n_veh)
    if set_seed_number 
        Random.seed!(1000)
    end
    
    ## decide the initial state of each vehicle 
    max_cap_veh = maximum(cap_veh)
    Q_max = ones(Int32, K_MAX) .* max_cap_veh
    veh_id = 0
    veh_ids_types= Vector{Vector{Int32}}(undef, n_type_veh)
    for i in 1:n_type_veh
        veh_ids_types[i]= Int32[]
        for j in 1:vec_n_veh[i]
            veh_id += 1
            Q_max[veh_id] = cap_veh[i]
            B[veh_id] = B_tmp[i]
            E_max[veh_id] = e_max[i]
            E_min[veh_id] = e_min[i]  
            E_init[veh_id]  = e_max[i]  # for all veh
        
            if is_electric_tmp[i] == true # for ev only
                if  percentage_full >0 
                    E_init[veh_id] = e_max[i] * percentage_full
                else # -1 # 30%,40%,...
                    E_init[veh_id] = e_max[i] * min(0.8, 0.3 + 0.1 * (j%6) )  # not used
                end
            end 
            beta_km[veh_id] = beta_km_tmp[i]
            beta[veh_id] = beta_km[veh_id] * v_k
            veh_type[veh_id] = veh_type_tmp[i]                
            is_electric[veh_id] = is_electric_tmp[i]
            max_ec[veh_id] = E_init[veh_id] - E_min[veh_id]
            e_price[veh_id]      = e_price_type[i]  
            fleet_cost[veh_id]   = fleet_cost_type[i] 
            co2_emission[veh_id] = co2_emission_type[i]   
            co2_emission_km[veh_id] = co2_emission_type_km[i]
            push!(veh_ids_types[i], veh_id)  
        end
    end  
 
    n_veh_added = max(0, K_MAX - n_k)  
    vec_n_charger = chargerInfo[:, "max_no_chargers"]# max num of chargers can be installed at the candidate cs sites (we assume that are the transit stations)
    max_n_charger = floor(Int32, sum(vec_n_charger))
    alpha_DC = chargerInfo[:, "alpha_DC"] # kw/min charging speed of DC chargers
    alpha = zeros(Float32, max_n_charger); id_chg = 0
    for i in 1:size(vec_n_charger)[1]
        for j in 1:vec_n_charger[i]
            id_chg += 1
            alpha[id_chg] = alpha_DC[i]
        end
    end
    pw_tmp = alpha
    # create additional dummy vehicles equally for all types of vehs
    seg=zeros(Int32,n_type_veh)
    sep_points = zeros(Int32,n_type_veh)
 
    if n_veh_added > 0 
        if global_parameter.is_gasoline_fleet
            n_gaso_veh = vec_n_veh[1]        
            veh_type[n_gaso_veh+1:end]    .= veh_type[1]
            is_electric[n_gaso_veh+1:end] .= is_electric[1]
            beta_km[n_gaso_veh+1:end] .= beta_km[1]
            beta[n_gaso_veh+1:end] .= beta[1]
            E_max[n_gaso_veh+1:end] .= E_max[1]
            E_min[n_gaso_veh+1:end] .= E_min[1]
            E_init[n_gaso_veh+1:end] .=E_init[1]
            B[n_gaso_veh+1:end] .= B[1]
            max_ec[n_gaso_veh+1:end] .= B[1]
            Q_max[n_gaso_veh+1:end] .= Q_max[1]
            e_price[n_gaso_veh+1:end] .= e_price[1]
            fleet_cost[n_gaso_veh+1:end] .= fleet_cost[1]
            co2_emission[n_gaso_veh+1:end] .= co2_emission[1]
            co2_emission_km[n_gaso_veh+1:end] .= co2_emission_km[1]
            veh_ids_types[1] = collect(1:K_MAX)
            for j in 2:n_type_veh
                veh_ids_types[j] = Int32[]
            end
        else
            ratio = 1.0/n_type_veh
            sep_points_tmp = [floor(Int32, n_veh_added * ratio* i) for i in 1:n_type_veh-1]
            sep_points[1], sep_points[2:end] = n_k,  n_k .+ sep_points_tmp 
            segs =[ [sep_points[i]+1, sep_points[i+1]]   for i in 1:n_type_veh-1]
            push!(segs, [sep_points[end]+1, n_k+n_veh_added])
            
            for tt in 1:n_type_veh # three types including gasoline
                veh_type[segs[tt][1]: segs[tt][2]] .= veh_type_tmp[tt];  is_electric[segs[tt][1]: segs[tt][2]] .= is_electric_tmp[tt];  
                beta_km[segs[tt][1]: segs[tt][2]] .= beta_km_tmp[tt]; beta[segs[tt][1]: segs[tt][2]] = beta_km[segs[tt][1]: segs[tt][2]] .* v_k
                E_max[segs[tt][1]: segs[tt][2]] .= e_max[tt];  E_min[segs[tt][1]: segs[tt][2]] .= e_min[tt];  
                E_init[segs[tt][1]: segs[tt][2]] .= e_max[tt] * percentage_full; 
                B[segs[tt][1]: segs[tt][2]] .= B_tmp[tt]; max_ec[segs[tt][1]: segs[tt][2]] .= (e_max[tt] * percentage_full -  e_min[tt])
                Q_max[segs[tt][1]: segs[tt][2]] .= cap_veh[tt]; 
                e_price[segs[tt][1]: segs[tt][2]] .= e_price_type[tt];  fleet_cost[segs[tt][1]: segs[tt][2]] .= fleet_cost_type[tt]; 
                co2_emission[segs[tt][1]: segs[tt][2]] .= co2_emission_type[tt]; 
                co2_emission_km[segs[tt][1]: segs[tt][2]] .= co2_emission_type_km[tt]; 
                veh_ids_types[tt] = vcat(veh_ids_types[tt],collect(segs[tt][1]: segs[tt][2]))  
            end
        end
    end
    # check veh_ids_types
    vec_count_veh_ids_type = zeros(Int,K_MAX)
    for j in 1:n_type_veh
        for k in veh_ids_types[j]
            vec_count_veh_ids_type[k] += 1
        end
    end
    flag =  findfirst(x->x!=1, vec_count_veh_ids_type)
    if !isnothing(flag)
        @show(vec_count_veh_ids_type, veh_ids_types)
        error("veh_ids_types error!")
    end
    # @show(K_MAX, veh_type, veh_ids_types, is_electric, E_max, E_init, E_min, max_ec, beta, Q_max, B, e_price, co2_emission, co2_emission_km )
    max_ec = E_init - E_min 
    degree_rand_chg_policy = global_parameter.degree_rand_chg_policy

 
    # beta has additional values as fictitious for additional vehs more than n_k 
    parameter_energy = Parameter_energy(pw_tmp, veh_type, is_electric, beta, beta_km,
                        E_max, E_min, E_init, B, max_ec, Q_max, e_price, 
                        fleet_cost, co2_emission, co2_emission_km, veh_info, max_num_chg, degree_rand_chg_policy, veh_ids_types) # used for E-DARP
    # global_parameter.threshold_advantage_ev= get_threshold_advantage_ev(parameter_energy, v_k) 
 
    # n_c is the num of requests
    coord_bs =zeros(Float32, n_bs,2)
    coord_bs[:,1],coord_bs[:,2] =  coordinates.x[1:n_bs] ,coordinates.y[1:n_bs] 
    coord_transit = zeros(Float32, n_ts, 2)
    for i in 1:n_ts
        coord_transit[i,:]  = Float32.([coordinates.x[n_bs+i],coordinates.y[n_bs+i]] )
    end
    coord_depot = zeros(Float32, 1,2)
    coord_depot[1,:] = Float32.([coordinates.x[end] ,coordinates.y[end]] )
    # n_chg_station = length(chargerInfo[:, "id_char_station"]) #there are several charging station sites 
    coord_chg = zeros(Float32, max_n_charger,2)
    v_chg_physical = zeros(Int32, max_n_charger) #physical node id of each charger
    coord_chg_station = zeros(Float32, n_chg_station,2)
    coord_chg_station[:,1] = Float32.(chargerInfo.x_char_station )
    coord_chg_station[:,2] = Float32.(chargerInfo.y_char_station)
    max_n_chg_by_site = 30 # a big number never exceeds
    count_chg=0; idxs_chg_staion = zeros(Int32, max_n_charger)
    list_chgs_idx_by_site=zeros(Int32,n_chg_station, max_n_chg_by_site)
    n_chg_by_site=zeros(Int32,n_chg_station)
    for i in 1:n_chg_station
        for _ in 1:vec_n_charger[i]
            count_chg += 1 
            coord_chg[count_chg,:]= coord_chg_station[i,:]
            idxs_chg_staion[count_chg]= i
            n_chg_by_site[i] += 1
            list_chgs_idx_by_site[i, n_chg_by_site[i]] = count_chg
        end
    end
    # initialization
    coord_all= vcat(coord_bs, coord_transit, coord_chg_station, coord_depot)   
    v_ts_physical    = collect(n_bs+1: n_bs + n_ts)
    # vector of physical nodes of the chargers
    v_cs_all_physical = collect(v_ts_physical[end]+1: v_ts_physical[end]+n_chg_station) # physical node id of all charging stations
    for s in 1:max_n_charger
        v_chg_physical[s] = v_ts_physical[end] + idxs_chg_staion[s] # physical node id of each charger
    end
    v_depot_physical = v_cs_all_physical[end]+1 
    # dist_all: travel time matrix for all bus
    dist_all = comp_dist(coord_all,v_k)
    dist_all = dist_all ./1000     
    nodes_info= Nodes_info(n_bs, n_ts, n_chg_station, coord_bs, coord_transit, coord_chg_station, coord_depot, coord_all,
                           v_ts_physical, v_cs_all_physical, v_depot_physical, dist_all);
    max_n_chg_by_site = vec_n_charger 
    fast_chg = init_chargers(parameter_energy, n_chg_station, max_n_charger, v_chg_physical, idxs_chg_staion,
                              max_n_chg_by_site, list_chgs_idx_by_site, cost_rapid_chg) ; # fast_chg are all the chargers
    for cs in 1:n_chg_station
        for _ in 1:vec_n_charger[cs]
            add_charger(fast_chg, cs) #DC fast 
        end 
    end 
    discrte_t_interval = 1/6 # 10 second 
    n_t_max = floor(Int32, (T_max+20)/discrte_t_interval) # 840 discrete charging occupancy state over T_max+20 minutes
    n_charger = fast_chg.n_fast_chg_installed
    # @show(fast_chg, n_charger )
    occ_state_chg = zeros(Int32, n_charger, n_t_max) 
    fast_chg.occ_state_chg = occ_state_chg
    scenarios_demand = Vector{Demand_data}(undef,n_scenario)
    demand_data = init_demand_scenario(n_c, n_ts, input_request, nodes_info, scenarios_demand)
    total_node_darp = 2*n_c+2
    # size(dist_all)
    # default is one station but decide how many chargers to install 
    
    # @show(parameter_energy)
    n_departure_train = size(timeTable)[1]
    instance = Instance(T_start, n_depot, n_c, n_bs, n_ts, n_departure_train, n_k, K_MAX, vec_n_veh, Q_max, v_k, t_s, T_max, parameter_energy, instance_name);
    obj_vals = zeros(Float32,10)
    output_solver = Output_solver(0,0,obj_vals,0,0,0,0,0,0,0,0,0);
    n_V= total_node_darp-1+max_n_charger 
    # output the darp instance for the first scenario for verification
    darp, data_darp = generate_darp_instance(demand_data, instance, global_parameter, nodes_info, timeTable, fast_chg, 0);
    if flag_dev  
    #    darp, data_darp = generate_darp_instance(demand_data, instance, global_parameter, nodes_info, timeTable, fast_chg, 0);
       generate_V_DF(data_darp, darp, nodes_info, instance, demand_data, n_c, n_V, fast_chg)
    end
    return instance, darp, scenarios_demand, nodes_info, timeTable, fast_chg, output_solver, demand_data
end

# adapted from Haruko's code
# function load_data_fleet_size(global_parameter::Global_parameter, input_tt::String, input_parameter::String, input_request::String, input_xy::String, input_vehicle::String, input_cs::String, is_gasoline_fleet)
    
#     tw_width, detour_factor, bigM = global_parameter.tw_width, global_parameter.detour_factor, global_parameter.bigM
#     timeTable = DataFrame(CSV.File(input_tt))
#     data_file = DataFrame(CSV.File(input_parameter)) # parameter
#     requestInfo = DataFrame(CSV.File(input_request))
#     coordinates = DataFrame(CSV.File(input_xy))
#     vehicleInfo = DataFrame(CSV.File(input_vehicle))
#     chargerInfo = DataFrame(CSV.File(input_cs))
 
#     K_MAX = max(vehicleInfo[1, "fleet_size"] + vehicleInfo[2, "fleet_size"], 50)
#     T_start =0; n_depot=1
#     T_max = convert(Float32, timeTable.arrival_time[end]) + 15 # time horizon7
#     #data_file = open(input_data)
#     n_c = nrow(requestInfo) #number of trip requests
#     n_bs = data_file[1, "no_bus_stop"] #number of bus stops
#     n_ts = data_file[1, "no_train_station"] #number of transit stations
#     n_type_veh = data_file[1, "no_bus_type"] #number of vehiclce type
#     t_s = data_file[1, "service_time"] #service time(getting on/off the bus)
#     v_k = data_file[1, "vehicle_speed"] # bus speed
#     n_chg_station = data_file[1, "n_char_station"] 

#     if is_gasoline_fleet
#        n_type_veh =1
#     end
 
#     n_tmp =  length(findall(x->x==0, chargerInfo.No_chargers))
#     n_chg_station_used = length(chargerInfo.No_chargers) - n_tmp
#     # if n_chg_station !=  n_chg_station_used
#     #     @show("warning!! n_chg_station !=  n_chg_station_used ! Please check data files! ",n_chg_station, n_chg_station_used)
#     # end
#     # n_tmp= length(findall(x->x==0, vehicleInfo.fleet_size))
#     # n_type_veh_used = length(vehicleInfo.fleet_size) - n_tmp
#     # if n_type_veh_used != n_type_veh
#     #     error("n_type_veh_used!=n_type_veh error! Please check data files! ",n_type_veh_used, n_type_veh)
#     # end

#     """
#     Vehicle info, set up the fleet
#     """
#     vec_n_veh = Int32[]
#     cap_veh = zeros(Int32, n_type_veh)
#     B = zeros(Float32, K_MAX)
#     E_max = zeros(Float32, K_MAX)
#     E_min = zeros(Float32, K_MAX)
#     E_init = zeros(Float32, K_MAX)
#     beta = zeros(Float32, K_MAX)
#     max_ec = zeros(Float32, K_MAX)
#     fleet_cost = zeros(Float32, K_MAX)
#     e_price = zeros(Float32, K_MAX)
#     co2_emission = zeros(Float32, K_MAX)
#     veh_type= zeros(Int32,K_MAX) 
#     is_electric= zeros(Bool,K_MAX) # 0 disel, 1 electric
#     beta_tmp = zeros(Float32, n_type_veh)
#     B_tmp = zeros(Float32, n_type_veh)
#     e_max = zeros(Float32, n_type_veh)
#     e_min = zeros(Float32, n_type_veh)
#     veh_type_tmp = zeros(Int32, n_type_veh)
#     is_electric_tmp = zeros(Bool, n_type_veh)
#     e_price_type = zeros(Float32, n_type_veh)
#     fleet_cost_type = zeros(Float32, n_type_veh)
#     co2_emission_type = zeros(Float32, n_type_veh)
    
#     ##### Vehicle type 1: Gasoline vehicle
#     for i in 1:n_type_veh
#         n_k_i = vehicleInfo[i, "fleet_size"]#fleet size of vehilce type 1
#         push!(vec_n_veh, n_k_i) 
#         cap_veh[i] = vehicleInfo[i, "cap_veh"]#vehicle capactiy 
#         beta_tmp[i] = vehicleInfo[i, "energy_consumption"] #energy consumption rate for bus type 1
#         e_max[i] = vehicleInfo[i, "cap_fuel"] #the max charging state. Big enough nubmer for gasoline vehilce 
#         e_min[i] = vehicleInfo[i, "min_fuel"] #the min charging state. 0 for gasoline vehicle 
#         B_tmp[i] = e_max[1]  #the initial buttery status 
#         veh_type_tmp[i] = vehicleInfo[i, "vehicle_type"] 
#         is_electric_tmp[i] = vehicleInfo[i, "is_electric"] 
#         ### new parameters ###
#         e_price_type[i] = vehicleInfo[i, "cost_energy"] #energy price for vehicle type 1 
#         fleet_cost_type[i] = vehicleInfo[i, "daily_purchase_cost"] #purchasing and maintenace cost of vhicle type 1 /day 
#         co2_emission_type[i] = vehicleInfo[i, "CO2"] #CO2 emission rate 
#     end
#     veh_info =  Veh_info(n_type_veh, veh_type_tmp, cap_veh, e_max, e_min, beta_tmp, e_price_type, fleet_cost_type, co2_emission_type, vec_n_veh, is_electric_tmp)
    
#     #total fleet size 
#     n_k = sum(vec_n_veh)
    
#     if set_seed_number
#         seed_num = n_c * 100 + n_bs
#         Random.seed!(seed_num)
#     end
    
#     ## decide the initial state of each vehicle 
#     max_cap_veh = maximum(cap_veh)
#     Q_max = ones(Int32, K_MAX) .* max_cap_veh
    
#     if is_gasoline_fleet
#         for veh_id in 1:n_k
#             Q_max[veh_id] = cap_veh[1]
#             B[veh_id] = B_tmp[1]
#             E_max[veh_id] = e_max[1]
#             E_min[veh_id] = e_min[1]
#             E_init[veh_id] = e_max[1]
#             beta[veh_id] = beta_tmp[1]
#             max_ec[veh_id] = E_init[veh_id] - E_min[veh_id]
#             veh_type[veh_id] = veh_type_tmp[1]
#             is_electric[veh_id] = is_electric_tmp[1]
#             e_price[veh_id]      = e_price_type[1]  
#             fleet_cost[veh_id]   = fleet_cost_type[1] 
#             co2_emission[veh_id] = co2_emission_type[1]    
#         end
#     else
#         veh_id = 0
#         for i in 1:n_type_veh
#             for j in 1:vec_n_veh[i]
#                 veh_id += 1
#                 Q_max[veh_id] = cap_veh[i]
#                 B[veh_id] = B_tmp[i]
#                 E_max[veh_id] = e_max[i]
#                 E_min[veh_id] = e_min[i]
#                 E_init[veh_id] = e_max[i]
#                 beta[veh_id] = beta_tmp[i]
#                 veh_type[veh_id] = veh_type_tmp[i]                
#                 is_electric[veh_id] = is_electric_tmp[i]
#                 max_ec[veh_id] = E_init[veh_id] - E_min[veh_id]
#                 e_price[veh_id]      = e_price_type[i]  
#                 fleet_cost[veh_id]   = fleet_cost_type[i] 
#                 co2_emission[veh_id] = co2_emission_type[i]   
#             end
#         end
#     end    
   
#     """
#     Charging station info, set up chargers
#     """
#     n_chg_station = length(chargerInfo[:, "id_char_station"])
#     n_charger = sum(chargerInfo[:, "No_chargers"])# number of chargers 
#     alpha_DC = data_file[1, "alpha_DC"] # kw/min charging speed of DC chargers
#     alpha = repeat([alpha_DC], n_charger)
#     pw_tmp = alpha
 
#     n_veh_added = max(0, K_MAX - n_k)
#     # set veh > n_k as gasoline vehs
#     is_electric[1] != 0 && error("input data error: please set the first veh type as gasoline veh!!")

#     veh_type[n_k+1:n_k+n_veh_added] .= veh_type[1];  is_electric[n_k+1:n_k+n_veh_added] .= is_electric[1];  beta[n_k+1:n_k+n_veh_added] .= beta[1]; 
#     E_max[n_k+1:n_k+n_veh_added] .= E_max[1];  E_min[n_k+1:n_k+n_veh_added] .= E_min[1];  E_init[n_k+1:n_k+n_veh_added] .= E_init[1]; 
#     B[n_k+1:n_k+n_veh_added] .= B[1];  max_ec[n_k+1:n_k+n_veh_added] .= max_ec[1];  Q_max[n_k+1:n_k+20] .= Q_max[1]; 
#     e_price[n_k+1:n_k+n_veh_added] .= e_price[1];  fleet_cost[n_k+1:n_k+n_veh_added] .= fleet_cost[1];  co2_emission[n_k+1:n_k+n_veh_added] .= co2_emission[1]; 
#     # beta has additional values as fictitious for additional vehs more than n_k 
#     parameter_energy = Parameter_energy(pw_tmp, veh_type, is_electric, beta, beta.*v_k,
#                         E_max, E_min, E_init, B, max_ec, Q_max, e_price, 
#                         e_price.*v_k, fleet_cost, co2_emission, co2_emission.*v_k, veh_info, 
#                         is_gasoline_fleet, 0) # used for E-DARP
#     parameter_energy.threshold_advantage_ev= get_threshold_advantage_ev(parameter_energy, v_k, is_gasoline_fleet)
#     # @show(parameter_energy)
#     """
#     create V_DF
#     """
#     #n_c is the num of requests
#     coord_bs =zeros(Float32, n_bs,2)
#     coord_bs[:,1],coord_bs[:,2] =  coordinates.x[1:n_bs] ,coordinates.y[1:n_bs] 
#     coord_depot = zeros(Float32, 1,2)
#     coord_depot[1,:] = Float32.([coordinates.x[end] ,coordinates.y[end]] )
#     coord_transit = zeros(Float32, n_ts, 2)
#     for i in 1:n_ts
#         coord_transit[i,:]  = Float32.([coordinates.x[n_bs+i],coordinates.y[n_bs+i]] )
#     end
#     # n_ts>1 && error("warning ! n_ts > 1 !")

#     coord_chg_station = zeros(Float32, n_chg_station,2)
#     coord_chg_station[:,1] = Float32.(chargerInfo.x_char_station)
#     coord_chg_station[:,2] = Float32.(chargerInfo.y_char_station)
#     set_node_requests =  collect(Set(requestInfo.id_bus_stop))
#     set_node_requests_sorted = sort(set_node_requests)
#     n_active_bs= length(set_node_requests_sorted)
#     node_bs= zeros(Int32,n_bs) # bs node index in the node physical set
#     node_bs[set_node_requests_sorted] = collect(1:n_active_bs) 
#     v_bs_physical = copy(node_bs) # get the active bs node idx in the the nodes_physical set
#     coord_physical_nodes= vcat(coord_bs[set_node_requests_sorted,:], coord_transit, coord_chg_station,coord_depot)
#     # set_physical: node_active_bs(1:n_c), node_transit_station, node_chg_station, node_depot
#     v_ts_physical          = collect(n_active_bs+1: n_active_bs + n_ts)
#     v_chg_station_physical = collect(v_ts_physical[end]+1: v_ts_physical[end]+n_chg_station)
#     v_depot_physical       = v_chg_station_physical[end]+1
#     # dist_all: travel time matrix for the set of nodes of the active physical bus nodes and all transit, charging station 
#     #           and the depot nodes(including bus stops, transit stops, charing stations, and the depot)
#     dist_all = comp_dist(coord_physical_nodes,v_k)
#     dist_all = dist_all ./1000
#     coord_chg = zeros(Float32, n_charger,2) # customers(bus stops)' y coordinate
#     count=0;  set_physical_charger= zeros(Int32, n_charger)
#     for i in 1:n_chg_station
#         for j in 1:chargerInfo[i, "No_chargers"]
#             count +=1
#             set_physical_charger[count]= i
#             coord_chg[count,:] = [chargerInfo[i, "x_char_station"], chargerInfo[i, "y_char_station"] ]#Charging station coordinate 
#         end
#     end
#     @show(n_c, n_bs, n_ts, n_type_veh,t_s, v_k,n_chg_station, n_chg_station_used, n_charger, parameter_energy.is_gasoline_fleet )
#     # @show(parameter_energy)
#     #0 if customer uses first mile service/ 1 if last mile
#     is_last_mile = floor.(Int32, requestInfo[:, "first_last_mile"]) # 0 first mile and 1 last mile
#     # create all nodes V_DF (DF means dataframe)
#     total_node_darp = 2 * n_c + 2 # total node of a DARP instance when you have _N requests (the depot nodes: 1, 2n+2, pickup: 2:n+1, dropoff: n+2: 2n+1)
#     n_col_darp = 13 # num of columns in data_darp 
#     data_darp=zeros(Float32, total_node_darp, n_col_darp)
#     # see Cordeau (2003) A Branch-and-Cut Algorithm for the Dial-a-Ride Problem
#     # Time window tightening
#     for idx in 1:n_c # idx is idx in DARP instance, pickup: 2:n+1, drop-off: n+2:2n+1, depot: 1,2n+2
#         qi= requestInfo.No_person[idx];bus_stop_id = requestInfo.id_bus_stop[idx]
#         if is_last_mile[idx] == 0 #first mile
#             p_physical=  v_bs_physical[requestInfo.id_bus_stop[idx]]
#             d_physical=  v_ts_physical[requestInfo.id_train_station[idx]]
#             Li_tmp = dist_all[p_physical, d_physical] * detour_factor
#             # Li_tmp =10
#             li_d= timeTable.arrival_time[requestInfo.id_train[idx]]
#             ei_d=li_d-tw_width
#             ei_p = max(0, ei_d - Li_tmp -t_s) 
#             li_p = min(li_d - dist_all[p_physical, d_physical] - t_s, T_max)  
#             bus_stop_p, bus_stop_d = bus_stop_id, -1
#         else
#             p_physical=  v_ts_physical[requestInfo.id_train_station[idx]]
#             d_physical=  v_bs_physical[requestInfo.id_bus_stop[idx]] 
#             Li_tmp = dist_all[p_physical, d_physical] * detour_factor
#             # Li_tmp =10
#             ei_p =  timeTable.arrival_time[requestInfo.id_train[idx]]
#             li_p =  ei_p + tw_width
#             ei_d =  max(0, ei_p + t_s + dist_all[p_physical, d_physical])
#             li_d =  min(li_p + t_s + Li_tmp, T_max)   
#             bus_stop_p, bus_stop_d = -1, bus_stop_id
#         end
#         train_id = requestInfo.id_train[idx]
#         train_station_id = timeTable.id_train_station[train_id]
#         # the first row is the depot
#         # node_id_darp, x, y, demand,ei,li,service time, max_ridetime, physical_bus_node,     train_id, bus_node_id
#         data_darp[idx+1,:]=   [idx+1,      coord_physical_nodes[p_physical,1],coord_physical_nodes[p_physical,2], qi,  ei_p, li_p, t_s, Li_tmp, p_physical, train_id, bus_stop_p, is_last_mile[idx],train_station_id]
#         # create the corresponding drop-off node
#         data_darp[idx+1+n_c,:]=[n_c+idx+1, coord_physical_nodes[d_physical,1],coord_physical_nodes[d_physical,2], -qi, ei_d, li_d, 0, Li_tmp, d_physical, train_id, bus_stop_d, is_last_mile[idx],train_station_id]
#     end
#     data_darp[1,:]   = [1,               coord_physical_nodes[v_depot_physical,1], coord_physical_nodes[v_depot_physical,2], 0, T_start, T_max, 0,0, v_depot_physical, -1,-1,-1,-1] 
#     data_darp[end,:] = [total_node_darp, coord_physical_nodes[v_depot_physical,1], coord_physical_nodes[v_depot_physical,2], 0, T_start, T_max, 0,0, v_depot_physical, -1,-1,-1,-1] 

#     K, n_cus, TH = n_k, n_c, T_max # n_cus here is the total_node_darp
#     start_depot, end_depot = 1, total_node_darp
#     # nodes = floor.(Int32, data_darp[:,1]) # lookup table to find the node id in V.xls, we have three types of node ids (node id for DARP, node id in V.xls, node id used for comparing with MILP results (not used))
#     qi, ei, li, s_t, Li, set_physical, is_last_mile_darp= floor.(Int32, data_darp[:,4]),  data_darp[:,5], data_darp[:,6], data_darp[:,7], data_darp[:,8], floor.(Int32,data_darp[:,9]), floor.(Int32,data_darp[:,12])
#     bus_stop_id, train_id, train_station_id = floor.(Int32, data_darp[:,11]), floor.(Int32, data_darp[:,10]),floor.(Int32, data_darp[:,end])
    
#     dist, dist_orig = set_dist( ei, li, s_t, Li, set_physical, dist_all, total_node_darp, bigM)
#     # preprocess 
#     preprocess(n_cus, li, dist, s_t)     # TW has been tightened already, used just for tighten li
#     darp=  DARP(n_cus, TH, start_depot, end_depot, ei, li, s_t, Li, qi, set_physical, is_last_mile_darp, bus_stop_id, train_id, train_station_id, dist, dist_orig,dist_all)
#     instance = Instance(T_start, n_depot, n_c, n_bs, n_ts, n_k, K_MAX, vec_n_veh, Q_max, v_k, t_s, T_max, n_charger, parameter_energy,
#         coord_bs, coord_transit, coord_depot, coord_chg)
#     n_V= total_node_darp-1+n_charger
   
#     V_DF = DataFrame(id=collect(Int32, 1:n_V), v_type=Vector{String}(undef, n_V), x=zeros(Float32,n_V), y=zeros(Float32,n_V),
#            qi = zeros(Float32, n_V), ei = zeros(Float32, n_V), li=zeros(Float32,n_V), s_t=zeros(Float32, n_V), Li=zeros(Float32, n_V),
#            chg_power=zeros(Float32,n_V), set_physical = zeros(Float32,n_V), bus_stop_id=zeros(Int32,n_V), is_last_mile=zeros(Int32,n_V),
#            id_train=zeros(Int32,n_V), train_station=zeros(Int32,n_V))
#     for idx in 1:n_c
#         if is_last_mile[idx] == 0
#            V_DF.v_type[idx] = "bus stop"
#            V_DF.v_type[idx+n_c] = "transit station"
#         else
#             V_DF.v_type[idx] = "transit station"
#             V_DF.v_type[idx+n_c] = "bus stop"
#         end
#     end
#     length_PD = 2*n_c
#     V_DF.x[1:length_PD]= data_darp[2:end-1,2]
#     V_DF.y[1:length_PD]= data_darp[2:end-1,3]
#     V_DF.qi[1:length_PD]= data_darp[2:end-1,4]
#     V_DF.ei[1:length_PD]= data_darp[2:end-1,5]
#     V_DF.li[1:length_PD]= data_darp[2:end-1,6]
#     V_DF.s_t[1:length_PD]= data_darp[2:end-1,7]
#     V_DF.Li[1:length_PD]= data_darp[2:end-1,8]
#     V_DF.chg_power[1:length_PD] .= 0
#     V_DF.set_physical[1:length_PD] = data_darp[2:end-1,9]
#     V_DF.id_train[1:length_PD]= data_darp[2:end-1,10]
#     V_DF.bus_stop_id[1:length_PD] =data_darp[2:end-1, 11]
#     V_DF.is_last_mile[1:length_PD]= data_darp[2:end-1,12]
#     V_DF.train_station[1:length_PD] =data_darp[2:end-1, 13]
#     # set up charger nodes
#     set_id_chg = collect(length_PD+1:length_PD+n_charger)
#     V_DF.x[set_id_chg]= coord_chg[:,1] ; V_DF.y[set_id_chg]= coord_chg[:,2]
#     V_DF.qi[set_id_chg] .=0;V_DF.ei[set_id_chg] .=0
#     V_DF.li[set_id_chg].= 0 ;V_DF.chg_power[set_id_chg] = parameter_energy.pw_charger 
#     V_DF.s_t[set_id_chg] .=0;V_DF.Li[set_id_chg] .=0
#     V_DF.set_physical[set_id_chg] = v_chg_station_physical[set_physical_charger]
#     V_DF.id_train[set_id_chg] .= -1 
#     V_DF.v_type[set_id_chg] =  ["charger" for i in 1:n_charger]
#     V_DF.bus_stop_id[set_id_chg] .=-1;     V_DF.is_last_mile[set_id_chg] .=-1
#     V_DF.train_station[set_id_chg] .=-1
#     #depot
#     V_DF.v_type[end] = "depot"
#     V_DF.x[end]= data_darp[end,2]
#     V_DF.y[end]= data_darp[end,3]
#     V_DF.qi[end]= data_darp[end,4]
#     V_DF.ei[end]= data_darp[end,5]
#     V_DF.li[1end]= data_darp[end,6]
#     V_DF.s_t[end]= data_darp[end,7]
#     V_DF.Li[end]= data_darp[end,8]
#     V_DF.chg_power[end] = 0
#     V_DF.set_physical[end] = data_darp[end,9]
#     V_DF.id_train[end]= data_darp[end,10]
#     V_DF.bus_stop_id[end] =data_darp[end, 11]
#     V_DF.is_last_mile[1end]= data_darp[end,12]
#     V_DF.train_station[end] =data_darp[end, 13]

#     if flag_dev == 1
#         XLSX.writetable("V_DF.xlsx", overwrite=true, V_DF)
#         res_data_darp = DataFrame(node_id_DARP=collect(Int32, 1:total_node_darp),  x=data_darp[:,2] , y=data_darp[:,3], qi=data_darp[:,4], 
#             ei=data_darp[:,5], li = data_darp[:,6], s_t = data_darp[:,7], Li=data_darp[:,8], node_physical=data_darp[:,9], id_train=data_darp[:,10],
#             bus_stop_id= data_darp[:,11], is_last_mile = data_darp[:,12], train_station = data_darp[:,13])
#         XLSX.writetable("res_data_darp.xlsx", overwrite=true, res_data_darp)
#         dist=DataFrame(dist , :auto)
#         XLSX.writetable("res_dist.xlsx",overwrite=true, dist)
#         dist_orig=DataFrame(dist_orig , :auto)
#         XLSX.writetable("res_dist_orig.xlsx", overwrite=true, dist_orig)
#         res_dist_all=DataFrame(dist_all , :auto)
#         XLSX.writetable("res_dist_all.xlsx", overwrite=true, res_dist_all)
#     end

#     return   darp, instance, V_DF, timeTable 
# end

# function get_lyrs_compatible(timeTable, dist_all, data_darp,set_l_r_new, key_lyr, is_last_mile, vec_used_lyr,tw_width,v_k )
    
#     n_c_lyr = length()
#     lyrs_compatible= falses(n_c_lyr,n_c_lyr) 
#     for i in 1:n_c_lyr
#         lyrs_compatible[i,i]= true
#     end
#     for i in 1:n_c_lyr-1 # lower lyer
#         t_arr_i= timeTable.arrival_time[vec_used_lyr[i]]
#         ts_i   = timeTable.id_train_station[vec_used_lyr[i]] 
#         for j in i+1:n_c_lyr # upper lyer
#             t_arr_j= timeTable.arrival_time[vec_used_lyr[j]]
#             ts_j   = timeTable.id_train_station[vec_used_lyr[j]] 
#             _temp = coord_transit[ts_i,:]-coord_transit[ts_j,:]
#             _dist_tt= sqrt(sum(_temp.^2)) / v_k
#             @show(_dist_tt)
#             # @show(vec_used_lyr[i],vec_used_lyr[j], timeTable[vec_used_lyr[i],3],timeTable[vec_used_lyr[j],2],_dist_tt)
#             # ei+tij =<lj && li+tij >= ej (the criteria if two first-mile layers are compatible)
#             if is_last_mile[i]==0 && is_last_mile[j]==0
#                 ei,li = t_arr_i - tw_width, t_arr_i
#                 ej,lj=  t_arr_j - tw_width, t_arr_j
#                 if  ei + _dist_tt <= lj && li + _dist_tt >= ej
#                     lyrs_compatible[i,j]=true; lyrs_compatible[j,i]=true
#                 end
#             end
            
#         end
#     end
#     return lyrs_compatible
    
# end