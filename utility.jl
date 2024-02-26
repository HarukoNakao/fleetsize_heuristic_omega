########################################################################################
# compute the dist matrix of physical bus stops, transit station, chargers and the depot (dist_all)  
# example: V_l = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 37, 38, 39] => physical bus nodes (1:10), physical transit nodes (11, 12), chargers (37, 38), and the depot (39)
########################################################################################

function comp_dist(coord_V::Matrix{Float32},v_k)
 
    n = size(coord_V)[1]
    dist=zeros(Float32, n, n) 
    for i in 1:n-1
        for j = i+1:n
             _temp = coord_V[i,:]-coord_V[j,:]
             dist[i,j] = sqrt(sum(_temp.^2)) / v_k
             dist[j,i] = dist[i,j] 
        end
    end 
    return   dist 
end

# efficient check charging events conflicts of two routes
function  chg_conflict_two_routes(solution::Solution, ri, rj, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval)
 
    flag_chg =  solution.vec_chg_event[ri, 1] * solution.vec_chg_event[rj, 1]
    chg_conflict = false
    if flag_chg > 0
        vec_chg_event = solution.vec_chg_event
        occ_state_chg_ri .*= 0; occ_state_chg_rj .*=  0
        for i in 1:floor(Int32, vec_chg_event[ri,2])
            idx_chg, t1, t2 = vec_chg_event[ri, 2+(i-1)*4+2], vec_chg_event[ri, 2+(i-1)*4+3],vec_chg_event[ri, 2+(i-1)*4+4]
            idx_chg = floor(Int32,idx_chg)
            h_t1, h_t2 = floor(Int32, t1/discrte_t_interval)+1, floor(Int32, t2/discrte_t_interval)+1
            occ_state_chg_ri[idx_chg, h_t1: h_t2] .+=  1
        end     
        for i in 1:floor(Int32, vec_chg_event[rj,2])
            idx_chg, t1, t2 = vec_chg_event[rj, 2+(i-1)*4+2], vec_chg_event[rj, 2+(i-1)*4+3],vec_chg_event[rj, 2+(i-1)*4+4]
            idx_chg = floor(Int32,idx_chg)
            h_t1, h_t2 = floor(Int32, t1/discrte_t_interval)+1, floor(Int32, t2/discrte_t_interval)+1
            occ_state_chg_rj[idx_chg, h_t1: h_t2] .+=  1
        end
        if sum(occ_state_chg_ri .* occ_state_chg_rj) > 0 
             chg_conflict = true 
        end
        return chg_conflict       
    else
         return false # no conflict
    end
end


# efficient check charging events conflicts of a solution
function check_chg_occ_constr(solution::Solution, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval)
    
    n_route = solution.n_route
    total_add_chg_time = sum(solution.vec_chg_event[1:n_route, 1]) 
    if total_add_chg_time > 0
        idxs = findall(x->x>0, solution.vec_chg_event[1:n_route, 2]) # list of routes have charging events
        n_route_to_check = length(idxs)
        if n_route_to_check > 1
            for u in 1:n_route_to_check-1
                r = idxs[u] 
                occ_state_chg_ri .*= 0
                for i in 1:floor(Int32, solution.vec_chg_event[r,2])
                    vec_tmp = solution.vec_chg_event[r, 2+(i-1)*4+2 : 2+(i-1)*4+4] #idx_chg, t1, t2
                    idx_chg = floor(Int32,vec_tmp[1])
                    h_t1, h_t2 = floor(Int32, vec_tmp[2]/discrte_t_interval)+1, floor(Int32, vec_tmp[3]/discrte_t_interval)+1
                    occ_state_chg_ri[idx_chg, h_t1: h_t2] .+=  1
                end                
                for v in u+1:n_route_to_check
                    r = idxs[v]
                    occ_state_chg_rj .*= 0
                    for i in 1:floor(Int32, solution.vec_chg_event[r,2]) 
                        vec_tmp = solution.vec_chg_event[r, 2+(i-1)*4+2 : 2+(i-1)*4+4] #idx_chg, t1, t2
                        idx_chg = floor(Int32,vec_tmp[1])
                        h_t1, h_t2 = floor(Int32, vec_tmp[2]/discrte_t_interval)+1, floor(Int32, vec_tmp[3]/discrte_t_interval)+1
                        occ_state_chg_rj[idx_chg, h_t1: h_t2] .+=  1
                    end
                    sum(occ_state_chg_ri .* occ_state_chg_rj)>0 && (return false)
                end
            end
            return true
        else    
            return true
        end                 
    else
        return true
    end
end
 

# for verification of solution
function check_chg_occ_constr(solution::Solution, max_n_fast_chg, instance::Instance)
    
    T_max = instance.T_max
    discrte_t_interval = 1/6 # 10 second
    n_t_max = floor(Int32, (T_max+20)/discrte_t_interval) # 840 discrete charging occupancy state over T_max+20 minutes
    occ_state_chg = zeros(Int32, max_n_fast_chg, n_t_max)
    count_tmp=0
    vec_chg_event = solution.vec_chg_event
    n_route = solution.n_route
    total_add_chg_time = sum(solution.vec_chg_event[1:n_route, 1])
    if total_add_chg_time > 0 
        for r in 1:solution.n_route
            for i in 1:floor(Int32, vec_chg_event[r,2])
                count_tmp += 1
                idx_chg, t1, t2 = vec_chg_event[r, 2+(i-1)*4+2], vec_chg_event[r, 2+(i-1)*4+3],vec_chg_event[r, 2+(i-1)*4+4]
                idx_chg = floor(Int32,idx_chg)
                h_t1, h_t2 = floor(Int32, t1/discrte_t_interval)+1, floor(Int32, t2/discrte_t_interval)+1
                occ_state_chg[idx_chg, h_t1: h_t2] .+=  1
            end
        end
        if count_tmp < 2
            return true
        else
            if maximum(occ_state_chg)>1 
                error("charging operation conflicts between routes")
                return false
            else
                return true
            end 
        end
    end
    return true
end

function get_cost_chg_infra(fast_chg, global_parameter, n_chg_station)

    total_chg_infra_cost =0
    for s = 1:n_chg_station
        if fast_chg.num_chg_installed_by_site[s] > 0
            total_chg_infra_cost += (global_parameter.parameter_chg_infra.vec_cs_site_cost[s] 
                                     + fast_chg.num_chg_installed_by_site[s]* fast_chg.cost_fast_chg[s])
        end 
    end
    return total_chg_infra_cost
end

function add_charger(fast_chg, chg_station_id)
    
    fast_chg.n_fast_chg_installed += 1 
    fast_chg.num_chg_installed_by_site[chg_station_id] += 1
    pos =  fast_chg.num_chg_installed_by_site[chg_station_id]
    idx_chg = fast_chg.list_chgs_idx_by_site[chg_station_id, pos]  
    push!(fast_chg.set_id_chg_installed, idx_chg)
    
end

function reset_cs(fast_chg, n_chg_station)

    fast_chg.n_fast_chg_installed = 0
    for s in 1:n_chg_station
        fast_chg.num_chg_installed_by_site[s] = 0
    end  
    empty!(fast_chg.set_id_chg_installed)
     
end



################################################################
# DARP solver using hybrid simulated annealing
################################################################
function solver(path_result, global_parameter, fast_chg, darp, instance, demand_data, nodes_info, timeTable, output_solver)
    
    parameter_energy = instance.parameter_energy
    Q, T_max, K_MAX= instance.Q_max, instance.T_max, instance.K_MAX
    ei, li, Li, qi = darp.ei, darp.li, darp.Li, darp.qi 
    n_cus, start_depot, end_depot, dist  = darp.n_cus, darp.start_depot, darp.end_depot,  darp.dist 
    n_run, penalty  =  global_parameter.n_run, global_parameter.penalty # penalty of one unserved customer (i.e. 40)
 
    n_init_sol = global_parameter.n_init_sol
    N_NODE = 2*n_cus+2 
    idx_sorted = zeros(Int32, K_MAX)
    # average distance between any two locations in the network under study 
    tmp = dist[findall(x->x>0,  dist)] 
    avg_dist_between_nodes= sum(tmp)/length(tmp) # need to use x->x>0
   
    ei_orig=copy(ei);li_orig=copy(li)
    
    e_r   = copy(ei)
    cap_r = zeros(Int32,N_NODE)
    l_r   = copy(li)
    maxcap_r = zeros(Int32,N_NODE)
    dt_r  = zeros(Float32,N_NODE) # travel time from the start_depot to node i on route r
    dt_r_back = zeros(Float32,N_NODE) # travel time from the end_depot back to node i
    
    # sort veh type by putting EV on the top of fleet. So the ev will always be used first    
    succ = ones(Int32, 2*n_cus+2)*-1
    pre =  ones(Int32, 2*n_cus+2)*-1
    RI =   zeros(Int32, K_MAX, 5) # a matrix contains route info: start node, end node, num of nodes, veh capacity (each route is a route of a veh), veh_id
   
    if !global_parameter.is_gasoline_fleet 
        idx_sorted = sortperm(parameter_energy.is_electric, rev=true)        
        RI[:,5] = idx_sorted # put the electric vehices on the top so they can be used first. the route exchange will be conducted when there are cost savings        
        RI[:,4] = Q[RI[:,5]] # vector of capacity of vehs
    else
        idx_sorted = sortperm(parameter_energy.is_electric )        
        RI[:,5] = idx_sorted # put the gasoline vehices on the top
        RI[:,4] = Q[RI[:,5]] # vector of capacity of vehs
    end
    unserved_users = Int32[]
   
    # set up charging event container
    max_n_chg_event_per_route = parameter_energy.max_num_chg; length_info_chg=4 # max num of chg operations allowed for each route, user-defined parameters
    # @show(max_n_chg_event_per_route)
    max_n_fast_chg = fast_chg.max_n_fast_chg
    vec_chg_event = zeros(Float32, K_MAX, 2+max_n_chg_event_per_route*length_info_chg) # additional time (related to the same route without charging), num of chg event, precedent node of the oute, fast(rapid) charger index in the list of all fast charger, t_start_chg, t_end_chg
    discrte_t_interval = 1/6 # 10 second
    n_t_max = floor(Int32, (T_max+20)/discrte_t_interval) # 840 discrete charging occupancy state over T_max+20 minutes
    occ_state_chg_ri = zeros(Int32, max_n_fast_chg, n_t_max)
    occ_state_chg_rj = copy(occ_state_chg_ri) 
    succ[start_depot]=end_depot
    pre[end_depot]=start_depot 
    solu_init = Solution(succ, pre, RI, unserved_users, 0, 0, 0, 0, 0, vec_chg_event, 0) 
    soluT     = deepcopy(solu_init) # temporary solution
    solu_best = deepcopy(solu_init) # best solution for one run of SA 
    best_solu = deepcopy(solu_init) # best solution for multiple runs of SA 
    obj_values=zeros(Float32, n_run); co2_velues=zeros(Float32, n_run)
    obj_best_sol = Inf;  charging_time_best_sol=0
    sum_obj = 0; sum_co2=0; co2_best_sol=0
    # idx_run=1
    for idx_run = 1:n_run 
        # @show(idx_run)
        ei=copy(ei_orig); li=copy(li_orig)
        solu    = deepcopy(solu_init)
        # gen init sol
         e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back = gen_init_solution(solu, solu_init, global_parameter, fast_chg, darp, instance,
                             n_init_sol, Q, N_NODE, start_depot, end_depot, ei, li, K_MAX, n_cus, qi, Li, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval) 
        # verify_solution(solu, instance, darp, fast_chg)
        solu_best= SA(solu, soluT, solu_init, demand_data,  instance, fast_chg, darp, global_parameter, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, avg_dist_between_nodes, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval)
        obj_val = solu_best.total_cost
       
        solu_best.total_chg_time = get_sol_chg_time(solu_best, fast_chg)
        #check co2 target
        # if !global_parameter.is_gasoline_fleet
        #     !check_co2_emission_target(solu_best, global_parameter) && error("check_co2_emission_target failed !")
        # end
        sum_obj += obj_val
        sum_co2 +=  solu_best.co2_emission 
        if  obj_val < obj_best_sol
            best_solu  = solu_best; obj_best_sol = obj_val
            charging_time_best_sol = solu_best.total_chg_time 
            co2_best_sol           =  solu_best.co2_emission
        end  
        obj_values[idx_run] = obj_val 
        co2_velues[idx_run] = solu_best.co2_emission 
    end
    avg_obj = sum_obj/n_run
    avg_co2 = sum_co2/n_run
    if n_run>1
        std_obj, std_co2 = std(obj_values), std(co2_velues)
        cv_obj, cv_co2 = std_obj/avg_obj, std_co2/avg_co2
    else
        std_obj=std_co2=cv_obj=cv_co2=0
    end
    output_solver.obj_values = obj_values
    output_solver.n_run, output_solver.avg_obj, output_solver.avg_co2 = n_run, avg_obj, avg_co2;
    output_solver.std_obj, output_solver.std_co2, output_solver.cv_obj= std_obj, std_co2, cv_obj;
    output_solver.cv_co2, output_solver.obj_best_sol                  = cv_co2, obj_best_sol;
    output_solver.co2_best_sol, output_solver.charging_time_best_sol  = co2_best_sol, charging_time_best_sol ;
    return best_solu 
end

function generate_V_DF(data_darp, darp, nodes_info, instance, demand_data, n_c, n_V, fast_chg) 
   
    is_last_mile = demand_data.first_last_mile  
    dist, dist_orig, dist_all = darp.dist, darp.dist_orig, darp.dist_all
    parameter_energy = instance.parameter_energy
    max_n_fast_chg, v_chg_physical = fast_chg.max_n_fast_chg, fast_chg.v_physical_chg
    coord_chg = nodes_info.coord_chg_station[fast_chg.idxs_chg_station,:]
    total_node_darp = 2*n_c +2

    V_DF = DataFrame(id=collect(Int32, 1:n_V), v_type=Vector{String}(undef, n_V), x=zeros(Float32,n_V), y=zeros(Float32,n_V),
           qi = zeros(Float32, n_V), ei = zeros(Float32, n_V), li=zeros(Float32,n_V), s_t=zeros(Float32, n_V), Li=zeros(Float32, n_V),
           chg_power=zeros(Float32,n_V), set_physical = zeros(Float32,n_V), bus_stop_id=zeros(Int32,n_V), is_last_mile=zeros(Int32,n_V),
           id_train=zeros(Int32,n_V), train_station=zeros(Int32,n_V))
    for idx in 1:n_c
        if is_last_mile[idx] == 0
           V_DF.v_type[idx] = "bus stop"
           V_DF.v_type[idx+n_c] = "transit station"
        else
            V_DF.v_type[idx] = "transit station"
            V_DF.v_type[idx+n_c] = "bus stop"
        end
    end
    length_PD = 2*n_c
    V_DF.x[1:length_PD]= data_darp[2:end-1,2]
    V_DF.y[1:length_PD]= data_darp[2:end-1,3]
    V_DF.qi[1:length_PD]= data_darp[2:end-1,4]
    V_DF.ei[1:length_PD]= data_darp[2:end-1,5]
    V_DF.li[1:length_PD]= data_darp[2:end-1,6]
    V_DF.s_t[1:length_PD]= data_darp[2:end-1,7]
    V_DF.Li[1:length_PD]= data_darp[2:end-1,8]
    V_DF.chg_power[1:length_PD] .= 0
    V_DF.set_physical[1:length_PD] = data_darp[2:end-1,9]
    V_DF.id_train[1:length_PD]= data_darp[2:end-1,10]
    V_DF.bus_stop_id[1:length_PD] =data_darp[2:end-1, 11]
    V_DF.is_last_mile[1:length_PD]= data_darp[2:end-1,12]
    V_DF.train_station[1:length_PD] =data_darp[2:end-1, 13]
    # set up charger nodes
    set_id_chg = collect(length_PD+1:length_PD+max_n_fast_chg)
    V_DF.x[set_id_chg]= coord_chg[:,1] ; V_DF.y[set_id_chg]= coord_chg[:,2]
    V_DF.qi[set_id_chg] .=0;V_DF.ei[set_id_chg] .=0
    V_DF.li[set_id_chg].= 0 ;V_DF.chg_power[set_id_chg] = parameter_energy.pw_charger 
    V_DF.s_t[set_id_chg] .=0;V_DF.Li[set_id_chg] .=0
    V_DF.set_physical[set_id_chg] .= v_chg_physical 
    V_DF.id_train[set_id_chg] .= -1 
    V_DF.v_type[set_id_chg] =  ["charger" for i in 1:max_n_fast_chg]
    V_DF.bus_stop_id[set_id_chg] .=-1;     V_DF.is_last_mile[set_id_chg] .=-1
    V_DF.train_station[set_id_chg] .=-1
    #depot
    V_DF.v_type[end] = "depot"
    V_DF.x[end]= data_darp[end,2]
    V_DF.y[end]= data_darp[end,3]
    V_DF.qi[end]= data_darp[end,4]
    V_DF.ei[end]= data_darp[end,5]
    V_DF.li[1end]= data_darp[end,6]
    V_DF.s_t[end]= data_darp[end,7]
    V_DF.Li[end]= data_darp[end,8]
    V_DF.chg_power[end] = 0
    V_DF.set_physical[end] = data_darp[end,9]
    V_DF.id_train[end]= data_darp[end,10]
    V_DF.bus_stop_id[end] =data_darp[end, 11]
    V_DF.is_last_mile[1end]= data_darp[end,12]
    V_DF.train_station[end] =data_darp[end, 13]

    path_result     = pwd()*"\\result\\"
    XLSX.writetable(path_result * "V_DF.xlsx", overwrite=true, V_DF)
    res_data_darp = DataFrame(node_id_DARP=collect(Int32, 1:total_node_darp),  x=data_darp[:,2] , y=data_darp[:,3], qi=data_darp[:,4], 
        ei=data_darp[:,5], li = data_darp[:,6], s_t = data_darp[:,7], Li=data_darp[:,8], node_physical=data_darp[:,9], id_train=data_darp[:,10],
        bus_stop_id= data_darp[:,11], is_last_mile = data_darp[:,12], train_station = data_darp[:,13])
    XLSX.writetable(path_result * "res_data_darp.xlsx", overwrite=true, res_data_darp)
    dist=DataFrame(dist , :auto)
    XLSX.writetable(path_result * "res_dist.xlsx",overwrite=true, dist)
    dist_orig=DataFrame(dist_orig , :auto)
    XLSX.writetable(path_result * "res_dist_orig.xlsx", overwrite=true, dist_orig)
    res_dist_all=DataFrame(dist_all , :auto)
    XLSX.writetable(path_result * "res_dist_all.xlsx", overwrite=true, res_dist_all)
 
end

function generate_darp_instance(demand_data, instance, global_parameter, nodes_info, timeTable, fast_chg, darp)
    
    max_n_charger, v_chg_physical = fast_chg.max_n_fast_chg, fast_chg.v_physical_chg
    is_last_mile, no_person, id_train = demand_data.first_last_mile, demand_data.no_person, demand_data.id_train
    id_bus_stop, id_train_station     = demand_data.id_bus_stop, demand_data.id_train_station
    v_ts_physical, v_depot_physical = nodes_info.v_ts_physical, nodes_info.v_depot_physical
    coord_all, dist_all = nodes_info.coord_all, nodes_info.dist_all
    detour_factor, tw_width, bigM = global_parameter.detour_factor, global_parameter.tw_width, global_parameter.bigM
    flag_init = global_parameter.flag_initialization
    n_c, t_s, T_max = instance.n_c, instance.t_s, instance.T_max
    
    T_start = 0
    total_node_darp = 2 * n_c + 2 # total node of a DARP instance when you have _N requests (the depot nodes: 1, 2n+2, pickup: 2:n+1, dropoff: n+2: 2n+1)
    n_col_darp = 13 # num of columns in data_darp 
    data_darp=zeros(Float32, total_node_darp, n_col_darp)
    set_physical = zeros(Int32,total_node_darp+max_n_charger)
    # see Cordeau (2003) A Branch-and-Cut Algorithm for the Dial-a-Ride Problem
    # Time window tightening
    for idx in 1:n_c # idx is idx in DARP instance, pickup: 2:n+1, drop-off: n+2:2n+1, depot: 1,2n+2
        qi= no_person[idx]; bus_stop_id = id_bus_stop[idx]
        if is_last_mile[idx] == 0 #first mile
            p_physical=  id_bus_stop[idx]
            d_physical=  v_ts_physical[id_train_station[idx]]
            Li_tmp = dist_all[p_physical, d_physical] * detour_factor
            # Li_tmp =10
            li_d= timeTable.arrival_time[id_train[idx]]
            ei_d=li_d-tw_width
            ei_p = max(0, ei_d - Li_tmp -t_s) 
            li_p = min(li_d - dist_all[p_physical, d_physical] - t_s, T_max)  
            bus_stop_p, bus_stop_d = bus_stop_id, -1
        else
            p_physical=  v_ts_physical[id_train_station[idx]]
            d_physical=  id_bus_stop[idx]
            Li_tmp = dist_all[p_physical, d_physical] * detour_factor
            # Li_tmp =10
            ei_p =  timeTable.arrival_time[id_train[idx]]
            li_p =  ei_p + tw_width
            ei_d =  max(0, ei_p + t_s + dist_all[p_physical, d_physical])
            li_d =  min(li_p + t_s + Li_tmp, T_max)   
            bus_stop_p, bus_stop_d = -1, bus_stop_id
        end
        train_id =  id_train[idx]
        train_station_id = timeTable.id_train_station[train_id]
        # the first row is the depot
        # node_id_darp, x, y, demand,ei,li,service time, max_ridetime, physical_bus_node,     train_id, bus_node_id
        data_darp[idx+1,:]=   [idx+1,      coord_all[p_physical,1],coord_all[p_physical,2], qi,  ei_p, li_p, t_s, Li_tmp, p_physical, train_id, bus_stop_p, is_last_mile[idx],train_station_id]
        # create the corresponding drop-off node
        data_darp[idx+1+n_c,:]=[n_c+idx+1, coord_all[d_physical,1],coord_all[d_physical,2], -qi, ei_d, li_d, 0, Li_tmp, d_physical, train_id, bus_stop_d, is_last_mile[idx],train_station_id]
    end
    data_darp[1,:]   = [1,               coord_all[v_depot_physical,1], coord_all[v_depot_physical,2], 0, T_start, T_max, 0,0, v_depot_physical, -1,-1,-1,-1] 
    data_darp[end,:] = [total_node_darp, coord_all[v_depot_physical,1], coord_all[v_depot_physical,2], 0, T_start, T_max, 0,0, v_depot_physical, -1,-1,-1,-1] 
    
    n_cus, TH =  n_c, T_max # n_cus here is the total_node_darp
    start_depot, end_depot = 1, total_node_darp
    qi, ei, li, s_t, Li, set_physical[1:total_node_darp], is_last_mile_darp= floor.(Int32, data_darp[:,4]),  data_darp[:,5], data_darp[:,6], data_darp[:,7], data_darp[:,8], floor.(Int32,data_darp[:,9]), floor.(Int32,data_darp[:,12])
    bus_stop_id, train_id, train_station_id = floor.(Int32, data_darp[:,11]), floor.(Int32, data_darp[:,10]),floor.(Int32, data_darp[:,end])
    set_physical[total_node_darp+1:total_node_darp+max_n_charger] = copy(v_chg_physical)
    # preprocess(n_cus, li, dist_orig, s_t) 
    dist, dist_orig = set_dist(ei, li, s_t, Li, set_physical, dist_all, total_node_darp, bigM)
    preprocess(n_cus, li, dist_orig, s_t) 
    if flag_init
       darp=  DARP(n_cus, TH, start_depot, end_depot, ei, li, s_t, Li, qi, set_physical, is_last_mile_darp, bus_stop_id, train_id, train_station_id, dist, dist_orig, dist_all)
    else
        darp.ei, darp.li, darp.s_t, darp.Li, darp.qi, darp.set_physical, darp.is_last_mile = ei, li, s_t, Li, qi, set_physical, is_last_mile_darp
        darp.bus_stop_id, darp.train_id, darp.train_station_id, darp.dist, darp.dist_orig = bus_stop_id, train_id, train_station_id, dist, dist_orig
    end 
    return darp, data_darp
end


function worst_user(solution::Solution, qi, route_i, start_depot, end_depot, n_cus, dist)
    
    user_selected=0; dist_saving_max =0; dist_saving=0 
    route = get_route_no_depot(solution, route_i)
    if length(route) == 2  
        user_selected = route[1]
    else
        users = route[findall(x->x<n_cus+2, route)]
        for v in users
            v_ni = v+n_cus
            pre_v = solution.pre[v]; suc_v = solution.succ[v]
            pre_v_ni =solution.pre[v_ni]; suc_v_ni = solution.succ[v_ni]
            dist_saving = dist[pre_v,v]+dist[v,suc_v]+dist[pre_v_ni,v_ni]+dist[v_ni,suc_v_ni]
                        - (dist[pre_v,suc_v] + dist[pre_v_ni,suc_v_ni] )
            if dist_saving >= dist_saving_max
                dist_saving_max = dist_saving
                user_selected  = v
            end
        end
    end
     
    return user_selected 
end

function show_result(avg_obj, obj_values_sol, best_obj_value, λ, solution::Solution, total_charging_time, start_depot, end_depot, cus_bus_assigned, fast_chg)
    
    # solution= best_solu
    n_route = solution.n_route 
    total_tt_bus =0; total_tt_walk =0    
    # tmp =solution.unserved_users .- 1 # index needs to reduce one 
    unserved_cus =solution.unserved_users
    # tmp = collect(solution.unserved_users)
    # unserved_bus_stops = nodes[tmp]
    used_veh= solution.RI[1:n_route, 5]
    # @show(cus_bus_assigned, used_veh)
    # unserved_node_DARP = collect(solution.unserved_users)
    # for r in 1:n_route
    #     # total_tt_walk  += solution.TT_walk_routes[r] 
    #     total_tt_bus += length_route(solution, r, start_depot, darp)
    #     route = get_route(solution, r, start_depot, end_depot)
    #     route_V = [nodes[node] for node in route]
    #     if solution.vec_chg_event[r,2]>0
    #         vec_chg_event = floor.(Int32, solution.vec_chg_event[r,:])
    #         for j in 1:vec_chg_event[2]
    #              idx1, idx2 = 2 + (j-1)*4 +1, 2 + (j-1)*4 +2
    #              idx_3 = findfirst(x->x==vec_chg_event[idx1], route)
    #              insert!(route,   idx_3+1, fast_chg.v_chg[vec_chg_event[idx2]] )
    #              insert!(route_V, idx_3+1, fast_chg.v_chg[vec_chg_event[idx2]] )
    #         end
    #     end 
    #     @show(r, route, route_V)
    # end 
    # @show(n_run, best_obj_value, avg_obj, std(obj_values_sol)/avg_obj*100, obj_values_sol, total_tt_bus, total_charging_time, total_tt_walk, solution.penalty_unserved, unserved_bus_stops, unserved_cus,
    #       solution.n_route)
          @show(n_run, best_obj_value, avg_obj, std(obj_values_sol)/avg_obj*100, obj_values_sol, total_tt_bus, total_charging_time, total_tt_walk, solution.penalty_unserved, unserved_bus_stops, unserved_cus,
          solution.n_route)
end


# setup the dist matrix for the DARP instance, containing only the two dummy depots and pickup and drop-off nodes 
# see Cordeau (2003) A Branch-and-Cut Algorithm for the Dial-a-Ride Problem for preprocessing 
function set_dist(ei, li, s_t, Li, set_physical, dist_all, total_node_darp, bigM)
    
    n_c= floor(Int32, (total_node_darp-2)/2)
 
    dist= zeros(Float32, total_node_darp, total_node_darp)
    
    for i in 2:total_node_darp-1
        n_i= set_physical[i]
        dist[1,i] = dist_all[n_i,end] 
    end
    dist[end,:] = dist[1,:]
    dist[:,1]   = dist[1,:]
    dist[:,end] = dist[:,1] 
    
    for i in 2:total_node_darp-2
        for j in i+1:total_node_darp-1
            n_i, n_j = set_physical[i], set_physical[j] 
            dist[i,j] = dist_all[n_i,n_j] 
            dist[j,i] = dist[i,j]
        end
    end
    
    dist_orig = copy(dist)
    # set large distance for infeasible arcsv, see Cordeau 2006, Operations Research 
    P=collect(2:n_c+1);D=collect(n_c+2 : 2*n_c+1)
    # Arc elimination
    PD=union(P,D)
    dist[end,:] .= bigM
    for d in D
        dist[1,d] = bigM
        dist[d, d-n_c] = bigM
    end
    for p in P
        dist[p,end] =  dist[p,1] = bigM
    end
    for i in PD 
        for j in PD
            if i!=j
                if ei[i]+s_t[i]+dist_orig[i,j] > li[j]
                    dist[i,j]= bigM
                end
            end
        end
    end
    for i in P 
        for j in PD
            if i!=j
                if dist_orig[i,j]+s_t[j]+dist_orig[j,i+n_c] > Li[i]
                    dist[i,j]=dist[j,n_c+i]=bigM
                end 
            end
        end
    end 
    return dist, dist_orig
end


#########################################
# efficient check whether there are double tours visiting the same transit dummmy node (we do not allow it)
# for example, if node 15 is a transit node (on a layered graph), a route: 1 5 15 17 8 15 is not allowed  
# this function is used to check a sequence for such situation
#########################################
# function no_double_tour(nodes_route, D′)
    
#     big_N = 3000 # sufficient large than n_V (total num. of nodes of the problem instance)
#     idxs = findall(x->x∈D′, nodes_route)
#     transit_nodes = nodes_route[idxs]
#     id =zeros(Int32, big_N)
#     for (idx, v) in enumerate(transit_nodes)
#         if id[v] == 0
#             id[v] = idxs[idx]
#         else
#             if idxs[idx] - id[v] > 1 
#                 return false
#             else
#                 id[v] = idxs[idx]
#             end
#         end
#     end
#     return true
# end

# function check_used_vehs(solution::Solution, K)

#     veh_used = zeros(Int32,K)
#     n_route= solution.n_route
#     vec_vehs = solution.RI[1:n_route,5]
#     for i in 1:n_route
#         veh_used[vec_vehs[i]] += 1
#     end
#     if !isempty(filter(x->x∉collect(1:K), vec_vehs)) || maximum(veh_used)> 1
#         @show(vec_vehs, veh_used, maximum(veh_used), filter(x->x∉collect(1:K), vec_vehs))
#         return false
#     end
#     return true
# end



# function compute_Li_walk_tt(Li, x_opt_array, cus_bus_assign, nodes_bus, coord_V, avg_wlk_speed, v_depot, flag_dev, R)

#     idxs = findall(x->x>0, cus_bus_assign) # idx of bus stops with positive requests
#     count = 1 ; total_tt_walk =0; count_cus=0; cus_bus_assigned = Any[]; walk_tt_bus_stop = zeros(Float32, length(idxs))
#     walk_tt_check = zeros(Float32, n_c, 5)
#     for idx in idxs

#         set_cus = findall(x->x>0, x_opt_array[:,idx]) 
#         n_cus = size(set_cus, 1)
#         push!(cus_bus_assigned, set_cus)
#         sum_tt_walk = 0
#         for k in 1:n_cus
#             count_cus += 1
#             g = nodes_bus[idx] 
#             # pay attention the node id of customers
#             vec = coord_V[v_depot + set_cus[k],:] - coord_V[g,:]
#             # @show(v_depot + set_cus[k], g)
#             tt_walk = (sqrt(sum(vec.^2)) / avg_wlk_speed)
#             total_tt_walk += tt_walk
#             sum_tt_walk += tt_walk
#             walk_tt_check[count_cus,:]=[set_cus[k], R[set_cus[k]], count+1, g, sqrt(sum(vec.^2)) / avg_wlk_speed]
#         end
#         walk_tt_bus_stop[count] = sum_tt_walk
#         count += 1
#     end 
    
#     # check
#     if flag_dev == 1
#         walk_tt_check_csv=DataFrame(cus_id = walk_tt_check[:,1], node_id_r = walk_tt_check[:,2], node_id_DARP = walk_tt_check[:,3], bus_stop = walk_tt_check[:,4],
#                        tt_walk = walk_tt_check[:,5] )
#         CSV.write("res_walk_tt_check.csv", walk_tt_check_csv)
#     end
#     return total_tt_walk, cus_bus_assigned, walk_tt_bus_stop
# end

# function compute_walk_tt_test(nodes_bus_assign,  set_cus, coord_V, avg_wlk_speed )

#     n_active_bs =  length(nodes_bus_assign)
#     walk_tt_bus_stop = zeros(Float32, n_active_bs)
 
#     for (idx, g) in enumerate(nodes_bus_assign)
#         sum_tt_walk = 0
#         for r in 1:length(set_cus[idx])
#             vec = coord_V[set_cus[idx][r],:] - coord_V[g,:]
#             tt_walk = (sqrt(sum(vec.^2)) / avg_wlk_speed)
#             sum_tt_walk += tt_walk
#         end
#         walk_tt_bus_stop[idx] = sum_tt_walk
#     end 

#     return  walk_tt_bus_stop
# end



# compute the dist matrix for customers and bus stop nodes
# function comp_dist_cus_bus(coord_V, nodes_bus, R, lyr_cus, lyr_bus, lyrs_compatible, avg_wlk_speed, bigM)

#     coord_r =coord_V[R, :]
#     coord_bus=coord_V[nodes_bus,:]
#     n_1 = size(coord_r,1)
#     n_2 = size(coord_bus,1)
#     dist=zeros(Float32, n_1, n_2) 
#     for i in 1:n_1
#         l_i = lyr_cus[i]        
#         for j = 1:n_2
#             vec = coord_r[i,:]-coord_bus[j,:]
#             dist[i,j] = sqrt(sum(vec.^2)) / avg_wlk_speed
#             l_j = lyr_bus[j]
#             if  l_i != l_j # need to handle each layer seperately for the cus assignment as the TW at customers destinations are different from one layer to another in general
#                 dist[i,j] = bigM 
#             end
#         end
#     end      
#     return dist  
# end


# function show_result_test(avg_obj, obj_values_sol, best_obj_value, λ, solution::Solution, total_charging_time, start_depot, end_depot, cus_bus_assigned, fast_chg)
    
#     start_depot =1;end_depot=2*n_cus+2
#     n_route = solution.n_route 
#     total_tt_bus =0; total_tt_walk =0    
#     # tmp = collect(solution.unserved_users) .- 1 # index needs to reduce one 
#     # unserved_cus = collect(cus_bus_assigned[tmp])
#     # tmp = collect(solution.unserved_users)
#     # unserved_bus_stops = nodes[tmp]
#     used_veh= solution.RI[1:n_route, 5]
#     # @show(cus_bus_assigned, used_veh)
#     # unserved_node_DARP = collect(solution.unserved_users)
#     for r in 1:n_route
#         total_tt_walk  += solution.TT_walk_routes[r] 
#         total_tt_bus += length_route(solution, r, start_depot, darp)
#         route = get_route(solution, r, start_depot, end_depot)
#         route_V = [nodes[node] for node in route]
#         if solution.vec_chg_event[r,2]>0
#             vec_chg_event = floor.(Int32, solution.vec_chg_event[r,:])
#             for j in 1:vec_chg_event[2]
#                  idx1, idx2 = 2 + (j-1)*4 +1, 2 + (j-1)*4 +2
#                  idx_3 = findfirst(x->x==vec_chg_event[idx1], route)
#                  insert!(route,   idx_3+1, fast_chg.v_chg[vec_chg_event[idx2]] )
#                  insert!(route_V, idx_3+1, fast_chg.v_chg[vec_chg_event[idx2]] )
#             end
#         end 
#         @show(r, route, route_V)
#     end 
#     @show(n_run, best_obj_value, avg_obj, std(obj_values_sol)/avg_obj*100, obj_values_sol, total_tt_bus, total_charging_time, total_tt_walk, solution.penalty_unserved,
#           solution.n_route)
# end


# # verification of the result of the bus stop assignment of customers 
# function verif_assign(x_opt_array, V, n_c, set_l_r_new, nodes_bus, lyrs_compatible)

#     temp = [findfirst(x->x>0, x_opt_array[i,:]) for i in 1:n_c]

#     if sum(V[nodes_bus[temp], 5] - set_l_r_new) > 0
#         temp = V[nodes_bus[temp], 5] - set_l_r_new
#         idxs = findall(x->x>0, temp)
#         for i in 1:length(idxs)
#             lyr1, lyr2= V[nodes_bus[i], 5], set_l_r_new[i]
#             if lyrs_compatible[lyr1, lyr2] == 0
#                 return false
#             end
#         end
#     end     
#     return true    
# end
