
# multiple dispatch function
9 

function co2_route(solution::Solution, route_i, start_depot, darp, parameter_energy)
    
    co2 = 0; dist_orig = darp.dist_orig # dist_orig is measured in travel time (minute)
    veh_id = solution.RI[route_i,5]
    vertex = solution.RI[route_i,1]
    co2 += dist_orig[start_depot, vertex] * parameter_energy.co2_emission[veh_id]
    for i in 1: solution.RI[route_i,3] 
        co2 += dist_orig[vertex, solution.succ[vertex]] * parameter_energy.co2_emission[veh_id]
        vertex = solution.succ[vertex]
    end
    return co2
end 

function co2_route_gv(solution::Solution, route_i, start_depot, darp, parameter_energy)
    
    veh_id_gv   = parameter_energy.veh_ids_types[1][1]
    co2 = 0; dist_orig = darp.dist_orig # dist_orig is measured in travel time (minute)
    vertex = solution.RI[route_i,1]
    co2 += dist_orig[start_depot, vertex] * parameter_energy.co2_emission[veh_id_gv]
    for i in 1: solution.RI[route_i,3] 
        co2 += dist_orig[vertex, solution.succ[vertex]] * parameter_energy.co2_emission[veh_id_gv]
        vertex = solution.succ[vertex]
    end
    return co2
end

function change_veh_type(solution::Solution, global_parameter, darp, parameter_energy, fast_chg, instance)

    co2_threshold, big_penalty_co2 = global_parameter.co2_threshold, global_parameter.big_penalty_co2
    penalty, penalty_per_veh_used  = global_parameter.penalty, global_parameter.penalty_per_veh_used
    sum_cost = 0;sum_co2 = 0;total_cost_with_penalty = 0;  qi = darp.qi;
    n_route =  solution.n_route 
    matrix_cost_v_type, matrix_co2_v_type = compute_matrix_cost_v_type(solution, darp, parameter_energy)
    for r in 1:n_route   
        cost, co2= cost_route(solution, r, darp, parameter_energy, fast_chg, instance)
        sum_cost += cost 
        sum_co2  += co2
    end
    unserved_users = solution.unserved_users
    total_penalty = sum(qi[unserved_users]) * penalty + n_route * penalty_per_veh_used #  penalty for num of used veh is not used (0)
    total_penalty_co2 = max(0, sum_co2 - co2_threshold) * big_penalty_co2
    solution.co2_emission = sum_co2 
    total_cost= sum_cost; total_cost_with_penalty= sum_cost + total_penalty + total_penalty_co2
  
    return total_cost, total_cost_with_penalty # total_cost is the obj_val for the lwer level problem

end


function check_capacity()
    cap_r[1]=0;maxcap_r[end_depot]=0
    v_i = pre_vi
    while v_j > 0
        e_r[v_j] = max(ei[v_j], e_r[v_i]+s_t[v_i]+dist[v_i,v_j])
        cap_r[v_j] =  cap_r[v_i] + qi[v_j]
        # dt_r[v_j]  =  dt_r[v_i]  + dist_orig[v_i,v_j]
        v_i,v_j = v_j, solution.succ[v_j]
    end

end

# inclduing the penalty of co2 emission when exceeding the targeted value
function cost_saving_exchange(solution::Solution, route_i, darp, parameter_energy, fast_chg)

    tt_route, dist_orig = 0, darp.dist_orig
    start_depot, end_depot, qi = darp.start_depot, darp.end_depot, darp.qi
    cap_r= zeros(end_depot)
    veh_id = solution.RI[route_i,5]
    vertex = solution.RI[route_i,1]
    cap_r[vertex] = qi[vertex]
   
    cost_minute = parameter_energy.beta[veh_id]*parameter_energy.e_price[veh_id]
    veh_id_gv   = parameter_energy.veh_ids_types[1][1]
    cost_minute_gv = parameter_energy.beta[veh_id_gv]*parameter_energy.e_price[veh_id_gv] # first veh is gv
    tt_route += dist_orig[start_depot, vertex] 
    for i in 1: solution.RI[route_i,3] 
        tt_route  += dist_orig[vertex, solution.succ[vertex]] 
        vertex_pre = vertex
        vertex     = solution.succ[vertex]
        cap_r[vertex] = cap_r[vertex_pre] + qi[vertex]
    end
    add_tt_acc_chg= 0 # get veh routing times for charging operations
    if solution.vec_chg_event[route_i,1]>0
        _, add_tt_acc_chg = get_add_tt_acc_chg(solution, darp, route_i, fast_chg, parameter_energy) # only include access time, not charging time, see the paper
    end 
    cost_ev = cost_minute * (tt_route + add_tt_acc_chg) + parameter_energy.fleet_cost[veh_id]
    cost_gv = cost_minute_gv * tt_route + parameter_energy.fleet_cost[veh_id_gv]
    cap_feasible = true
    cap_gv =  parameter_energy.cap_passenger[veh_id_gv] 
    if cap_gv < parameter_energy.cap_passenger[veh_id] # check capacity feasiblity
        ! isnothing(findfirst(x->x>cap_gv, cap_r)) && (cap_feasible = false)
    end
    if cap_feasible 
         return max(0, cost_ev - cost_gv)
    else
        # @show("gaso cap infeaible !", cap_feasible)
        return 0
    end
end

# inclduing the penalty of co2 emission when exceeding the targeted value
function cost_route(solution::Solution, route_i, darp, parameter_energy, fast_chg, instance)

    sum_co2_route, tt_route, dist_orig = 0, 0, darp.dist_orig
    start_depot, end_depot, dist_all = darp.start_depot, darp.end_depot, darp.dist_all
    set_physical = darp.set_physical 
    T_max, K = instance.T_max, instance.K 
    n_charger_installed=0; success=false;success_chg=false; pos_insert_chg=0; forward_time_slack=0; vec_t_start=0; route=0

    veh_id = solution.RI[route_i,5]
    vertex = solution.RI[route_i,1]
   
    cost_minute = parameter_energy.beta[veh_id]*parameter_energy.e_price[veh_id]
    tt_route += dist_orig[start_depot, vertex] 
    for i in 1: solution.RI[route_i,3] 
        tt_route += dist_orig[vertex, solution.succ[vertex]] 
        vertex = solution.succ[vertex]
    end
    add_tt_acc_chg= 0 # get veh routing times for charging operations
    if solution.vec_chg_event[route_i,1]>0
        success, add_tt_acc_chg = get_add_tt_acc_chg(solution, darp, route_i, fast_chg, parameter_energy) # only include access time, not charging time, see the paper
    end
 
    if !parameter_energy.is_electric[veh_id] 
        sum_co2_route = tt_route * parameter_energy.co2_emission[veh_id]
    end 
    # return cost_minute * (tt_route + add_tt_chg) + parameter_energy.fleet_cost[veh_id], sum_co2_route
    return cost_minute * (tt_route + add_tt_acc_chg) + parameter_energy.fleet_cost[veh_id], sum_co2_route
end


# inclduing the penalty of co2 emission when exceeding the targeted value
function cost_route_check(solution::Solution, route_i, darp, parameter_energy, fast_chg, instance)

    sum_co2_route, tt_route, dist_orig = 0, 0, darp.dist_orig
    start_depot, end_depot, dist_all = darp.start_depot, darp.end_depot, darp.dist_all
    set_physical = darp.set_physical 
    T_max, K = instance.T_max, instance.K 
    n_charger_installed=0; success=false;success_chg=false; pos_insert_chg=0; forward_time_slack=0; vec_t_start=0; route=0

    veh_id = solution.RI[route_i,5]
    vertex = solution.RI[route_i,1]
   
    cost_minute = parameter_energy.beta[veh_id]*parameter_energy.e_price[veh_id]
    tt_route += dist_orig[start_depot, vertex] 
    route=zeros(Int,solution.RI[route_i,3]+2)
    route[1:2]=[start_depot,vertex];pos_route =3
    for i in 1: solution.RI[route_i,3] 
        tt_route += dist_orig[vertex, solution.succ[vertex]] 
        vertex = solution.succ[vertex]
        route[pos_route] = vertex
        pos_route +=1
    end
    add_tt_acc_chg= 0 # get veh routing times for charging operations
    if solution.vec_chg_event[route_i,1]>0
        success, add_tt_acc_chg = get_add_tt_acc_chg(solution, darp, route_i, fast_chg, parameter_energy) # only include access time, not charging time, see the paper
    end
 
    if !parameter_energy.is_electric[veh_id] 
        sum_co2_route = tt_route * parameter_energy.co2_emission[veh_id]
    end 
    @show(route_i, route, cost_minute, tt_route, add_tt_acc_chg, parameter_energy.fleet_cost[veh_id] )
    return cost_minute * (tt_route + add_tt_acc_chg) + parameter_energy.fleet_cost[veh_id], sum_co2_route
end
 

function length_route(route, darp)
    
    tt_route = 0 
    dist = darp.dist
    for idx in 1:length(route)-1        
        tt_route += dist[route[idx], route[idx+1]]
    end  
    return tt_route
end

function show_route_detail(solution::Solution, route_i, start_depot, darp, parameter_energy)
    
    end_depot= darp.end_depot
    vertex = solution.RI[route_i,1]
    cost_minute = parameter_energy.beta[veh_id]*parameter_energy.e_price[veh_id]
    tt_route = 0; vec_dist=Float32[]
    tt_route += dist_orig[start_depot, vertex] 
    push!(vec_dist,dist_orig[start_depot, vertex] )
    for i in 1: solution.RI[route_i,3] 
        tt_route += dist_orig[vertex, solution.succ[vertex]] 
        push!(vec_dist, dist_orig[vertex, solution.succ[vertex]])
        vertex = solution.succ[vertex]
    end
    vec_dist_km = vec_dist .* v_k 
    tt_route_km = tt_route * v_k
    @show(get_route(solution, 1, start_depot, end_depot))
    @show(vec_dist, vec_dist_km, sum(vec_dist), tt_route, tt_route_km)
end

# function using real travel time (non penalized) for energy consumption calculation
function length_route_energy_check(route, darp)
    
    tt_route = 0 
    dist_orig = darp.dist_orig
    for idx in 1:length(route)-1        
        tt_route += dist_orig[route[idx], route[idx+1]]
    end  
    return tt_route
end

function check_co2_emission_target(solution::Solution, global_parameter)

  
    if solution.co2_emission > global_parameter.co2_threshold
        @show(solution.co2_emission , global_parameter.co2_threshold)
        return false
    else
        return true
    end

end


function show_co2_detail(solution::Solution, darp, instance, parameter_energy)

    start_depot = darp.start_depot; sum_co2 = 0
    for r in 1:solution.n_route
        vehid = solution.RI[r,5]
        dist_route_km = length_route(solution, r, darp) * instance.v_k
        co2 = co2_route(solution, r, start_depot, darp, parameter_energy)
        sum_co2 += co2
        @show(r, dist_route_km, parameter_energy.is_electric[vehid], parameter_energy.co2_emission_km[vehid], co2)
    end   
    @show(sum_co2, solution.co2_emission )
end

function get_route_chg_time_energy(solution::Solution, r, fast_chg::Fast_chargers)
   
    total_route_chg_time =0; total_route_chg_energy = 0
    for j = 1:floor(Int32, solution.vec_chg_event[r, 2])
        idxs = collect(2 +(j-1)*4+2 : 2 + j*4)
        chg_id = floor(Int32, solution.vec_chg_event[r, idxs[1]])
        dur_chg= solution.vec_chg_event[r, idxs[3]] - solution.vec_chg_event[r, idxs[2]]
        total_route_chg_time += dur_chg
        total_route_chg_energy += (dur_chg * fast_chg.pw[chg_id])
    end
    return total_route_chg_time, total_route_chg_energy
end

function get_sol_chg_time(solution::Solution, fast_chg::Fast_chargers)
    
    total_chg_time= 0
    for r in 1:solution.n_route  
        dur_chg, _ = get_route_chg_time_energy(solution, r, fast_chg)   
        total_chg_time +=  dur_chg      
    end
    return total_chg_time
end


#function to get all acc time to and from chg
function get_add_tt_acc_chg(solution::Solution, darp, r, fast_chg, parameter_energy)
 
    set_physical, dist_all = darp.set_physical, darp.dist_all
    start_depot, end_depot = darp.start_depot, darp.end_depot
    add_tt = 0 ; max_num_chg = parameter_energy.max_num_chg
    info_chg = zeros(Float32, 2+max_num_chg*4) # store at most 4 charging operations
    route = get_route(solution, r, start_depot, end_depot)
    for i in 1:floor(Int32, solution.vec_chg_event[r,2])
        v1    = floor(Int32, solution.vec_chg_event[r, 2+(i-1)*4+1])
        idx = findfirst(x->x==v1, route)
        if isnothing(idx) # charging schedule walks off, need to reschedule it
           return false, 0
        end 
        v2= route[idx+1]
        v1, v2   = set_physical[v1], set_physical[v2]
        chg_id   = floor(Int32,solution.vec_chg_event[r, 2+(i-1)*4+2]) # chg id on the list of used fast chargers
        v_physi_chg       =  fast_chg.v_physical_chg[chg_id]
        t_access_1 = dist_all[v1, v_physi_chg] 
        t_access_2 = dist_all[v_physi_chg, v2]  
        dist_v1_v2 = dist_all[v1, v2]
        add_tt += (t_access_1 + t_access_2 - dist_v1_v2) #correct 
    end
    return true, add_tt
end

function cost_solution(solution::Solution, global_parameter, darp, parameter_energy, fast_chg, instance)
    
    co2_threshold, big_penalty_co2 = global_parameter.co2_threshold, global_parameter.big_penalty_co2
    penalty, penalty_per_veh_used  = global_parameter.penalty, global_parameter.penalty_per_veh_used
    sum_cost = 0;sum_co2 = 0;total_cost_with_penalty = 0;  qi = darp.qi;
    n_route =  solution.n_route 
    for r in 1:n_route   
        cost, co2= cost_route(solution, r, darp, parameter_energy, fast_chg, instance)
        sum_cost += cost 
        sum_co2  += co2
    end
    unserved_users = solution.unserved_users
    total_penalty = sum(qi[unserved_users]) * penalty + n_route * penalty_per_veh_used #  penalty for num of used veh is not used (0)
    total_penalty_co2 = max(0, sum_co2 - co2_threshold) * big_penalty_co2
    solution.co2_emission = sum_co2 
    total_cost= sum_cost; total_cost_with_penalty= sum_cost + total_penalty + total_penalty_co2
    solution.total_cost, solution.total_cost_with_penalty = total_cost, total_cost_with_penalty
    # return total_cost, total_cost_with_penalty # total_cost is the obj_val for the lwer level problem
  
end


function cost_solution_check(solution::Solution, global_parameter, darp, parameter_energy, fast_chg, instance)
    
    # solution =solu_best
    start_depot = darp.start_depot
    co2_threshold, big_penalty_co2 = global_parameter.co2_threshold, global_parameter.big_penalty_co2
    penalty, penalty_per_veh_used  = global_parameter.penalty, global_parameter.penalty_per_veh_used
    sum_cost = 0;sum_co2 = 0;total_cost_with_penalty = 0;  qi = darp.qi;
    n_route =  solution.n_route 
    for r in 1:n_route   
        cost, co2= cost_route_check(solution, r, darp, parameter_energy, fast_chg, instance)
        sum_cost += cost 
        sum_co2  += co2
    end
    unserved_users = solution.unserved_users
    total_penalty = sum(qi[unserved_users]) * penalty + n_route * penalty_per_veh_used #  penalty for num of used veh is not used (0)
    total_penalty_co2 = max(0, sum_co2 - co2_threshold) * big_penalty_co2
    solution.co2_emission = sum_co2 
    total_cost= sum_cost; total_cost_with_penalty= sum_cost + total_penalty + total_penalty_co2
    @show( sum(qi[unserved_users]))
    return total_cost, total_cost_with_penalty # total_cost is the obj_val for the lwer level problem
  
end

 

function get_route(solution::Solution, route_i, start_depot, end_depot)
    
    n_nodes = solution.RI[route_i,3] + 2
    route = zeros(Int32,n_nodes)
    route[1]= start_depot; route[end]= end_depot
    route[2] = solution.RI[route_i,1]
    for i in 2:n_nodes-2
        route[i+1] = solution.succ[route[i]]
    end
    return route 
    
end

function get_route_no_depot(solution::Solution, route_i)
    
    n_nodes = solution.RI[route_i,3]  
    route = zeros(Int32,n_nodes)
    route[1] = solution.RI[route_i,1]
    for i in 1:n_nodes-1
        route[i+1] = solution.succ[route[i]]
    end
    return route 
    
end

function compute_dt_r(route::Vector{Int32}, dt_r, dist_orig)
    
    dt_r[route[1]]=0    
    for i in 1:length(route)-1
        dt_r[route[i+1]] = dt_r[route[i]] + dist_orig[route[i], route[i+1]]
    end
end


# v_i,v_j, j is the inserted node, pos_insert is the position index after v_i
function update_auxiliary_var(solution::Solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
    
    n_cus, ei, li, qi = darp.n_cus, darp.ei, darp.li, darp.qi
    s_t, dist = darp.s_t, darp.dist

    pre_vi =best_pos[1];user=best_pos[2];suc_vi=best_pos[3]
    pre_v_ni = best_pos[4]; v_ni = best_pos[5]; suc_v_ni=best_pos[6]
    v_j= user # pickup node
    end_depot = 2*n_cus+2
    # forward updates
    cap_r[1]=0;maxcap_r[end_depot]=0
    v_i = pre_vi
    while v_j > 0
        e_r[v_j] = max(ei[v_j], e_r[v_i]+s_t[v_i]+dist[v_i,v_j])
        cap_r[v_j] =  cap_r[v_i] + qi[v_j]
        # dt_r[v_j]  =  dt_r[v_i]  + dist_orig[v_i,v_j]
        v_i,v_j = v_j, solution.succ[v_j]
    end
    # backward updates
    v_j = suc_v_ni 
    v_i = v_ni
    while v_i > 0 
        l_r[v_i] = min(li[v_i], l_r[v_j] - dist[v_j,v_i] - s_t[v_i])
        # l_r[v_i] = li[v_i]
        maxcap_r[v_i] = max(0, maxcap_r[v_j] + qi[v_i])
        # dt_r_back[v_i] = dt_r_back[v_j] + dist_orig[v_i,v_j]
        v_i,v_j = solution.pre[v_i], v_i
    end
    
end

function route_update(solution::Solution, route_i::Vector{Int32}, ri, n_cus)
    
    start_depot = 1
    end_depot = 2*n_cus+2
    route = copy(route_i[2:end-1]) 
    n_nodes = length(route)
    solution.RI[ri,1] = route[1] ; solution.RI[ri,2] = route[end] 
    solution.RI[ri,3] = n_nodes 
    solution.pre[route[1]] = start_depot
    solution.succ[route[1]] = route[2]
    for i in 2:n_nodes - 1
        solution.pre[route[i]] = route[i-1]
        solution.succ[route[i]] = route[i+1]
    end
    solution.pre[route[end]] = route[end-1]
    solution.succ[route[end]] = end_depot 

end

function route_auxiliary_update(solution::Solution, darp, route_i::Vector{Int32}, ri, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
    
    start_depot, end_depot = darp.start_depot, darp.end_depot
    route = copy(route_i[2:end-1])
    # @show(route)
    n_nodes = length(route)
    solution.RI[ri,1] = route[1] ; solution.RI[ri,2] = route[end] 
    solution.RI[ri,3] = n_nodes
    solution.pre[route[1]] = start_depot
    solution.succ[route[1]] = route[2]
    for i in 2:n_nodes - 1
        solution.pre[route[i]] = route[i-1]
        solution.succ[route[i]] = route[i+1]
    end
    solution.pre[route[end]] = route[end-1]
    solution.succ[route[end]] = end_depot
    reset_auxiliary_route(solution, darp, ri, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)

end


function route_auxiliary_update_no_depot(solution::Solution, route::Vector{Int32}, ri, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
    
    start_depot =1
    end_depot = 2 * n_cus + 2
    # @show(route)
    n_nodes = length(route)
    solution.RI[ri,1] = route[1] ; solution.RI[ri,2] = route[end] 
    solution.RI[ri,3] = n_nodes
    solution.pre[route[1]] = start_depot
    solution.succ[route[1]] = route[2]
    for i in 2:n_nodes - 1
        solution.pre[route[i]] = route[i-1]
        solution.succ[route[i]] = route[i+1]
    end
    solution.pre[route[end]] = route[end-1]
    solution.succ[route[end]] = end_depot
    reset_auxiliary_route(solution, darp, ri, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)

end

function reset_cap_r(solution::Solution, route_i, cap_r, darp)
 
    start_depot = 1
    cap_r[start_depot] = 0 
    v_i = start_depot
    v_j = solution.RI[route_i,1]
    while v_j > 0
        cap_r[v_j] =  cap_r[v_i] + darp.qi[v_j]
        v_i,v_j = v_j, solution.succ[v_j]
    end
end

function reset_auxiliary_route(solution::Solution, darp, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)

    start_depot, end_depot, ei, li, qi, s_t, dist = darp.start_depot, darp.end_depot, darp.ei, darp.li, darp.qi, darp.s_t, darp.dist 
     
    e_r[start_depot]= 0
    l_r[end_depot]  = li[end_depot]
    cap_r[start_depot] = 0 ; maxcap_r[end_depot] = 0
    # dt_r[start_depot] = 0; dt_r_back[end_depot] = 0
    v_i = start_depot
    v_j = solution.RI[route_i,1]
    while v_j > 0
        e_r[v_j] = max(ei[v_j], e_r[v_i]+s_t[v_i]+dist[v_i,v_j])
        cap_r[v_j] =  cap_r[v_i] + qi[v_j]
        # dt_r[v_j] = dt_r[v_i] + dist_orig[v_i,v_j]
        v_i,v_j = v_j, solution.succ[v_j]
    end
    # backward updates
    v_j = end_depot# start from the next node of the drop off node
    v_i = solution.RI[route_i,2]
    while v_i > 0 
        l_r[v_i] = min(li[v_i], l_r[v_j] - dist[v_j,v_i]-s_t[v_i])
        maxcap_r[v_i] = max(0, maxcap_r[v_j] + qi[v_i])
        # dt_r_back[v_i] = dt_r_back[v_j] + dist_orig[v_i,v_j]
        v_i,v_j = solution.pre[v_i], v_i
    end

end


function insert_consecutive(solution::Solution, route_i,  start_depot, user, end_depot, n_cus, qi) # user is pickup node is from 2,...,n+1
    
    solution.RI[route_i,1] = user
    solution.RI[route_i,2] = user+n_cus
    solution.RI[route_i,3] += 2 # pickup and dropoff nodes
    solution.succ[start_depot]=user
    solution.succ[user] = user+n_cus
    solution.succ[user+n_cus] = end_depot
    solution.pre[end_depot]=user+n_cus
    solution.pre[user+n_cus] = user
    solution.pre[user] = start_depot
    
end

function insert_user_route(solution::Solution, route_i, best_pos, start_depot, end_depot, qi) # user is pickup node is from 2,...,n+1
    
    pre_vi =best_pos[1];user=best_pos[2];suc_vi=best_pos[3]
    pre_v_ni = best_pos[4]; v_ni =best_pos[5];suc_v_ni=best_pos[6]

    solution.RI[route_i,3] += 2 # pickup and dropoff nodes 
    if pre_vi == start_depot
        solution.RI[route_i,1] = user
    end
    if  suc_v_ni  == end_depot
        solution.RI[route_i,2] = v_ni
    end
    solution.succ[pre_vi]=user
    solution.succ[user]=suc_vi
    solution.succ[pre_v_ni] = v_ni
    solution.succ[v_ni] = suc_v_ni
    
    solution.pre[user]= pre_vi
    solution.pre[suc_vi]= user
    solution.pre[v_ni] = pre_v_ni
    solution.pre[suc_v_ni] =  v_ni
    
end

# correct e_r,cap_r,l_r,maxcap_r at the depot nodes
function init_route_info(e_r,cap_r,l_r,maxcap_r, li, start_depot, end_depot, dt_r, dt_r_back)
    
    e_r[start_depot]=0; cap_r[start_depot]=0
    l_r[end_depot] = li[end_depot]; maxcap_r[end_depot] = 0
    # dt_r[start_depot] = 0; dt_r_back[end_depot] = 0

end

# if veh_type is infeasibe for route r, the cost = bigM
function check_cap_r(solution::Solution, darp, r, parameter_energy, type_veh, matrix_cost_v_type)

    bigM = 1e6
    start_depot = darp.start_depot
    n_node = 2 + 2 * darp.n_cus
    cap_r_tmp = zeros(n_node)
    cap_type_veh =  parameter_energy.veh_info.cap_veh_type[type_veh]
    v_i = start_depot;  v_j = solution.RI[r,1] 
    while v_j > 0 
        cap_r_tmp[v_j] =  cap_r_tmp[v_i] + darp.qi[v_j]
        v_i,v_j = v_j, solution.succ[v_j]
    end
    flag= findfirst(x->x>cap_type_veh, cap_r_tmp)
    if !isnothing(flag)
        if solution.n_route ==1
            matrix_cost_v_type[type_veh] = bigM
        else
            matrix_cost_v_type[r, type_veh] = bigM
        end
    end
end

function init_e_r_cap_r(solution::Solution, darp, route_i, start_depot, e_r, cap_r, ei)

    e_r[start_depot]=0; cap_r[start_depot]=0
    v_i = start_depot
    v_j = solution.RI[route_i,1]
    while v_j > 0
        e_r[v_j] = max(ei[v_j], e_r[v_i]+darp.s_t[v_i]+darp.dist[v_i,v_j])
        cap_r[v_j] =  cap_r[v_i] + darp.qi[v_j]
        v_i,v_j = v_j, solution.succ[v_j]
    end
end


# get the least cost and position to insert a request on a route
function cost_greedy_insert(solution::Solution, fast_chg, global_parameter, instance, darp, user, route_i, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back, flag_init_sol)
      
    bigM = global_parameter.bigM
    Q, parameter_energy, K = instance.Q_max, instance.parameter_energy, instance.K
    ei, li, s_t, Li = darp.ei, darp.li, darp.s_t, darp.Li
    n_cus, TH, start_depot, end_depot = darp.n_cus, darp.TH, darp.start_depot, darp.end_depot
    qi, dist, dist_orig = darp.qi, darp.dist, darp.dist_orig
    dist_all = darp.dist_all 

    init_route_info(e_r,cap_r,l_r,maxcap_r, li, start_depot, end_depot, dt_r,dt_r_back)
    init_e_r_cap_r(solution, darp, route_i, start_depot, e_r, cap_r, ei)
    inserted = false
    min_dist_increase = bigM 
    best_pos = Int32[]
    route = get_route(solution, route_i, start_depot, end_depot)
    # route = route_with_depots[2:end-1]   
    n_nodes = solution.RI[route_i,3]  
    vi    = user
    v_ni = user + n_cus
    # @show(route_i, route_with_depots, n_nodes, vi, v_ni)
    # insert pickup node, idx and idx_2 are i_th insert positions after the depot 
    idx=1 
    pre_vi = start_depot 
    suc_vi = solution.RI[route_i,1] 
    suc_v_ni = solution.RI[route_i,1]  
    copy_route = Int32[]; best_route_temp = Int32[]      
    while idx < n_nodes +2
        #check tw and cap at user's pickup location
        eh_i = max(ei[vi], e_r[pre_vi]+s_t[pre_vi]+dist[pre_vi,vi])
        if  eh_i <= li[vi] && cap_r[pre_vi]+qi[vi] <= Q[route_i] && qi[vi] + maxcap_r[suc_vi] <= Q[route_i]
            # insert v,v_ni immidiately after the depot
            dist_increase = dist[pre_vi, user] + dist[user, v_ni] + dist[v_ni, suc_vi]
                - dist[pre_vi,suc_vi] 
            if dist_increase < min_dist_increase
                eh_v_ni = max(ei[v_ni],eh_i+s_t[vi]+dist[vi,v_ni])
                eh_suc_v_ni = max(ei[suc_vi],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_vi])
                if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
                    if eh_v_ni - li[user] - s_t[user] <= Li[user] 
                        copy_route = copy(route)
                        insert!(copy_route, idx+1, user); insert!(copy_route, idx+2 , user + n_cus)
                        if  eight_step(copy_route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi) 
                            min_dist_increase = dist_increase
                            best_pos = [pre_vi,vi,v_ni, vi, v_ni,suc_vi]
                            best_route_temp = copy(copy_route)
                            inserted = true
                        end
                    end
                end
            end
            # insert vi (user), v_ni after the first request
            idx_2=idx+1
            eh_pre_v_ni = max(ei[suc_vi], eh_i+s_t[vi]+dist[vi,suc_vi])
            pre_v_ni = suc_vi
            suc_v_ni = solution.succ[pre_v_ni]
            while idx_2 < n_nodes +2 
                eh_v_ni = max(ei[v_ni],eh_pre_v_ni+s_t[pre_v_ni]+dist[pre_v_ni,v_ni])
                eh_suc_v_ni = max(ei[suc_v_ni],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_v_ni])
                eh_pre_v_ni = eh_v_ni
                dist_increase = dist[pre_vi, user] + dist[user, suc_vi] + dist[pre_v_ni, v_ni] + dist[v_ni, suc_v_ni]
                              - ( dist[pre_vi, suc_vi] + dist[pre_v_ni, suc_v_ni] )

                if dist_increase < min_dist_increase 
                    if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
                        if eh_v_ni - li[user] - s_t[user] <= Li[user]  
                            copy_route = copy(route)
                            insert!(copy_route, idx+1, user); insert!(copy_route, idx_2+2 , user + n_cus)
                            if  eight_step(copy_route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi) 
                                # check ride_time and route duration
                                min_dist_increase = dist_increase
                                best_pos = [pre_vi,vi,suc_vi, pre_v_ni, v_ni,suc_v_ni]
                                best_route_temp = copy(copy_route) 
                                inserted = true
                            end
                        end                   
                    end
                end                
                idx_2 += 1 
                pre_v_ni = suc_v_ni
                suc_v_ni = solution.succ[suc_v_ni] 
            end
        end
        idx += 1
        pre_vi = suc_vi
        suc_vi = solution.succ[suc_vi]
        suc_v_ni = suc_vi
    end 
    return min_dist_increase, best_pos, best_route_temp
     
end

# insert the user to the best positions of route_i
function greedy_insert(solution::Solution, fast_chg, global_parameter, instance, darp, user, route_i, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back, flag_init_sol)
      
    bigM = global_parameter.bigM
    Q, parameter_energy, K, T_max = instance.Q_max, instance.parameter_energy, instance.K, instance.T_max
    ei, li, s_t, Li = darp.ei, darp.li, darp.s_t, darp.Li
    n_cus, TH, start_depot, end_depot = darp.n_cus, darp.TH, darp.start_depot, darp.end_depot
    qi, dist, dist_orig = darp.qi, darp.dist, darp.dist_orig
    dist_all, set_physical = darp.dist_all, darp.set_physical

    init_route_info(e_r,cap_r,l_r,maxcap_r, li, start_depot, end_depot, dt_r,dt_r_back)
    init_e_r_cap_r(solution, darp, route_i, start_depot, e_r, cap_r, ei)
    inserted = false; success_chg=false
    min_dist_increase = bigM 
    best_pos = Int32[]
    route = get_route(solution, route_i, start_depot, end_depot)
    # route = route_with_depots[2:end-1]   
    n_nodes = solution.RI[route_i,3]  
    vi    = user
    v_ni = user + n_cus
    # @show(route_i, route_with_depots, n_nodes, vi, v_ni)
    # insert pickup node, idx and idx_2 are i_th insert positions after the depot 
    idx=1 
    pre_vi = start_depot 
    suc_vi = solution.RI[route_i,1] 
    suc_v_ni = solution.RI[route_i,1]  
    copy_route = Int32[]; best_route_temp = Int32[]      
    while idx < n_nodes +2
        #check tw and cap at user's pickup location
        eh_i = max(ei[vi], e_r[pre_vi]+s_t[pre_vi]+dist[pre_vi,vi])
        if  eh_i <= li[vi] && cap_r[pre_vi]+qi[vi] <= Q[route_i] && qi[vi] + maxcap_r[suc_vi] <= Q[route_i]
            # insert v,v_ni immidiately after the depot
            dist_increase = dist[pre_vi, user] + dist[user, v_ni] + dist[v_ni, suc_vi]
                - dist[pre_vi,suc_vi] 
            if dist_increase < min_dist_increase
                eh_v_ni = max(ei[v_ni],eh_i+s_t[vi]+dist[vi,v_ni])
                eh_suc_v_ni = max(ei[suc_vi],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_vi])
                if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
                    if eh_v_ni - li[user] - s_t[user] <= Li[user] 
                        copy_route = copy(route)
                        insert!(copy_route, idx+1, user); insert!(copy_route, idx+2 , user + n_cus)
                        if  eight_step(copy_route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi) 
                            min_dist_increase = dist_increase
                            best_pos = [pre_vi,vi,v_ni, vi, v_ni,suc_vi]
                            best_route_temp = copy(copy_route)
                            inserted = true
                        end
                    end
                end
            end
            # insert vi (user), v_ni after the first request
            idx_2=idx+1
            eh_pre_v_ni = max(ei[suc_vi], eh_i+s_t[vi]+dist[vi,suc_vi])
            pre_v_ni = suc_vi
            suc_v_ni = solution.succ[pre_v_ni]
            while idx_2 < n_nodes +2 
                eh_v_ni = max(ei[v_ni],eh_pre_v_ni+s_t[pre_v_ni]+dist[pre_v_ni,v_ni])
                eh_suc_v_ni = max(ei[suc_v_ni],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_v_ni])
                eh_pre_v_ni = eh_v_ni
                dist_increase = dist[pre_vi, user] + dist[user, suc_vi] + dist[pre_v_ni, v_ni] + dist[v_ni, suc_v_ni]
                              - ( dist[pre_vi, suc_vi] + dist[pre_v_ni, suc_v_ni] )

                if dist_increase < min_dist_increase 
                    if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
                        if eh_v_ni - li[user] - s_t[user] <= Li[user]  
                            copy_route = copy(route)
                            insert!(copy_route, idx+1, user); insert!(copy_route, idx_2+2 , user + n_cus)
                            if  eight_step(copy_route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi) 
                                # check ride_time and route duration
                                min_dist_increase = dist_increase
                                best_pos = [pre_vi,vi,suc_vi, pre_v_ni, v_ni,suc_v_ni]
                                best_route_temp = copy(copy_route) 
                                inserted = true
                            end
                        end                       
                    end
                end                
                idx_2 += 1 
                pre_v_ni = suc_v_ni
                suc_v_ni = solution.succ[suc_v_ni] 
            end
        end
        idx += 1
        pre_vi = suc_vi
        suc_vi = solution.succ[suc_vi]
        suc_v_ni = suc_vi
    end
    if inserted == true
        # route = get_route(soluT, route_i, start_depot, end_depot)
    #     n_charger_installed, success_chg, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(solution, darp, route_i, best_route_temp, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
    #     if ! success_chg && (n_charger_installed > 0)                            
    #         success_chg = repair_charging(solution, darp, route_i, best_route_temp, dist_all, fast_chg, parameter_energy, set_physical, dt_r, pos_insert_chg, forward_time_slack, vec_t_start)
    #     end
    #     inserted = success_chg
    # end
        # @show(route_i, route, best_route_temp, ec_route, parameter_energy.max_ec[route_i], length_route(best_route_temp)*parameter_energy.beta[route_i] )
        veh_id = solution.RI[route_i,5]
        if parameter_energy.is_electric[veh_id] == true # EV 
            compute_dt_r(best_route_temp, dt_r, dist_orig)
            ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
            if ec_route > parameter_energy.max_ec[veh_id]
                # @show(veh_id, ec_route, parameter_energy.max_ec[veh_id], get_route(solution, route_i, start_depot, end_depot))
                if fast_chg.n_fast_chg_installed == 0
                    return false
                else
                    if insert_charging_route(solution, darp, route_i, best_route_temp, dist_all, ec_route, fast_chg, parameter_energy, dt_r, dt_r_back, flag_init_sol)[1]
                        insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
                        update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
                    else
                        inserted = false
                    end 
                end
            else
                solution.vec_chg_event[route_i,:] .*= 0           
                insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
                update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
            end             
        else
            solution.vec_chg_event[route_i,:] .*= 0           
            insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
            update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
        end         
    end
    return inserted
end


function gen_solution(sol::Solution, solu_init::Solution, global_parameter, fast_chg, darp, instance)

    # sol= solution
    Q = instance.Q_max 
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    ei, li, qi  =  darp.ei, darp.li, darp.qi 
    N_NODE = 2*n_cus +2
    penalty  =  global_parameter.penalty
    parameter_energy = instance.parameter_energy

    flag_init_sol = true
    sol.succ .= -1
    sol.pre  .= -1
    sol.RI .*=  0 
    sol.RI[:,4] = Q  # Q is a vector of capacity of veh
    sol.RI[:, 5] = copy(solu_init.RI[:,5])# veh_id
    sol.unserved_users = Int32[]
    sol.penalty_unserved = 0
    sol.succ[start_depot]=end_depot
    sol.pre[end_depot]=start_depot
    sol.n_route = 0; sol.total_cost=0; sol.total_cost_with_penalty=0; sol.total_chg_time =0;
    sol.vec_chg_event .*= 0 
    sol.co2_emission = 0
    
    # init auxiliary info  
    e_r   = copy(ei)
    cap_r = zeros(Int32,N_NODE)
    l_r   = copy(li)
    maxcap_r = zeros(Int32,N_NODE)
    dt_r  = zeros(Float32,N_NODE)  
    dt_r_back = zeros(Float32,N_NODE)  
    
    e_r[1]=ei[1]; e_r[end_depot]=ei[end_depot]
    l_r[1]=li[1]; l_r[end_depot]=li[end_depot]

    unserved = collect(2:n_cus+1) 
    shuffle!(unserved)  
    # insert first request
    user = unserved[1]
    sol.n_route +=1  
    insert_consecutive(sol, 1, start_depot, user, end_depot, n_cus, qi) 
    update_auxiliary_var(sol, darp, 1, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem 
    setdiff!(unserved, user)
    inserted = false
    for user in unserved
        n_route = sol.n_route 
        route_rand= shuffle!(collect(1:n_route))  
        for r in route_rand
            inserted = greedy_insert(sol, fast_chg, global_parameter, instance, darp, user, r, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back, flag_init_sol)
            inserted && break
        end  
        if !inserted # insert the user on a new route
            sol.n_route +=1  
            insert_consecutive(sol, sol.n_route, start_depot, user, end_depot, n_cus, qi) 
            update_auxiliary_var(sol, darp, sol.n_route, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem 
        end        
    end
    if  sol.n_route > instance.K_MAX
        error("sol.n_route > instance.K_MAX error ! in gen_solution()")
    end
    # set_unserved_users(sol, qi, penalty, K, n_cus)
    cost_solution(sol, global_parameter, darp, parameter_energy, fast_chg, instance) # it does not matter for neglecting the computation of waiting time for the init solution 
    
    return sol.total_cost_with_penalty, sol, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back
end


function gen_init_solution(solution::Solution, solu_init::Solution, global_parameter, fast_chg, darp, instance, n_init_sol, Q, N_NODE, start_depot, end_depot, ei, li, K_MAX, n_cus, qi, Li, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval)
    
    # solution =solu
    soluT = deepcopy(solu_init)
    e_r   = copy(ei)
    cap_r = zeros(Int32,N_NODE)
    l_r   = copy(li)
    maxcap_r = zeros(Int32,N_NODE)
    dt_r  = zeros(Float32,N_NODE) # travel time from the start_depot to node i on route r
    dt_r_back = zeros(Float32,N_NODE) # travel time from the end_depot back to node i
    flag_success = false
    min_val = Inf
    for i= 1:n_init_sol   
        val_current, sol_current, e_r_sol,
        cap_r_sol, l_r_sol, maxcap_r_sol, dt_r_sol, dt_r_back_sol = gen_solution(solution, solu_init, global_parameter, fast_chg, darp, instance)
        if  val_current < min_val && check_chg_occ_constr(sol_current, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval) 
            min_val = val_current             
            update_sol(soluT, sol_current)
            e_r, cap_r, l_r, maxcap_r = copy(e_r_sol), copy(cap_r_sol), copy(l_r_sol), copy(maxcap_r_sol)
            dt_r, dt_r_back = copy(dt_r_sol), copy(dt_r_back_sol)
            flag_success = true             
        end
    end    
    # n_route= soluT.n_route; veh_ids= soluT.RI[1:n_route,5]
    # @show(flag_success, n_route, instance.parameter_energy.veh_type[veh_ids], soluT.co2_emission)
    if flag_success 
       update_sol(solution, soluT)
    else 
        repair_init_sol(solution, solu_init, n_init_sol, global_parameter, instance, darp, fast_chg, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
    end
    return e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back 
end


# gen ev routes no charging while energy feasible
function repair_init_sol(solution::Solution, solu_init::Solution, n_init_sol, global_parameter, instance, darp, fast_chg, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)

    # setting n_charger as 0 temporarly
    copy_n_fast_chg_installed = fast_chg.n_fast_chg_installed
    fast_chg.n_fast_chg_installed = 0
    flag_success = false
    min_val = Inf
    soluT = deepcopy(solu_init) 

    for i= 1:n_init_sol   
        val_current, sol_current, e_r_sol,
        cap_r_sol, l_r_sol, maxcap_r_sol, dt_r_sol, dt_r_back_sol = gen_solution(solution, solu_init, global_parameter, fast_chg, darp, instance)
        if  val_current < min_val && check_chg_occ_constr(sol_current, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval) 
            min_val = val_current  
            update_sol(soluT, sol_current)
            e_r, cap_r, l_r, maxcap_r = copy(e_r_sol), copy(cap_r_sol), copy(l_r_sol), copy(maxcap_r_sol)
            dt_r, dt_r_back = copy(dt_r_sol), copy(dt_r_back_sol)
            flag_success = true             
        end
    end    
    if flag_success 
        update_sol(solution, soluT)
    else
        error("repair_init_sol failed, please check !!")
    end
    fast_chg.n_fast_chg_installed =  copy_n_fast_chg_installed # restore it

end 


# # function to compute the threshold to switch from  EV to gasoline vehicle based on the cost savings
# # given a travel distance (time), when the travel cost of a ev is smaller than the threshold, exchange it with a gasoline veh would save costs
# function get_threshold_advantage_ev(parameter_energy, v_k)

#     n_type = parameter_energy.veh_info.n_type_veh 
#     big_value = 10000
#     vec_thresholds = zeros(n_type)
#     veh_gasoline = findfirst(x->x==0, parameter_energy.veh_info.is_electric)
#     for type_veh in 1:n_type      
#         if  type_veh == veh_gasoline
#              vec_thresholds[type_veh] = -1
#         else 
#             diff_gasoline_veh_cost = parameter_energy.veh_info.daily_purchase_cost[veh_gasoline] - parameter_energy.veh_info.daily_purchase_cost[type_veh]
#             diff_gasoline_cost     =  parameter_energy.veh_info.cost_energy[veh_gasoline] * parameter_energy.veh_info.energy_consumption_km[veh_gasoline]
#                                         - parameter_energy.veh_info.cost_energy[type_veh] * parameter_energy.veh_info.energy_consumption_km[type_veh]
#             if diff_gasoline_veh_cost * diff_gasoline_cost < 0
#                 threshold_dist = diff_gasoline_veh_cost / (- diff_gasoline_cost)
#                 vec_thresholds[type_veh] = threshold_dist/v_k
#             end
#             if diff_gasoline_veh_cost > 0 && diff_gasoline_cost > 0
#                 vec_thresholds[type_veh] = big_value
#             end  
#         end
#     end
#     return vec_thresholds
# end


# #function to check whether the charged energy is sufficient for the route
# function get_add_tt_chg(solution::Solution, darp, r, fast_chg, parameter_energy)
 
#     set_physical, dist_all = darp.set_physical, darp.dist_all
#     start_depot, end_depot = darp.start_depot, darp.end_depot
#     add_tt = 0 ; max_num_chg = parameter_energy.max_num_chg
#     info_chg = zeros(Float32, 2+max_num_chg*4) # store at most 4 charging operations
#     route = get_route(solution, r, start_depot, end_depot)
#     for i in 1:floor(Int32, solution.vec_chg_event[r,2])
#         v1    = floor(Int32, solution.vec_chg_event[r, 2+(i-1)*4+1])
#         idx = findfirst(x->x==v1, route)
#         if isnothing(idx) # charging schedule walks off, need to reschedule it
#            return false, 0
#         end 
#         v2= route[idx+1]
#         v1, v2   = set_physical[v1], set_physical[v2]
#         chg_id   = floor(Int32,solution.vec_chg_event[r, 2+(i-1)*4+2]) # chg id on the list of used fast chargers
#         v_physi_chg       =  fast_chg.v_physical_chg[chg_id]
#         t_access_1 = dist_all[v1, v_physi_chg] 
#         t_access_2 = dist_all[v_physi_chg, v2]
#         dist_v1_v2 = dist_all[v1, v2]
#         add_tt += (t_access_1 + t_access_2 - dist_v1_v2)
#     end
#     return true, add_tt
# end
 


# # rand insert user in a feasible and acceptable position of a route or in the unserved user pool
# function rand_insert(solution::Solution, fast_chg, global_parameter, instance, darp, user, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, start_depot, end_depot, Li)
    
#     TH, n_cus  = darp.TH, darp.n_cus
#     ei, li, s_t, qi, dist, dist_orig = darp.ei, darp.li, darp.s_t, darp.qi, darp.dist, darp.dist_orig
#      dist_all =  darp.dist_all   
#     parameter_energy, K = instance.parameter_energy, instance.K
#     bigM = global_parameter.bigM
#     flag_init_sol = false
#     # insert a user to the best positions
#     init_route_info(e_r,cap_r,l_r,maxcap_r, li, start_depot, end_depot, dt_r, dt_r_back)
#     init_e_r_cap_r(solution, darp, route_i, start_depot, e_r, cap_r, ei)
#     inserted = false
#     min_dist_increase = bigM 
#     best_pos = Int32[] 
#     n_nodes = solution.RI[route_i,3] 
#     route = get_route(solution, route_i, start_depot, end_depot)
#     # @show(user,route_i,route)
#     vi    = user
#     v_ni = user + n_cus
#     # insert pickup node, idx and idx_2 are i_th insert positions after the depot 
#     idx=1 
#     pre_vi = start_depot 
#     suc_vi = solution.RI[route_i,1] 
#     suc_v_ni = solution.RI[route_i,1]  
#     copy_route = Int32[]; best_route_temp = Int32[]
#     # copy_route =  zeros(Int32, length(route)+2)

#     while idx < n_nodes +2
#         #check tw and cap at user's pickup location
#         eh_i = max(ei[vi], e_r[pre_vi]+s_t[pre_vi]+dist[pre_vi,vi])
#         # @show(Q, Q[route_i])
#         if  eh_i <= li[vi] && cap_r[pre_vi]+qi[vi] <= Q[route_i] && qi[vi] + maxcap_r[suc_vi] <= Q[route_i]
#             dist_increase = dist[pre_vi, user] + dist[user, v_ni] + dist[v_ni, suc_vi]
#                 - dist[pre_vi,suc_vi] 
#             if dist_increase < min_dist_increase
#                 eh_v_ni = max(ei[v_ni],eh_i+s_t[vi]+dist[vi,v_ni])
#                 eh_suc_v_ni = max(ei[suc_vi],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_vi])
#                 if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
#                     if eh_v_ni - li[user] - s_t[user] <= Li[user]  
#                         copy_route = copy(route)
#                         insert!(copy_route, idx+1, user); insert!(copy_route, idx+2 , user + n_cus)
#                         if eight_step(copy_route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)
#                             min_dist_increase = dist_increase
#                             best_pos = [pre_vi,vi,v_ni, vi, v_ni,suc_vi]
#                             best_route_temp = copy(copy_route)
#                             inserted = true
#                         end
#                     end
#                 end
#             end
#             # insert dropoff later
#             idx_2=idx+1
#             eh_pre_v_ni = max(ei[suc_vi], eh_i+s_t[vi]+dist[vi,suc_vi])
#             pre_v_ni = suc_vi
#             suc_v_ni = solution.succ[pre_v_ni]
#             while idx_2 < n_nodes +2
#                 # check tw for pre_v_ni if idx_2 > idx
#                 # @show(idx_2,pre_v_ni, v_ni, suc_v_ni)
#                 eh_v_ni = max(ei[v_ni],eh_pre_v_ni+s_t[pre_v_ni]+dist[pre_v_ni,v_ni])
#                 eh_suc_v_ni = max(ei[suc_v_ni],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_v_ni])
#                 eh_pre_v_ni = eh_v_ni
#                 dist_increase = dist[pre_vi, user] + dist[user, suc_vi] + dist[pre_v_ni, v_ni] + dist[v_ni, suc_v_ni]
#                               - ( dist[pre_vi, suc_vi] + dist[pre_v_ni, suc_v_ni] )

#                 if dist_increase < min_dist_increase 
#                     if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
#                         if eh_v_ni - li[user] - s_t[user] <= Li[user]  
#                             copy_route = copy(route)
#                             insert!(copy_route, idx+1, user); insert!(copy_route, idx_2+2 , user + n_cus)
#                             if eight_step(copy_route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)
#                                 # check ride_time and route duration
#                                 min_dist_increase = dist_increase
#                                 best_pos = [pre_vi,vi,suc_vi, pre_v_ni, v_ni,suc_v_ni]
#                                 best_route_temp = copy(copy_route)
#                                 inserted = true
#                                 break
#                             end
#                         end                     
#                     end
#                 end                
#                 idx_2 += 1 
#                 pre_v_ni = suc_v_ni; 
#                 suc_v_ni = solution.succ[suc_v_ni] 
#             end
#         end
#         idx += 1
#         pre_vi = suc_vi
#         suc_vi = solution.succ[suc_vi]
#         suc_v_ni = suc_vi
#     end
#     if inserted == true
#         veh_id = solution.RI[route_i,5]
#         if parameter_energy.is_electric[veh_id] == true # EV 
#             compute_dt_r(best_route_temp, dt_r, dist_orig)
#             ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
#             if ec_route > parameter_energy.max_ec[veh_id]
#                 if fast_chg.n_fast_chg_installed == 0
#                     return false
#                 else
#                     if insert_charging_route(solution, darp, route_i, best_route_temp, dist_all, ec_route, fast_chg, parameter_energy, dt_r, dt_r_back, flag_init_sol)[1]
#                         insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
#                         update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#                     else
#                         inserted = false
#                     end    
#                 end        
#             else
#                 solution.vec_chg_event[route_i,:] .*= 0           
#                 insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
#                 update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) 
#             end
#         else
#             insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
#             update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#         end
#     end
#     return inserted
# end

# function repair_init_sol(solution::Solution, sol_current_best::Solution, instance, darp, fast_chg, occ_state_chg_init, discrte_t_interval, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
    
#     # n_route = sol_current_best.n_route
#     # start_depot,end_depot = darp.start_depot,darp.end_depot
#     # for r in 1:n_route
#     #     @show(get_route(sol_current_best,r, start_depot,end_depot))
#     # end
#     # verify_solution(sol_current_best, instance, darp, fast_chg)
#     flag_silent = false
#     @show("repair_init_sol call !!")
#     success= false; iter_big_enough = 50
#     for count_split_route in 1:iter_big_enough
#         @show(count_split_route, verify_solution_silent(sol_current_best, instance, darp, fast_chg) )
#         split_route(sol_current_best, instance, darp, fast_chg, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, occ_state_chg_init, discrte_t_interval)
#         if  check_chg_occ_constr(sol_current_best, occ_state_chg_init, discrte_t_interval) 
#             if ! verify_solution_silent(sol_current_best, instance, darp, fast_chg) 
#                error("verify_solution_silent failed !! in repair_init_sol(), please check !! ")
#             else
#                 update_sol(solution, sol_current_best) 
#                 for r in 1:solution.n_route
#                     reset_auxiliary_route(solution, darp, r, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
#                 end
#                 success= true ;break
#             end
#         end
#     end
#     return success
# end 

# function set_unserved_users(solution::Solution, qi, penalty, K, n_cus)
    
#     if solution.n_route >K # move all users on unserved pool 
#         for r in K+1:solution.n_route 
#             route =get_route_no_depot(solution, r)
#             n_nodes = length(route)
#             solution.pre[route] = ones(Int32,n_nodes)*-1
#             solution.succ[route] = ones(Int32,n_nodes)*-1
#             solution.RI[r, 1:3]= zeros(Int32,3) 
#             users = route[findall(x->x<n_cus+2, route)]
#             union!(solution.unserved_users, users)
#             solution.penalty_unserved += sum(qi[users]) * penalty
#         end
#         solution.n_route = K
#     end
# end


# function gen_solution(sol::Solution, global_parameter, fast_chg, darp, instance, Q, start_depot,end_depot,ei,li,N_NODE,K_MAX,n_cus, qi, Li, occ_state_chg, discrte_t_interval)

#     # sol = solu
#     K =instance.K
#     penalty  =  global_parameter.penalty
#     s_t, dist= darp.s_t, darp.dist
#     parameter_energy = instance.parameter_energy

#     flag_init_sol = true
#     sol.succ .= -1
#     sol.pre  .= -1
#     sol.RI .*=  0 
#     sol.RI[:,4] = Q  # Q is a vector of capacity of veh
#     sol.RI[:, 5] = collect(1:K_MAX)# veh_id
#     sol.unserved_users = Set{Int32}()
#     sol.penalty_unserved = 0
#     sol.succ[start_depot]=end_depot
#     sol.pre[end_depot]=start_depot
#     sol.n_route = 0; sol.total_cost=0; sol.total_chg_time =0;
#     sol.vec_chg_event .*= 0 
#     sol.co2_emission = 0
    
#     # init auxiliary info  
#     e_r   = copy(ei)
#     cap_r = zeros(Int32,N_NODE)
#     l_r   = copy(li)
#     maxcap_r = zeros(Int32,N_NODE)
#     dt_r  = zeros(Float32,N_NODE)  
#     dt_r_back = zeros(Float32,N_NODE)  
    
#     e_r[1]=ei[1]; e_r[end_depot]=ei[end_depot]
#     l_r[1]=li[1]; l_r[end_depot]=li[end_depot]

#     unserved = collect(2:n_cus+1) 
#     shuffle!(unserved) 

#     n_min = min(n_cus,K)
#     for r in 1:n_min
#         for user in unserved
#             route = [start_depot, user, user+n_cus, end_depot] 
#             if  check_tw(route, 0, ei, li, s_t, dist)
#                 sol.n_route += 1
#                 insert_consecutive(sol, r, start_depot, user, end_depot, n_cus, qi) 
#                 update_auxiliary_var(sol, darp, r, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem 
#                 setdiff!(unserved, user)
#                 break
#             end                
#         end
#     end
#     if K < n_cus
#         for user in unserved     
#             inserted = false   
#             for r in 1:K 
#                 inserted = greedy_insert(sol, fast_chg, global_parameter, instance, darp, user, r, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back, flag_init_sol)
#                 inserted && break
#             end
#             if ! inserted  
#                 sol.pre[user]= -1; sol.succ[user]= -1  
#                 sol.pre[user+n_cus]= -1; sol.succ[user+n_cus]= -1   
#                 union!(sol.unserved_users, user)
#                 sol.penalty_unserved += sum(qi[user]) * penalty
#             end
#         end
#     end 
#     set_unserved_users(sol, qi, penalty, K, n_cus)
#     sol.total_cost = cost_solution(sol, global_parameter, darp, parameter_energy, penalty) # it does not matter for neglecting the computation of waiting time for the init solution 
#     return sol.total_cost, sol, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back
# end

# insert the user to the best positions of route_i
# function greedy_insert(solution::Solution, fast_chg, global_parameter, instance, darp, user, route_i, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back, flag_init_sol)
      
#     bigM = global_parameter.bigM
#     Q, parameter_energy, K = instance.Q_max, instance.parameter_energy, instance.K
#     ei, li, s_t, Li = darp.ei, darp.li, darp.s_t, darp.Li
#     n_cus, TH, start_depot, end_depot = darp.n_cus, darp.TH, darp.start_depot, darp.end_depot
#     qi, dist, dist_orig = darp.qi, darp.dist, darp.dist_orig
#     dist_all = darp.dist_all 

#     init_route_info(e_r,cap_r,l_r,maxcap_r, li, start_depot, end_depot, dt_r,dt_r_back)
    # init_e_r_cap_r(solution, darp, route_i, start_depot, e_r, cap_r, ei)
#     inserted = false
#     min_dist_increase = bigM *0.7
#     best_pos = ()    
#     route = get_route(solution, route_i, start_depot, end_depot)
#     # route = route_with_depots[2:end-1]   
#     n_nodes = solution.RI[route_i,3]  
#     vi    = user
#     v_ni = user + n_cus
#     # @show(route_i, route_with_depots, n_nodes, vi, v_ni)
#     # insert pickup node, idx and idx_2 are i_th insert positions after the depot 
#     idx=1 
#     pre_vi = start_depot 
#     suc_vi = solution.RI[route_i,1] 
#     suc_v_ni = solution.RI[route_i,1]  
#     copy_route = Int32[]; best_route_temp = Int32[]      
#     while idx < n_nodes +2
#         #check tw and cap at user's pickup location
#         eh_i = max(ei[vi], e_r[pre_vi]+s_t[pre_vi]+dist[pre_vi,vi])
#         if  eh_i <= li[vi] && cap_r[pre_vi]+qi[vi] <= Q[route_i] && qi[vi] + maxcap_r[suc_vi] <= Q[route_i]
#             dist_increase = dist[pre_vi, user] + dist[user, v_ni] + dist[v_ni, suc_vi]
#                 - dist[pre_vi,suc_vi] 
#             if dist_increase < min_dist_increase
#                 eh_v_ni = max(ei[v_ni],eh_i+s_t[vi]+dist[vi,v_ni])
#                 eh_suc_v_ni = max(ei[suc_vi],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_vi])
#                 if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
#                     if eh_v_ni - li[user] - s_t[user] <= Li[user] 
#                         copy_route = copy(route)
#                         insert!(copy_route, idx+1, user); insert!(copy_route, idx+2 , user + n_cus)
#                         if  eight_step(copy_route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi) 
#                             min_dist_increase = dist_increase
#                             best_pos = (pre_vi,vi,v_ni, vi, v_ni,suc_vi) 
#                             best_route_temp = copy(copy_route)
#                             inserted = true
#                         end
#                     end
#                 end
#             end
#             # insert dropoff later
#             idx_2=idx+1
#             eh_pre_v_ni = max(ei[suc_vi], eh_i+s_t[vi]+dist[vi,suc_vi])
#             pre_v_ni = suc_vi
#             suc_v_ni = solution.succ[pre_v_ni]
#             while idx_2 < n_nodes +2 
#                 eh_v_ni = max(ei[v_ni],eh_pre_v_ni+s_t[pre_v_ni]+dist[pre_v_ni,v_ni])
#                 eh_suc_v_ni = max(ei[suc_v_ni],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_v_ni])
#                 eh_pre_v_ni = eh_v_ni
#                 dist_increase = dist[pre_vi, user] + dist[user, suc_vi] + dist[pre_v_ni, v_ni] + dist[v_ni, suc_v_ni]
#                               - ( dist[pre_vi, suc_vi] + dist[pre_v_ni, suc_v_ni] )

#                 if dist_increase < min_dist_increase 
#                     if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
#                         if eh_v_ni - li[user] - s_t[user] <= Li[user]  
#                             copy_route = copy(route)
#                             insert!(copy_route, idx+1, user); insert!(copy_route, idx_2+2 , user + n_cus)
#                             if  eight_step(copy_route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi) 
#                                 # check ride_time and route duration
#                                 min_dist_increase = dist_increase
#                                 best_pos = (pre_vi,vi,suc_vi, pre_v_ni, v_ni,suc_v_ni)
#                                 best_route_temp = copy(copy_route) 
#                                 inserted = true
#                             end
#                         end
#                     else
#                         break                        
#                     end
#                 end                
#                 idx_2 += 1 
#                 pre_v_ni = suc_v_ni
#                 suc_v_ni = solution.succ[suc_v_ni] 
#             end
#         end
#         idx += 1
#         pre_vi = suc_vi
#         suc_vi = solution.succ[suc_vi]
#         suc_v_ni = suc_vi
#     end
#     if inserted == true
#         # @show(route_i, route, best_route_temp, ec_route, parameter_energy.max_ec[route_i], length_route(best_route_temp)*parameter_energy.beta[route_i] )
#         veh_id = solution.RI[route_i,5]
#         if parameter_energy.is_electric[veh_id] == true # EV 
#             compute_dt_r(best_route_temp, dt_r, dist_orig)
#             ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
#             if ec_route > parameter_energy.max_ec[veh_id]
#                 # @show(veh_id, ec_route, parameter_energy.max_ec[veh_id], get_route(solution, route_i, start_depot, end_depot))
#                 if fast_chg.n_fast_chg_installed == 0
#                     return false
#                 else
#                     if insert_charging_route(solution, darp, route_i, best_route_temp, dist_all, ec_route, fast_chg, parameter_energy, dt_r, dt_r_back, flag_init_sol)[1]
#                         insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
#                         update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#                     else
#                         inserted = false
#                     end 
#                 end
#             else
#                 solution.vec_chg_event[route_i,:] .*= 0           
#                 insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
#                 update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#             end             
#         else
#             solution.vec_chg_event[route_i,:] .*= 0           
#             insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
#             update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#         end         
#     end
#     return inserted
# end


# function check_dt_r(solution::Solution, r, dist_route)
    
#     length_r = length_route(solution, r, 1)
#     if round(length_r - dist_route, digits = DIGITS) == 0
#         return true
#     else
#         return false
#     end

# end
#get energy consumption (ec) of route in O(1)
# function get_ec_route(solution::Solution, r, v, dt_r, dt_r_back, parameter_energy)

#     dist_route = dt_r[v]+ dt_r_back[v]
#     if !check_dt_r(solution, r, dist_route)
#         @show(r, length_route(solution, r, 1), dist_route, v,  dt_r[v], dt_r_back[v])
#         @show(route, dt_r[route], dt_r_back[route])
#         error("check_dt_r failed!")
#     else        
#         return dist_route *  parameter_energy.beta[r] 
#     end
# end


# function gen_solution(sol::Solution, Q, start_depot,end_depot,ei,li,N_NODE,K_MAX,n_cus, qi, Li, occ_state_chg, discrte_t_interval)

#     # sol = solu
#     sol.succ = ones(Int32, 2*n_cus+2)*-1
#     sol.pre =  ones(Int32, 2*n_cus+2)*-1
#     sol.RI =   zeros(Int32, K_MAX, 5) 
#     sol.TT_walk_routes .*= 0
#     sol.unserved_users = Set{Int32}()
#     sol.penalty_unserved = 0
#     sol.RI[:,4] = Q  # Q is a vector of capacity of veh
#     sol.RI[:, 5] = collect(1:K_MAX)# veh_id
#     sol.succ[start_depot]=end_depot
#     sol.pre[end_depot]=start_depot
#     sol.n_route = 0; sol.total_cost=0 
#     sol.vec_chg_event .*= 0

#     # init auxiliary info  
#     e_r   = copy(ei)
#     cap_r = zeros(Int32,N_NODE)
#     l_r   = copy(li)
#     maxcap_r = zeros(Int32,N_NODE)
#     dt_r  = zeros(Float32,N_NODE)  
#     dt_r_back = zeros(Float32,N_NODE)  
    
#     e_r[1]=ei[1]; e_r[end_depot]=ei[end_depot]
#     l_r[1]=li[1]; l_r[end_depot]=li[end_depot]

#     users = collect(2:n_cus+1) 
#     # create the first route with least eng consumption route
#     route_i = 1
#     tmp = [dist_orig[start_depot, i] + dist_orig[i, i+n_cus] + dist_orig[i+n_cus,end_depot] for i in users]
#     user =argmin(tmp)
#     route = [start_depot, user, user+n_cus, end_depot]
#     compute_dt_r(route, dt_r)
#     veh_id = sol.RI[1,5]
#     ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
#     if ec_route > parameter_energy.max_ec[veh_id]
#         if insert_charging_simple(sol, route_i, route, dist_all, ec_route, fast_chg, parameter_energy, set_physical, nodes, dt_r, dt_r_back)
#             sol.n_route +=1 
#             insert_consecutive(sol, route_i, start_depot, user, end_depot, n_cus, qi, walk_tt_bus_stop) 
#             update_auxiliary_var(sol, route_i, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#         else
#             sol.pre[user]=sol.succ[user]=-1
#             sol.pre[user+n_cus]=sol.succ[user+n_cus]=-1 
#             union!(sol.unserved_users, user)
#             sol.penalty_unserved += sum(qi[user]) * penalty
#         end                     
#     else
#         sol.n_route +=1  
#         insert_consecutive(sol, route_i, start_depot, user, end_depot, n_cus, qi, walk_tt_bus_stop) 
#         update_auxiliary_var(sol, route_i, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#     end

#     setdiff!(users, user)
#     unserved = shuffle(users)
  
#     for user in unserved         
#         # insert the user on route_i
#         inserted = greedy_insert(sol, user, route_i, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back, Q, ei, li, start_depot, end_depot, Li, layer_nodes, lyrs_compatible)
#         if ! inserted # insert the user on a new route
#             route_i_tmp = route_i + 1
#             if route_i_tmp < K+1
#                 route = [start_depot, user, user+n_cus, end_depot]
#                 compute_dt_r(route, dt_r)
#                 veh_id = sol.RI[route_i_tmp,5]
#                 ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
#                 if ec_route > parameter_energy.max_ec[veh_id]
#                     if insert_charging_simple(sol, route_i_tmp, route, dist_all, ec_route, fast_chg, parameter_energy, set_physical, nodes, dt_r, dt_r_back)
#                         route_i += 1
#                         sol.n_route +=1  
#                         insert_consecutive(sol, route_i, start_depot, user, end_depot, n_cus, qi, walk_tt_bus_stop) 
#                         update_auxiliary_var(sol, route_i, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#                     else
#                         sol.pre[user]= -1; sol.succ[user]= -1  
#                         sol.pre[user+n_cus]= -1; sol.succ[user+n_cus]= -1   
#                         union!(sol.unserved_users, user)
#                         sol.penalty_unserved += sum(qi[user]) * penalty
#                     end                        
#                 else
#                     route_i += 1
#                     sol.n_route +=1  
#                     insert_consecutive(sol, route_i, start_depot, user, end_depot, n_cus, qi, walk_tt_bus_stop) 
#                     update_auxiliary_var(sol, route_i, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#                 end
#             else
#                 route_i += 1
#                 sol.n_route +=1  
#                 insert_consecutive(sol, route_i, start_depot, user, end_depot, n_cus, qi, walk_tt_bus_stop) 
#                 update_auxiliary_var(sol, route_i, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#             end
#         end
#     end
#     set_unserved_users(sol, qi, penalty)
#     sol.total_cost = length_solution(sol, λ, penalty, start_depot)
#     return sol.total_cost, sol, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back
# end


# # rand insert a user in a feasible position within the extended comptible layers
# function rand_insert_extend(solution::Solution, user, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, start_depot, end_depot, Li, layer_nodes, lyrs_compatible, lyrs_compatible_extend)
    
#     # insert a user to the best positions
#     # init_route_info(e_r,cap_r,l_r,maxcap_r, li, start_depot, end_depot, dt_r, dt_r_back)
#     # init_e_r_cap_r(solution, route_i, start_depot, e_r, cap_r, ei)
#     best_pos = ()  
#     inserted = false
#     min_dist_increase = bigM  
#     n_nodes = solution.RI[route_i,3] 
#     route_with_depots = get_route(solution, route_i, start_depot, end_depot)
#     vi    = user
#     v_ni  = user + n_cus    
#     idx=1 
#     pre_vi = start_depot 
#     suc_vi = solution.RI[route_i,1] 
#     suc_v_ni = solution.RI[route_i,1]  
#     copy_route = Int32[]; nodes_route = Int32[]; best_route_temp = Int32[]
#     while idx < n_nodes +2
#         eh_i = max(ei[vi], e_r[pre_vi]+s_t[pre_vi]+dist[pre_vi,vi])
#         if eh_i <= li[vi] && cap_r[pre_vi]+qi[vi] <= Q[route_i] && qi[vi] + maxcap_r[suc_vi] <= Q[route_i]
#             dist_increase = dist[pre_vi, user] + dist[user, v_ni] + dist[v_ni, suc_vi]
#             - dist[pre_vi,suc_vi] 
#             if dist_increase < min_dist_increase
#                 eh_v_ni = max(ei[v_ni],eh_i+s_t[vi]+dist[vi,v_ni])
#                 eh_suc_v_ni = max(ei[suc_vi],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_vi])
#                 if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
#                     if eh_v_ni - li[user] - s_t[user] <= Li[user] 
#                         # copy_route[1:idx-1], copy_route[idx+2:end] = route_with_depots[1:idx-1], route_with_depots[idx:end]
#                         # copy_route[idx:idx+1] = [user, user+n_cus]
#                         copy_route = copy(route_with_depots)
#                         insert!(copy_route, idx+1, user); insert!(copy_route, idx+2 , user + n_cus)
#                         nodes_route = nodes[copy_route[2:end-1]] 
#                         # if  no_double_tour(nodes_route, D′) && check_tw_ride_time(copy_route[idx_pre_vi_init:end], eh_pre_vi_init, ei, li, s_t, dist) 
#                         if  no_double_tour(nodes_route, D′) && eight_step(copy_route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH) 
#                             min_dist_increase = dist_increase
#                             best_pos = (pre_vi,vi,v_ni, vi, v_ni,suc_vi) 
#                             best_route_temp = copy(copy_route)
#                             inserted = true   
#                             # @show("test 13, ", user, best_route_temp, inserted)                              
#                         end
#                     end
#                 end
#             end
#             # insert dropoff later
#             idx_2=idx+1
#             eh_pre_v_ni = max(ei[suc_vi], eh_i+s_t[vi]+dist[vi,suc_vi])
#             pre_v_ni = suc_vi
#             suc_v_ni = solution.succ[pre_v_ni]
#             while  idx_2 < n_nodes +2
#                 eh_v_ni = max(ei[v_ni],eh_pre_v_ni+s_t[pre_v_ni]+dist[pre_v_ni,v_ni])
#                 eh_suc_v_ni = max(ei[suc_v_ni],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_v_ni])
#                 eh_pre_v_ni = eh_v_ni
#                 dist_increase = dist[pre_vi, user] + dist[user, suc_vi] + dist[pre_v_ni, v_ni] + dist[v_ni, suc_v_ni]
#                 - ( dist[pre_vi, suc_vi] + dist[pre_v_ni, suc_v_ni])                        
#                 if dist_increase < min_dist_increase 
#                     if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
#                         if eh_v_ni - li[user] - s_t[user] <= Li[user] 
#                             copy_route = copy(route_with_depots)
#                             insert!(copy_route, idx+1, user); insert!(copy_route, idx_2+2, user + n_cus)
#                             nodes_route = nodes[copy_route[2:end-1]]
#                             # if  no_double_tour(nodes_route, D′) && check_tw_ride_time(copy_route[idx_pre_vi_init:end], eh_pre_vi_init, ei, li, s_t, dist) 
#                             if  no_double_tour(nodes_route, D′) && eight_step(copy_route ,route_i, ei, li, s_t, dist, Q, Li, n_cus, TH)  
#                                 min_dist_increase = dist_increase
#                                 best_pos = (pre_vi,vi,suc_vi, pre_v_ni, v_ni,suc_v_ni) 
#                                 best_route_temp = copy(copy_route)
#                                 inserted = true
#                                 # @show("test 14, ", user, best_route_temp, inserted)    
#                                 break                                 
#                             end
#                         end
#                     else
#                         break                        
#                     end
#                 end                
#                 idx_2 += 1 
#                 pre_v_ni = suc_v_ni; 
#                 suc_v_ni = solution.succ[suc_v_ni]                         
#             end
#         end       
#         # end      
#         # inserted == true && break #added on 10-26
#         idx += 1
#         pre_vi = suc_vi
#         suc_vi = solution.succ[suc_vi]
#         suc_v_ni = suc_vi
#     end         
    
#     if inserted == true
#         # @show(route_i, route, best_route_temp, ec_route, parameter_energy.max_ec[route_i], length_route(best_route_temp)*parameter_energy.beta[route_i] )
#         compute_dt_r(best_route_temp, dt_r)
#         veh_id = solution.RI[route_i,5]
#         ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
#         # @show(route_i, route, best_route_temp, ec_route, parameter_energy.max_ec[route_i], length_route(best_route_temp)*parameter_energy.beta[route_i] )
#         if route_i < K+1 && (ec_route > parameter_energy.max_ec[veh_id]) 
#             if insert_charging_route(solution, route_i, best_route_temp, dist_all, ec_route, fast_chg, parameter_energy, nodes, dt_r, dt_r_back)[1]
#                 insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi, walk_tt_bus_stop)
#                 update_auxiliary_var(solution, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#                 # @show("test 15 : ")
#             else
#                 inserted = false
#             end            
#         else
#             insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi, walk_tt_bus_stop)
#             update_auxiliary_var(solution, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#             # @show("test 16 : ")
#         end
#     end
#     return inserted
# end


# # rand insert a user in a feasible position within the extended comptible layers
# function rand_insert_extend(solution::Solution, user, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, start_depot, end_depot, Li, layer_nodes, lyrs_compatible, lyrs_compatible_extend)
    
#     # insert a user to the best positions
#     # init_route_info(e_r,cap_r,l_r,maxcap_r, li, start_depot, end_depot, dt_r, dt_r_back)
#     # init_e_r_cap_r(solution, route_i, start_depot, e_r, cap_r, ei)
#     flag_init_sol = false
#     inserted = false
#     min_dist_increase = bigM * 0.7
#     best_pos = ()    
#     route_with_depots = get_route(solution, route_i, start_depot, end_depot)
#     route = route_with_depots[2:end-1] # no depot
#     vi    = user; v_ni  = user + n_cus
#     lyr_vi= layer_nodes[vi]
#     users = route[findall(x->x<n_cus+2, route)]
#     lyrs_compa_vi = findall(x->x>0, lyrs_compatible_extend[lyr_vi,:])
#     compa_users   = users[findall(x->x∈lyrs_compa_vi, layer_nodes[users])]
#     n_compa_users = length(compa_users)  
#     copy_route = zeros(Int32, length(route_with_depots)+2); nodes_route = Int32[]; best_route_temp = Int32[]
#     if isempty(compa_users)
#         idx= findfirst(x->x>lyr_vi, layer_nodes[route_with_depots]) # the depot layer is 1
#         if isnothing(idx)
#             # inset after the last user of the route
#             best_pos = (route[end],vi,v_ni, vi, v_ni, end_depot) 
#             t0 = e_r[route[end]]
#             idx_insert = length(route_with_depots) 
#             copy_route[1:idx_insert-1] = route_with_depots[1:idx_insert-1]
#             copy_route[idx_insert:end] =[vi, v_ni, end_depot]
#             best_route_temp = copy_route
#             check_tw(best_route_temp[idx_insert-1:end], t0, ei, li, s_t, dist) && (inserted = true) 
#             # @show("test 11, ", user, best_route_temp, inserted)
#             # if !inserted && eight_step(best_route_temp, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH )
#             #     error("test 1: ", " if !inserted && eight_step( best_route_temp, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH )")
#             # end
#         else
#             # insert before route[idx]
#             pre_vi = route_with_depots[idx-1]            
#             suc_v_ni = route_with_depots[idx]
#             best_pos = (pre_vi,vi,v_ni, vi, v_ni, suc_v_ni) 
#             idx_insert = idx   
#             copy_route[1:idx_insert-1] = route_with_depots[1:idx_insert-1]  
#             copy_route[idx_insert:idx_insert+1] =[vi, v_ni]
#             copy_route[idx_insert+2:end] = route_with_depots[idx_insert:end]
#             best_route_temp = copy_route
#             # only need to check tw after the remainig route      
#             # check_tw_ride_time(best_route_temp[idx-1:end], e_r[pre_vi], ei, li, s_t, dist) && (inserted = true) 
#             check_tw_ride_time(best_route_temp[idx-1:end], e_r[pre_vi], ei, li, s_t, dist) && (inserted = true) 
#             # @show("test 12, ", user, best_route_temp, inserted)
#             # if !inserted && eight_step( best_route_temp, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH )
#             #     error("test 2: ", " if !inserted && eight_step( best_route_temp, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH )")
#             # end
#         end
#     else
#         idxs_compa_users = findall(x->x∈compa_users, route_with_depots)
#         idx0= idxs_compa_users[1]
#         # idx_pre_vi_init = idx0-1 
#         # eh_pre_vi_init = e_r[route_with_depots[idx0-1]]
#         # @show(route_with_depots, route_with_depots[idx0-1], compa_users[1], compa_users, layer_nodes[route_with_depots])
#         idx_end  = idx0 + 2*n_compa_users 
#         shuffle!(idxs_compa_users)         
#         idx = idxs_compa_users[1]
#         pre_vi = solution.pre[route_with_depots[idx]] 
#         suc_vi   = route_with_depots[idx]
#         suc_v_ni = route_with_depots[idx]
#         while idx < idx0 + n_compa_users +1
#             eh_i = max(ei[vi], e_r[pre_vi]+s_t[pre_vi]+dist[pre_vi,vi])
#             if eh_i <= li[vi] && cap_r[pre_vi]+qi[vi] <= Q[route_i] && qi[vi] + maxcap_r[suc_vi] <= Q[route_i]
#                 dist_increase = dist[pre_vi, user] + dist[user, v_ni] + dist[v_ni, suc_vi]
#                 - dist[pre_vi,suc_vi] 
#                 if dist_increase < min_dist_increase
#                     eh_v_ni = max(ei[v_ni],eh_i+s_t[vi]+dist[vi,v_ni])
#                     eh_suc_v_ni = max(ei[suc_vi],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_vi])
#                     if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
#                         if eh_v_ni - li[user] - s_t[user] <= Li[user] 
#                             copy_route[1:idx-1], copy_route[idx+2:end] = route_with_depots[1:idx-1], route_with_depots[idx:end]
#                             copy_route[idx:idx+1] = [user, user+n_cus]
#                             nodes_route = nodes[copy_route[2:end-1]]
#                             # if  no_double_tour(nodes_route, D′) && check_tw_ride_time(copy_route[idx_pre_vi_init:end], eh_pre_vi_init, ei, li, s_t, dist) 
#                                 if  no_double_tour(nodes_route, D′) && eight_step(copy_route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH) 
#                                 min_dist_increase = dist_increase
#                                 best_pos = (pre_vi,vi,v_ni, vi, v_ni,suc_vi) 
#                                 best_route_temp = copy(copy_route)
#                                 inserted = true   
#                                 # @show("test 13, ", user, best_route_temp, inserted)                              
#                             end
#                         end
#                     end
#                 end
#                 # insert dropoff later
#                 idx_2=idx+1
#                 eh_pre_v_ni = max(ei[suc_vi], eh_i+s_t[vi]+dist[vi,suc_vi])
#                 pre_v_ni = suc_vi
#                 suc_v_ni = solution.succ[pre_v_ni]
#                 while idx_2 < idx_end && (suc_v_ni > 0)
#                     eh_v_ni = max(ei[v_ni],eh_pre_v_ni+s_t[pre_v_ni]+dist[pre_v_ni,v_ni])
#                     eh_suc_v_ni = max(ei[suc_v_ni],eh_v_ni+s_t[v_ni]+dist[v_ni,suc_v_ni])
#                     eh_pre_v_ni = eh_v_ni
#                     dist_increase = dist[pre_vi, user] + dist[user, suc_vi] + dist[pre_v_ni, v_ni] + dist[v_ni, suc_v_ni]
#                     - ( dist[pre_vi, suc_vi] + dist[pre_v_ni, suc_v_ni])                        
#                     if dist_increase < min_dist_increase 
#                         if  eh_v_ni <= li[v_ni] && eh_suc_v_ni <= li[suc_v_ni] # feasible
#                             if eh_v_ni - li[user] - s_t[user] <= Li[user] 
#                                 copy_route = copy(route_with_depots)
#                                 insert!(copy_route, idx, user); insert!(copy_route, idx_2+1, user + n_cus)
#                                 nodes_route = nodes[copy_route[2:end-1]]
#                                 # if  no_double_tour(nodes_route, D′) && check_tw_ride_time(copy_route[idx_pre_vi_init:end], eh_pre_vi_init, ei, li, s_t, dist) 
#                                     if  no_double_tour(nodes_route, D′) && eight_step(copy_route ,route_i, ei, li, s_t, dist, Q, Li, n_cus, TH)  
#                                     min_dist_increase = dist_increase
#                                     best_pos = (pre_vi,vi,suc_vi, pre_v_ni, v_ni,suc_v_ni) 
#                                     best_route_temp = copy(copy_route)
#                                     inserted = true
#                                     # @show("test 14, ", user, best_route_temp, inserted)    
#                                     break                                 
#                                 end
#                             end
#                         else
#                             break                        
#                         end
#                     end                
#                     idx_2 += 1 
#                     pre_v_ni = suc_v_ni; 
#                     suc_v_ni = solution.succ[suc_v_ni]                         
#                 end
#             end       
#             # end      
#             inserted == true && break #added on 10-26
#             idx += 1
#             pre_vi = suc_vi
#             suc_vi = solution.succ[suc_vi]
#             suc_v_ni = suc_vi
#         end         
#     end
#     if inserted == true
#         # @show(route_i, route, best_route_temp, ec_route, parameter_energy.max_ec[route_i], length_route(best_route_temp)*parameter_energy.beta[route_i] )
#         compute_dt_r(best_route_temp, dt_r)
#         veh_id = solution.RI[route_i,5]
#         ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
#        if route_i < K+1 && (ec_route > parameter_energy.max_ec[veh_id]) 
#             if insert_charging_route(solution, route_i, best_route_temp, dist_all, ec_route, fast_chg, parameter_energy, dt_r, dt_r_back, flag_init_sol)[1]
#                 insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
#                 update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#             else
#                 inserted = false
#             end            
#         else
#             insert_user_route(solution, route_i, best_pos, start_depot, end_depot, qi)
#             update_auxiliary_var(solution, darp, route_i, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
#         end
#     end
#     return inserted
# end


# # to insert charging for a route with only one user (need to recharge after the start_depot)
# function insert_charging_simple(solution::Solution, r, route::Vector{Int32}, dist_all, ec_route, fast_chg, parameter_energy, set_physical, nodes, dt_r, dt_r_back )
    
#     veh_id = solution.RI[r,5]
#     v= (route[1], route[2])
#     idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg)
#     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
#     vec_delta_e  =  parameter_energy.E_min[veh_id] - (parameter_energy.E_init[veh_id] - (ec_route + parameter_energy.beta[veh_id] * (t_access_1+ t_access_2 - t_v_v1)))
#     t_chg = vec_delta_e / fast_chg.pw[idx_chg]
#     forward_time_slack, vec_t_start = get_forward_time_slack(route, ei, li, s_t, dist, n_cus, qi, [1], dist_orig)
#     if t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[1]
#         additional_time = t_access_1 + t_access_2 - t_v_v1 +  t_chg
#         t_start_chg = vec_t_start[1] + t_access_1 + rand(1)[1]* (forward_time_slack[1]-t_chg) # random start time
#         solution.vec_chg_event[r,1:6] = [additional_time, 1, v[1], idx_chg, t_start_chg, t_start_chg+t_chg]
#         # check_energy(solution, route, r, parameter_energy, fast_chg)
#         return true
#     else 
#         return false
#     end
# end

# function reject_cus(solution::Solution, r, route, veh_id, parameter_energy, e_r, cap_r,l_r,maxcap_r, dt_r, dt_r_back)

#     max_ec = parameter_energy.max_ec[veh_id]
#     dt_r_eng = dt_r[route] .* parameter_energy.beta[veh_id]
#     i_ec_violate = findfirst(x->x-max_ec>0, dt_r_eng)
#     route_tmp = route[1:i_ec_violate-1]; route_tmp_new=[start_depot]
#     for v in route_tmp
#         if v < n_cus+2 && (v+n_cus ∈ route_tmp)
#             union!(route_tmp_new, v, v+n_cus)
#         end
#     end
    
#     n_nodes_new= length(route_tmp_new)
#     if n_nodes_new > 2
#         route_seg1, route_seg2 = route_tmp_new, route[n_nodes_new+1:end]
#         solution.succ[route_seg1[end]] = end_depot 
#         rejected_users = route_seg2[findall(x->x<n_cus+2, route_seg2)]
#         for user in rejected_users          
#             solution.penalty_unserved += qi[user]*penalty
#             solution.TT_walk_routes[r] -= walk_tt_bus_stop[user-1]
#             union!(solution.unserved_users, user)
#             solution.pre[user]=solution.succ[user]=-1
#             solution.pre[user+n_cus]=solution.succ[user+n_cus]=-1
#         end
#         #update solution
#         solution.RI[r,2:3] = [route_seg1[end], n_nodes_new-1] 
#         reset_auxiliary_route(solution, r, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
#     else
#         # remove the route
#         rejected_users = route[findall(x->1<x<n_cus+2, route)]
#         for user in rejected_users          
#             solution.penalty_unserved += qi[user]*penalty
#             solution.TT_walk_routes[r] -= walk_tt_bus_stop[user-1]
#             union!(solution.unserved_users, user)
#             solution.pre[user]=solution.succ[user]=-1
#             solution.pre[user+n_cus]=solution.succ[user+n_cus]=-1
#         end
#         #update solution
#         n_route = solution.n_route
#         solution.n_route -= 1
#         for i in r+1:n_route
#             solution.RI[i-1,:] = solution.RI[i,:]
#             solution.TT_walk_routes[i-1]   = solution.TT_walk_routes[i]
#             solution.vec_chg_event[i-1,:]  = solution.vec_chg_event[i,:]
#         end  
#     end
# end


# function split_route(solution::Solution, instance, darp, fast_chg, parameter_energy, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back) # user is pickup node is from 2,...,n+1
    
#     # solution = sol_current_best
#     start_depot, end_depot = darp.start_depot, darp.end_depot
#     dist_all, dist_orig, set_physical = darp.dist_all, darp.dist_orig, darp.set_physical
#     T_max, K, K_MAX = instance.T_max, instance.K, instance.K_MAX

#     n_route = solution.n_route
#     id_long_route = 0
#     conflicts =zeros(Int32, n_route, n_route)
#     for i in 1:n_route-1        
#         for j in i+1:n_route 
#             conflicts[i,j]= chg_conflict_two_routes(solution, i, j, occ_state_chg_init, discrte_t_interval)
#             conflicts[j,i] = conflicts[i,j]                 
#         end
#     end

#     vec_tmp =vec(sum(conflicts, dims=2))
#     val_max = maximum(vec_tmp)
#     candidate= findall(x->x==val_max, vec_tmp)
#     if length(candidate) >1 
#         vec_route_length = solution.RI[candidate,3]
#         id_long_route= candidate[argmax(vec_route_length)]
#     else
#         id_long_route = candidate[1]
#     end
#     route_i = id_long_route

#     route = get_route(solution, route_i, start_depot, end_depot)
#     # ! flag_silent && @show(conflicts, id_long_route, route)
#     users = route[findall(x->1<x<n_cus+2, route)]
#     n_user_remove = floor(Int32, length(users)/2)
#     shuffle!(users)
#     unserved= users[1:n_user_remove] 
#     user = unserved[1]
#     # get veh_id for this new ev
#     veh_ids_sol= solution.RI[1:n_route,5]
#     veh_ids_candidate = collect(1:K_MAX) 
#     setdiff!(veh_ids_candidate, veh_ids_sol) 
#     idx_tmp = findfirst(x->x==true, parameter_energy.is_electric[veh_ids_candidate])
#     if isnothing(idx_tmp) 
#         error("not sufficient size of ev fleet, please increase K_MAX !!") 
#     end
#     veh_id_new = veh_ids_candidate[idx_tmp]  # list gv
#     solution.n_route +=1  
#     r = solution.n_route  
#     solution.RI[r,5] = veh_id_new # find unused ev 
    
#     insert_consecutive(solution, r, start_depot, user, end_depot, n_cus, qi) 
#     update_auxiliary_var(solution, darp, r, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem 
    
#     filter!(x->x∉[user, user + n_cus], route) 
#     route_auxiliary_update(solution, darp, route , route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back) 
#     setdiff!(unserved, user)
#     flag_init_sol = true # using rand charger insert policy
#     routes = setdiff!(collect(1:r), [route_i, r])
#     shuffle!(routes)    
#     for user in unserved
#         done = false
#        if rand() > 0.5  
#             if greedy_insert(solution, fast_chg, global_parameter, instance, darp, user, r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, flag_init_sol) 
#                 filter!(x->x∉[user, user + n_cus], route); done = true
#                 break
#             end
#         else
#             if rand_insert(solution, fast_chg, global_parameter, instance, darp, user, r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, start_depot, end_depot, Li) 
#                 filter!(x->x∉[user, user + n_cus], route); done = true
#                 break
#             end
#         end  
#         if ! done  
#             error("cannot inser user in split_route() error !! please create a new route !!")
#         end        
#     end
#     if  r > K_MAX
#         error("solution.n_route > instance.K_MAX error ! in split_route(), please increase the available EVs !!")
#     end 
#     if solution.vec_chg_event[route_i, 1] > 0 
#         _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(solution, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
#         if  !chg_succ 
#             repair_charging(solution, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start)
#         end
#     end
     
#     route_r= get_route(solution, r, start_depot, end_depot)
#     update_charging(solution, darp, r, route_r, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
#     reset_auxiliary_route(solution, darp, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
#     reset_auxiliary_route(solution, darp, r,       e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
# end


# function split_route(solution::Solution, instance, darp, fast_chg, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, occ_state_chg_init, discrte_t_interval) # user is pickup node is from 2,...,n+1
    
#     # solution = sol_current_best
#     start_depot, end_depot = darp.start_depot, darp.end_depot
#     dist_all, set_physical = darp.dist_all, darp.set_physical
#     T_max, K, K_MAX = instance.T_max, instance.K, instance.K_MAX
#     parameter_energy = instance.parameter_energy

#     n_route = solution.n_route
#     id_long_route = 0
#     conflicts =zeros(Int32, n_route, n_route)
#     for i in 1:n_route-1        
#         for j in i+1:n_route 
#             conflicts[i,j]= chg_conflict_two_routes(solution, i, j, occ_state_chg_init, discrte_t_interval)
#             conflicts[j,i] = conflicts[i,j]                 
#         end
#     end
#     vec_tmp =vec(sum(conflicts, dims=2))
#     # solution.vec_chg_event[1:n_route,:]
#     val_max = maximum(vec_tmp)
#     candidate= findall(x->x==val_max, vec_tmp)
#     if length(candidate) >1 
#         vec_route_length = solution.RI[candidate,3]
#         id_long_route= candidate[argmax(vec_route_length)]
#     else
#         id_long_route = candidate[1]
#     end
#     r_i = id_long_route

#     route = get_route(solution, r_i, start_depot, end_depot) 
#     reset_cap_r(solution, r_i, cap_r, darp)
#     cap_r_ri = [cap_r[vi] for vi in route] 
#     idx_ri   = findall(x->x==0, cap_r_ri) # find location of zero load  
#     if isnothing(idx_ri) 
#         error("idx_ri empty error in split_route() !! please check !!") 
#     end  
#     idx_middle = solution.RI[r_i,3] / 2+1 
#     vec_diff =  idx_ri .- idx_middle
#     abs_vec_diff = abs.(vec_diff)
#     id_end_route_ri_new = idx_ri[argmin(abs_vec_diff)]
#     length_new_route = solution.RI[r_i,3] - (id_end_route_ri_new-1)
#     # update the current route
#     vi = route[id_end_route_ri_new]
#     v_route_end_1 = route[end-1]
#     suc_vi = solution.succ[vi]
#     solution.succ[vi] = end_depot
#     solution.RI[r_i,2]=vi; solution.RI[r_i,3] = id_end_route_ri_new-1
#     # create another route
#     # get veh_id for this new ev
#     veh_ids_sol= solution.RI[1:n_route,5]
#     veh_ids_candidate = collect(1:K_MAX) 
#     setdiff!(veh_ids_candidate, veh_ids_sol) 
#     idx_tmp = findfirst(x->x==true, parameter_energy.is_electric[veh_ids_candidate])
#     if isnothing(idx_tmp) 
#         error("not sufficient size of ev fleet, please increase K_MAX !!") 
#     end
#     veh_id_new = veh_ids_candidate[idx_tmp]  # list gv
#     solution.n_route +=1  
#     r = solution.n_route  
#     solution.RI[r,5] = veh_id_new # find unused ev 
#     solution.pre[suc_vi] = start_depot
#     solution.RI[r,1] = suc_vi ; solution.RI[r,2] = v_route_end_1;
#     solution.RI[r,3] = length_new_route;
#     if  r > K_MAX
#         error("solution.n_route > instance.K_MAX error ! in split_route(), please increase the available EVs !!")
#     end 

#     route = get_route(solution, r_i, start_depot, end_depot)
#     _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(solution, darp, r_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
#     ! chg_succ &&  repair_charging(solution, darp, r_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start)
#     route_r= get_route(solution, r, start_depot, end_depot)
#     _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(solution, darp, r, route_r, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
#     ! chg_succ &&  repair_charging(solution, darp, r, route_r, dist_all, fast_chg, parameter_energy, set_physical, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start)

# end