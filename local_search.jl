# local search operators

# 2-opt :https://en.wikipedia.org/wiki/2-opt
# i->i+1->...j->j+1 => i->j-...->i+1->j+1
function two_opt(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
 
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    TH, ei, li, s_t, qi, dist=  darp.TH, darp.ei, darp.li, darp.s_t, darp.qi, darp.dist
 
    n_route = solution.n_route
    if n_route < 1  
       return false
    end
    route_i = rand(1:n_route)  
    n_nodes = solution.RI[route_i,3]
    route = get_route(solution, route_i, start_depot, end_depot)
    new_route = zeros(Int32,2*n_cus+2) ; r_new=zeros(Int32,2*n_cus+2)
    nodes_route = Int32[]
    dist_min = 0 ; done = false 
    vj=0;vj_1=0   
    if n_nodes > 3
        max_dist_segs = min(n_nodes-2, 4) # max number of users to be reversed between i and j, including j, 4 is a parameter
        for i in 2:n_nodes-max_dist_segs
            vi = route[i];vi_1= route[i+1]
            for j in i+2: i+max_dist_segs 
                vj = route[j]; vj_1 = route[j+1]
                seg = reverse!(route[i+1:j]) 
                n_seg = length(seg)
                dist_saving = dist[vi,vi_1] + dist[vj,vj_1] - (dist[vi, vj] + dist[vi_1,vj_1])
                if dist_saving > dist_min                 
                    if check_pre_seg(seg, n_cus)   # ok
                        # check tw feasiblity of changed segments
                        feasible = false
                        eh_vs = e_r[vi]
                        eh_vs_1 = max(ei[vj], eh_vs+s_t[vi]+dist[vi,vj])
                        cap_r_temp  =  cap_r[vi] + qi[vj]
                        if eh_vs_1 <= li[vj] || cap_r_temp  <= Q[route_i]
                            feasible = true
                            eh_vs = eh_vs_1
                            # check capacity 
                            for k in 1:length(seg)-1
                                vs_1=seg[k+1];vs=seg[k]
                                cap_r_temp += qi[seg[k+1]]
                                eh_vs_1 = max(ei[vs_1], eh_vs+s_t[vs]+dist[vs,vs_1])
                                if  eh_vs_1 > li[vs_1] ||  cap_r_temp > Q[route_i]
                                    feasible = false
                                    break
                                end
                                eh_vs = eh_vs_1
                            end
                            #check vj_1
                            if feasible
                                eh_vs_1 = max(ei[vj_1], eh_vs+s_t[vi_1]+dist[vi_1,vj_1])
                                if  eh_vs_1 > li[vj_1] ||  cap_r_temp + qi[vj_1]  > Q[route_i]
                                    feasible = false
                                end
                            end
                        end
                        if  feasible
                            new_route[1:i], new_route[i+1:i+n_seg] = route[1:i], seg
                            new_route[i+n_seg+1:n_nodes+2]= route[j+1:end]
                            nodes_route = new_route[2:n_seg+1]
                            # if no_double_tour(nodes_route, D′) && eight_step(new_route[1:n_nodes+2],route_i, ei, li, s_t, dist, Q, Li, n_cus, TH)
                            if eight_step(new_route[1:n_nodes+2],route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)
                                r_new= copy(new_route[1:n_nodes+2])
                                dist_min = dist_saving
                                done = true
                            end
                        end
                    end
                end              
            end
        end
        # exchange swap them
        done && route_auxiliary_update(solution, darp, r_new, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back) 
    end
    return false
end


# compute relatedness for pair of requests
function compute_relatedness(darp, global_parameter, relatedness)

    phi = 9; chi = 3; psi = 2
    dist, li, qi, n_cus = darp.dist, darp.li, darp.qi, darp.n_cus
    bigM = global_parameter.bigM
    max_qi= maximum(qi)

    dist_normalized = dist./bigM
    li_normalized = li./li[end]
    qi_normalized = qi./max_qi

    for i in collect(2:n_cus)
        for j in collect(i+1:n_cus+1)
            relatedness.dist[i,j] = dist_normalized[i,j] + dist_normalized[i+n_cus, j+n_cus]
            relatedness.tw[i,j]   = abs(li_normalized[i] - li_normalized[j]) + abs(li_normalized[i+n_cus] - li_normalized[j+n_cus])
            tmp = abs(qi_normalized[i] - qi_normalized[j])
            relatedness.shaw[i,j] = phi * relatedness.dist[i,j] + chi * relatedness.tw[i,j] + psi * tmp
            relatedness.dist[j,i] = relatedness.dist[i,j]
            relatedness.tw[j,i]   = relatedness.tw[i,j]
            relatedness.shaw[j,i] = relatedness.shaw[i,j]
        end
    end 
end  




# update charging after remove requests
function update_charging_af_remove(soluT::Solution, fast_chg, instance, darp, dt_r, dt_r_back)

    start_depot, end_depot = darp.start_depot, darp.end_depot
    dist_all, set_physical = darp.dist_all, darp.set_physical 
    K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy
    n_charger_installed= 0; success=false; pos_insert_chg =0; forward_time_slack=0; vec_t_start=0
    for r_i in collect(1:soluT.n_route)
        route = get_route(soluT, r_i, start_depot, end_depot)
        n_charger_installed, success, pos_insert_chg, forward_time_slack, vec_t_start =update_charging(soluT, darp, r_i, route, dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
        if ! success && (n_charger_installed > 0)                            
            repair_charging(soluT, darp, r_i, route, fast_chg, parameter_energy, dt_r, pos_insert_chg, forward_time_slack, vec_t_start)
        end
        # check_energy(soluT, darp, route, r_i, parameter_energy, fast_chg)
    end       
end

# rand remove n_remove requests
function rand_remove(soluT::Solution, darp, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness::Relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)   
  
    # @show("rand_remove call!")
    r_i = 0 ; p = 6 # mean : 0.12
    n_route= soluT.n_route
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    route_of_requests = zeros(Int32, 2*n_cus+2) 
    vec_routes = Vector{Vector{Int32}}(undef, n_route)
    length(unserved_users) == (n_cus-1) && return
    for r in collect(1:n_route)
        vec_routes[r] = get_route(soluT, r, start_depot, end_depot)
        for node in vec_routes[r]
            route_of_requests[node] = r
        end  
    end
    users= shuffle(candidate_users_remove)[1:n_remove]
    for user in users 
        length(unserved_users) == (n_cus-1) && break
        r_i = route_of_requests[user]
        union!(unserved_users, user) 
        filter!(x->x∉[user, user + n_cus], vec_routes[r_i])
        update_tmp_sol(soluT, user, route_of_requests, r_i, vec_routes, darp, e_r,cap_r,l_r, maxcap_r, dt_r, dt_r_back)
        # if n_route==1 && (length(vec_routes[r_i]) == 4)
        #     break
        # end
    end
    update_charging_af_remove(soluT, fast_chg, instance, darp, dt_r, dt_r_back)
    soluT.unserved_users =  unserved_users
   
end

# remove n_remove worst requests according their cost 
# p_worst = 3, based on Ropke and Pisinger (2016)  mean 0.21,
function worst_remove(soluT::Solution, darp, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness::Relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
    
    p_worst = 3   
    n_route= soluT.n_route
    n_cus, start_depot, end_depot, dist = darp.n_cus, darp.start_depot, darp.end_depot, darp.dist   
    bigM=1000
    length(unserved_users) == (n_cus-1) && return
    route_of_requests = zeros(Int32, 2*n_cus+2)  
    vec_routes = Vector{Vector{Int32}}(undef, n_route)
    vec_cost= zeros(n_cus+1) 

    for r in collect(1:n_route)
        vec_routes[r] = get_route(soluT, r, start_depot, end_depot)
        for node in  vec_routes[r]
            route_of_requests[node] = r
        end  
    end
    users_routes= Int32[]
    for r in collect(1:n_route) 
        route = get_route_no_depot(soluT, r)        
        users = route[findall(x->x<n_cus+2, route)]
        union!(users_routes, users)
        for v in users
            v_ni = v+n_cus
            pre_v = soluT.pre[v]; suc_v = soluT.succ[v]
            pre_v_ni =soluT.pre[v_ni]; suc_v_ni = soluT.succ[v_ni]
            vec_cost[v] = dist[pre_v,v]+dist[v,suc_v]+dist[pre_v_ni,v_ni]+dist[v_ni,suc_v_ni]
                        - (dist[pre_v,suc_v] + dist[pre_v_ni,suc_v_ni] )
        end
    end
    vec_cost_users  = vec_cost[users_routes]
    idx_sorted_cost = sortperm(vec_cost_users, rev=true)  
    candidate_users_remove = candidate_users_remove[idx_sorted_cost] # in descending order
    n_candidate = length(candidate_users_remove)

    for idx in collect(1:n_remove)   
        length(unserved_users) == (n_cus-1) && break
        y_p = rand()^p_worst
        idx = max(1, floor(Int32, y_p * n_candidate))   
        user_remove = candidate_users_remove[idx]
        setdiff!(candidate_users_remove, user_remove)
        union!(unserved_users, user_remove)
        r_remove = route_of_requests[user_remove] 
        filter!(x->x∉[user_remove, user_remove + n_cus], vec_routes[r_remove])
        update_tmp_sol(soluT, user_remove, route_of_requests, r_remove, vec_routes, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)   
        n_candidate -= 1
        n_candidate == 0 && break  
        # if n_route==1 && (length(vec_routes[r_remove]) == 4)
        #     break
        # end
    end 
    update_charging_af_remove(soluT, fast_chg, instance, darp, dt_r, dt_r_back)
    soluT.unserved_users = unserved_users
    
end 

# remove n_remove requests based distance relatedness
# p = 6 based on Ropke and Pisinger (2016) (mean : 0.12)
function dist_remove(soluT::Solution, darp, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness::Relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
    
    # @show("dist_remove call!")
    r_i = 0 ; p = 6 # mean : 0.12
    n_route= soluT.n_route; user_remove = 0
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    route_of_requests = zeros(Int32, 2*n_cus+2) 
    vec_routes = Vector{Vector{Int32}}(undef, n_route)
    length(unserved_users) == (n_cus-1) && return
    for r in collect(1:n_route)
        vec_routes[r] = get_route(soluT, r, start_depot, end_depot)
        for node in  vec_routes[r]
            route_of_requests[node] = r
        end  
    end
    user = rand(candidate_users_remove)
    r_i = route_of_requests[user]
    union!(unserved_users, user)
    setdiff!(candidate_users_remove, user)
    filter!(x->x∉[user, user + n_cus], vec_routes[r_i])
    update_tmp_sol(soluT, user, route_of_requests, r_i, vec_routes, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
    n_candidate= length(candidate_users_remove)
    for _ in collect(1:n_remove-1)
        length(unserved_users) == (n_cus-1) && break
        user = rand(unserved_users)
        r_i = route_of_requests[user] 
        if length(candidate_users_remove) >1        
            rel = relatedness.dist[user, candidate_users_remove] 
            idx_relatedness_sorted = sortperm(rel)
            y_p = rand()^p
            idx = max(1, floor(Int32, y_p * n_candidate))        
            # @show(idx, idx_relatedness_sorted, idx_relatedness_sorted[idx], candidate_users_remove, n_candidate)
            user_remove = candidate_users_remove[idx_relatedness_sorted[idx]]
        else
            user_remove = candidate_users_remove[1]
        end
        union!(unserved_users, user_remove)
        r_remove = route_of_requests[user_remove]
        setdiff!(candidate_users_remove, user_remove)
        filter!(x->x∉[user_remove, user_remove + n_cus], vec_routes[r_remove])
        update_tmp_sol(soluT, user_remove, route_of_requests, r_remove, vec_routes, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
        n_candidate -= 1 
        n_candidate == 0 && break  
        # if n_route==1 && (length(vec_routes[r_remove]) == 4)
        #     break
        # end
    end  
    update_charging_af_remove(soluT, fast_chg, instance, darp, dt_r, dt_r_back)
    soluT.unserved_users = unserved_users

end

# used li instead of arrival time 
# remove n_remove requests based time window relatedness
function tw_remove(soluT::Solution, darp, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness::Relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
    
    # @show("tw_remove call!")
    r_i = 0 ; p = 6 #based on Ropke and Pisinger (2016)
    n_route= soluT.n_route; user_remove = 0
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    route_of_requests = zeros(Int32, 2*n_cus+2) 
    vec_routes = Vector{Vector{Int32}}(undef, n_route)
    length(unserved_users) == (n_cus-1) && return

    for r in collect(1:n_route)
        vec_routes[r] = get_route(soluT, r, start_depot, end_depot)
        for node in  vec_routes[r]
            route_of_requests[node] = r
        end  
    end
    user = rand(candidate_users_remove)
    r_i = route_of_requests[user]
    union!(unserved_users, user)
    setdiff!(candidate_users_remove, user)
    filter!(x->x∉[user, user + n_cus], vec_routes[r_i])
    update_tmp_sol(soluT, user, route_of_requests, r_i, vec_routes, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
    n_candidate= length(candidate_users_remove)
    for _ in collect(1:n_remove-1)
        length(unserved_users) == (n_cus-1) && break
        user = rand(unserved_users)
        r_i = route_of_requests[user]
        if length(candidate_users_remove) >1         
            rel = relatedness.tw[user, candidate_users_remove] 
            idx_relatedness_sorted = sortperm(rel)
            y_p = rand()^p
            idx = max(1, floor(Int32, y_p * n_candidate))
            # @show(idx, idx_relatedness_sorted, idx_relatedness_sorted[idx], candidate_users_remove, n_candidate)
            user_remove = candidate_users_remove[idx_relatedness_sorted[idx]]
        else
            user_remove = candidate_users_remove[1]
        end
        union!(unserved_users, user_remove)
        r_remove = route_of_requests[user_remove]
        setdiff!(candidate_users_remove, user_remove)
        filter!(x->x∉[user_remove, user_remove + n_cus], vec_routes[r_remove])
        update_tmp_sol(soluT, user_remove, route_of_requests, r_remove, vec_routes, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
        n_candidate -= 1  
        n_candidate == 0 && break   
        # if n_route==1 && (length(vec_routes[r_remove]) == 4)
        #     break
        # end
    end 
    update_charging_af_remove(soluT, fast_chg, instance, darp, dt_r, dt_r_back)
    soluT.unserved_users = unserved_users
end

# difference with Ropke and Pisinger (2016), I use li instead of arrival time, and do not consider the last term (|Ki|-|Kj|/min(|Ki|,|Kj|))(computational expensive)
# remove n_remove based on adapted Shaw relatedness
function shaw_remove(soluT::Solution, darp, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness::Relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
      
    # @show("shaw_remove call!")
    
    r_i = 0 ; p = 6 # based on Ropke and Pisinger (2016) mean : 0.12
    n_route= soluT.n_route; user_remove = 0
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    route_of_requests = zeros(Int32, 2*n_cus+2) 
    vec_routes = Vector{Vector{Int32}}(undef, n_route)
    length(unserved_users) == (n_cus-1) && return

    for r in collect(1:n_route)
        vec_routes[r] = get_route(soluT, r, start_depot, end_depot)
        for node in  vec_routes[r]
            route_of_requests[node] = r
        end  
    end
    user = rand(candidate_users_remove)
    r_i = route_of_requests[user]
    union!(unserved_users, user)
    
    setdiff!(candidate_users_remove, user)
    filter!(x->x∉[user, user + n_cus], vec_routes[r_i])
    update_tmp_sol(soluT, user, route_of_requests, r_i, vec_routes, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
    n_candidate= length(candidate_users_remove)
    for _ in collect(1:n_remove-1)
        length(unserved_users) == (n_cus-1) && break
        user = rand(unserved_users)
        r_i = route_of_requests[user]     
        if length(candidate_users_remove) >1
            rel = relatedness.shaw[user, candidate_users_remove] 
            idx_relatedness_sorted = sortperm(rel)
            y_p = rand()^p
            idx = max(1, floor(Int32, y_p * n_candidate))        
            # @show(idx, idx_relatedness_sorted, idx_relatedness_sorted[idx], candidate_users_remove, n_candidate)
            user_remove = candidate_users_remove[idx_relatedness_sorted[idx]]
        else
            user_remove = candidate_users_remove[1]
        end
        union!(unserved_users, user_remove)
        r_remove = route_of_requests[user_remove]
        setdiff!(candidate_users_remove, user_remove) 
        filter!(x->x∉[user_remove, user_remove + n_cus], vec_routes[r_remove]) 
        update_tmp_sol(soluT, user_remove, route_of_requests, r_remove, vec_routes, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
        n_candidate -= 1   
        n_candidate == 0 && break   
        # if n_route==1 && (length(vec_routes[r_remove]) == 4)
        #     break
        # end
    end 
    update_charging_af_remove(soluT, fast_chg, instance, darp, dt_r, dt_r_back)
    soluT.unserved_users = unserved_users

end

# wrapper for remove operators
function RMOP(selected_rm_op, soluT::Solution, darp, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness::Relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back) 
   
    selected_rm_op(soluT::Solution, darp, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness::Relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
end

# wrapper for insert operators
function INOP(selected_in_op, soluT::Solution, fast_chg, global_parameter, recorder_routes, recorder_lscost_routes, instance, darp, unserved_users, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back) 
   
    selected_in_op(soluT::Solution, fast_chg, global_parameter, recorder_routes, recorder_lscost_routes, instance, darp, unserved_users, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back)
end


# randomly run one of the three relocation LS operators
 function relocate_ensemble(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)

    ls_relocation= [relocate, relocate_worst] 
    ls_idx= rand(collect(1:length(ls_relocation)))
    ls = ls_relocation[ls_idx]
    return LS(ls, solution, soluT, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, T, Li)

end

# insert unserved requests to the routes of the current solution
function insert_unserved(soluT::Solution, recorder_routes, recorder_lscost_routes, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
          
    unserved_users = soluT.unserved_users     
    if length(unserved_users) >0         
        rand() > 0.5 ? idx_selected = 1 : idx_selected = 2
        vec_insert_op = [greedy_insert_unserved, regret_insert_unserved]
        selected_in_op = vec_insert_op[idx_selected]
        INOP(selected_in_op, soluT, fast_chg, global_parameter, recorder_routes, recorder_lscost_routes, instance, darp, unserved_users, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back) 
    end
end

# relocate a worst user to the least cost position of the other routes
function relocate_worst(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    qi = darp.qi
    dist_all, set_physical, dist  = darp.dist_all, darp.set_physical, darp.dist
    K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy  
    penalty = global_parameter.penalty
    flag_init_sol = false; success = false
    update_sol(soluT, solution) #copy solution as temporary solution soluT
    n_route = soluT.n_route; route = Int32[]; routes = Int32[]
    
    route_i = rand(1:n_route)  
    route = get_route(soluT, route_i, start_depot, end_depot)
    user = worst_user(soluT, qi, route_i, start_depot, end_depot, n_cus, dist) 
    success = greedy_insert_request(soluT, fast_chg, global_parameter, recorder_routes, instance, darp, user, route_i, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back)
 
    if success  
        if length(route)==4  # remove this route
            soluT.RI[end, 4:5] =  soluT.RI[route_i,4:5] # move this unused veh (capacity and veh_id) to the end of the list
            soluT.n_route -= 1
            for i in route_i+1:n_route
                soluT.RI[i-1,:] = soluT.RI[i,:]
                soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
            end     
            update_sol(solution, soluT)            
            return true  
        else #update this route         
            filter!(x->x∉[user, user + n_cus], route) 
            route_auxiliary_update(soluT, darp, route , route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)  
            # if eight_step(route, route_i, ei, li, s_t, dist_orig, Q, Li, n_cus, TH)
            if soluT.vec_chg_event[route_i, 1] > 0 
                _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
                if  chg_succ
                    update_sol(solution, soluT)  
                else
                    repair_charging(soluT, darp, route_i, route, fast_chg, parameter_energy, dt_r, pos_insert_chg, forward_time_slack, vec_t_start) && update_sol(solution, soluT)
                end
            end 
            return true 
        end 
    end  
    return false
end 


# regret insert all unserved customers to the other routes based on Ropke and Pisinger (2016)
 function regret_insert_unserved(soluT::Solution, fast_chg, global_parameter, recorder_routes, recorder_lscost_routes, instance, darp, unserved_users, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back)
    
    # @show("regret_insert_unserved call!")
    n_route = soluT.n_route # num of candidate routes
    if n_route < 2
       greedy_insert_unserved(soluT, fast_chg, global_parameter, recorder_routes, recorder_lscost_routes, instance, darp, unserved_users, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back)
       return
    end 
    n_cus, start_depot, end_depot, qi = darp.n_cus, darp.start_depot, darp.end_depot, darp.qi 
    dist_all, dist_orig, set_physical = darp.dist_all, darp.dist_orig, darp.set_physical
    T_max, K = instance.T_max, instance.K
    parameter_energy     =  instance.parameter_energy   
    penalty              = global_parameter.penalty
    
    best_pos, best_route_temp= Int32[], Int32[]
    
    flag_init_sol = false; 
    routes= collect(1:n_route)
    
    route_r = 0; idx_min_cost = 0
    vec_cost_routes = zeros(n_route) # cost_increase, route_i, best_position
    for user in unserved_users
        for (idx, _r) in enumerate(routes)
            vec_cost_routes[idx], best_pos, best_route_temp = cost_greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, flag_init_sol) 
            recorder_routes[idx].r, recorder_routes[idx].best_pos, recorder_routes[idx].new_route = _r, best_pos, best_route_temp
        end
        vec_idx_sorted_cost= sortperm(vec_cost_routes) # sorted index set
        idx_min_cost = vec_idx_sorted_cost[1] # index of least cost row     
        min_cost = vec_cost_routes[idx_min_cost]
        # k_regret 
        sum_regret = 0 
        n_route_regret = min(3, n_route) # we compute 2-regret and 3-regret only
        for idx in collect(2:n_route_regret) # need to have at least 2 routes for computing regret, see the paper of Ropke and Pisinger (2006)
            sum_regret += (vec_cost_routes[vec_idx_sorted_cost[idx]] - min_cost)
            recorder_lscost_routes[user].regret[idx-1]=  sum_regret # first column 2-regret, second column: 3-regret
        end 
        recorder_lscost_routes[user].insert_cost = vec_cost_routes[idx_min_cost] 
        recorder_lscost_routes[user].r           = recorder_routes[idx_min_cost].r 
        recorder_lscost_routes[user].best_pos    = recorder_routes[idx_min_cost].best_pos 
        recorder_lscost_routes[user].new_route   = recorder_routes[idx_min_cost].new_route
    end
    if n_route > 3
        k_regret = rand(collect(1:2)) # 2-regret or 3-regret
    else
        k_regret = 1  #2-regret
    end
    vec_regret= [recorder_lscost_routes[user].regret[k_regret] for user in unserved_users]
    idx_vec_regret = sortperm(vec_regret, rev= true) #
    sorted_unserved_users = unserved_users[idx_vec_regret]
    # @show(unserved_users, sorted_unserved_users, length(unserved_users), length(sorted_unserved_users))
    flag_route_changed = zeros(Bool, n_route)
    copy_unserved = Int32[]
    for user in sorted_unserved_users
        route_r = recorder_lscost_routes[user].r  # route to insert the rejected request
        veh_id = soluT.RI[route_r,5]
        best_pos, best_route_temp = recorder_lscost_routes[user].best_pos, recorder_lscost_routes[user].new_route
        # @show(user,route_r, veh_id, best_pos, best_route_temp, flag_route_changed[route_r] )
        if length(best_pos) == 0
            push!(copy_unserved, user)
        else
            if ! flag_route_changed[route_r] 
                if parameter_energy.is_electric[veh_id] == true # EV 
                    compute_dt_r(best_route_temp, dt_r, dist_orig)
                    ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
                    if ec_route > parameter_energy.max_ec[veh_id]
                        if fast_chg.n_fast_chg_installed > 0
                            if insert_charging_route(soluT, darp, route_r, best_route_temp, dist_all, ec_route, fast_chg, parameter_energy, dt_r, dt_r_back, flag_init_sol)[1]
                                insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                                update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem 
                                flag_route_changed[route_r] = true
                            else
                                push!(copy_unserved, user)
                            end 
                        else
                            push!(copy_unserved, user)
                        end
                    else    
                        soluT.vec_chg_event[route_r,:] .*= 0           
                        insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                        update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem  
                        flag_route_changed[route_r] = true
                    end
                else
                    soluT.vec_chg_event[route_r,:] .*= 0 
                    insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                    update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
                    flag_route_changed[route_r] = true
                end  
            else  # route has been changed due to inserting precedent unserved customers, so we insert cus to the best pos of this route
                if !greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, route_r, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back, flag_init_sol)
                    push!(copy_unserved, user)
                end
            end
        end
    end

    # try again 
    remaining_unserved = Int32[]
    for user in copy_unserved
        success = greedy_insert_request(soluT, fast_chg, global_parameter, recorder_routes, instance, darp, user, 0, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back)
        if ! success
           push!(remaining_unserved, user)
        end
    end
    # @show(remaining_unserved)
    soluT.unserved_users = remaining_unserved         
    if length(remaining_unserved) > 0
        # create a new route 
        used_vehs = soluT.RI[1:n_route, end]
        unused_vehs = setdiff!(collect(1:K), used_vehs)
        set_type_veh = parameter_energy.is_electric[used_vehs]
        if global_parameter.co2_reduc_target < 1 && (sum(set_type_veh) == n_route) # all EVs
            veh_id = findfirst(x->x==0, parameter_energy.is_electric[unused_vehs]) # use a disel veh
            vec_tmp =  [parameter_energy.cap_passenger[veh_id], veh_id]
            # @show(set_type_veh, vec_tmp )
        else
            vec_tmp =  soluT.RI[end, 4:5]
        end
        
        soluT.RI[n_route+2:end, :] = soluT.RI[n_route+1:end-1,: ] # shift a row
        soluT.RI[n_route+1,1:3] = [0,0,0]; soluT.RI[n_route+1,4:5] = vec_tmp
        soluT.n_route += 1
        
        route_i = soluT.n_route ; user = remaining_unserved[1]  
        soluT.vec_chg_event[route_i, :] .*= 0
        route = [start_depot, user, user+n_cus, end_depot]
        insert_consecutive(soluT, route_i, start_depot, user, end_depot, n_cus, qi) 
        update_auxiliary_var(soluT, darp, route_i, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back) # need to revise for flexbus problem
        setdiff!(soluT.unserved_users, user)
        for user in remaining_unserved[2:end]        
            inserted = greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, route_i, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back, flag_init_sol)
            inserted  && setdiff!(soluT.unserved_users, user)
        end 
        # # update charging
        # route = get_route(soluT, route_i, start_depot, end_depot)
        # n_charger_installed, success, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
        # if ! success && (n_charger_installed > 0)                            
        #     repair_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, pos_insert_chg, forward_time_slack, vec_t_start)
        # end
        # length(soluT.unserved_users)>0 &&  @show("regret_insert_unserved :", length(soluT.unserved_users), soluT.n_route, soluT.RI[1:soluT.n_route,end]) 
    end
end
  
# create a veh: if gasoline_fleet = true, create a gasoline veh, otherwise , create a electric veh of the first ev type
function create_route(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    
    n_type_veh, penalty = global_parameter.n_type_veh, global_parameter.penalty
    parameter_energy    = instance.parameter_energy
    veh_ids_types       = parameter_energy.veh_ids_types #set of veh_ids for each type
    qi, dist_all                  = darp.qi, darp.dist_all
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    ei, li, s_t, Li, dist, TH     = darp.ei, darp.li, darp.s_t, darp.Li, darp.dist, darp.TH
    update_sol(soluT, solution) #copy solution as temporary solution soluT  
    n_route = soluT.n_route 
    veh_id = 0
    if length(soluT.unserved_users) > 0
        # @show(" create_route call ! ", length(soluT.unserved_users))
        user = rand(collect(soluT.unserved_users)) # user denotes a bus stop
        route = [start_depot, user, user+n_cus, end_depot]
        if !eight_step(route, n_route+1, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)
            return false
        end         
        route_i = n_route+1 
        veh_used = soluT.RI[1:n_route, 5]
        if global_parameter.is_gasoline_fleet == true 
            veh_type_chosed = 1 
        else
            veh_type_chosed = rand(collect(1:n_type_veh))     
            if veh_type_chosed == 1 
                if soluT.co2_emission + length_route(route, darp) * parameter_energy.veh_info.CO2_type[1] > global_parameter.co2_threshold 
                    # @show(soluT.co2_emission,veh_used, parameter_energy.is_electric[veh_used], length_route(route, darp),length_route(route, darp) * parameter_energy.veh_info.CO2_type[1], global_parameter.co2_threshold  )
                    veh_type_chosed = rand(collect(2:n_type_veh))
                    # @show(veh_type_chosed)
                end
            end
        end 
        # @show(veh_type_chosed)
        vehs = setdiff(veh_ids_types[veh_type_chosed], veh_used)
        veh_id = vehs[1]
        # @show(vehs, veh_id, parameter_energy.is_electric[veh_id])
        soluT.RI[route_i, 4] = parameter_energy.cap_passenger[veh_id]
        soluT.RI[route_i, 5] = veh_id   
        soluT.n_route += 1  
        soluT.RI[route_i,3]= 0 # make it clean 
        insert_consecutive(soluT, route_i, start_depot, user, end_depot, n_cus, qi)  
        update_auxiliary_var(soluT, darp, route_i, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back) # need to revise for flexbus problem
        # @show(route_i, user, get_route(solution, route_i, start_depot, end_depot), "unserved_users (BE): ", solution.unserved_users)
        setdiff!(soluT.unserved_users, user)
        soluT.penalty_unserved -= qi[user]* penalty
        _, success,  _, _, _ = update_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
        if success 
            # @show(success)
            update_sol(solution, soluT)
            return true
        else 
            return false
        end
    end
    return false 
       
end

# destroy a random route and some part of other routes and re-insert them into the solution, a set of remove and insert heuristics are randomly paired
function destroy_repair(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)                                                         
    update_sol(soluT, solution) #copy solution as temporary solution soluT  
    n_cus  = darp.n_cus ; n_route = soluT.n_route
    max_degree_destruction = 0.275 # based on the tuned value from Roman Lutz (2014) Adaptive Large Neighborhood Search: A heuristic for the Rich Pickup
                                   # and Delivery Problem with Time Windows. Bachelor thesis at Ulm University.
    unserved_users =  soluT.unserved_users
    candidate_users_remove = setdiff!(collect(2:n_cus+1), unserved_users) # pay attention here, user is from 2:n_cus+1
     
    tmp = floor(Int, length(candidate_users_remove) * max_degree_destruction) 
    degree_destruct = min(60, tmp)
    degree_destruct == 0 && (degree_destruct += 1)
    n_remove = rand(collect(1:degree_destruct)) # num of requests to remove from the other routes
 
    remove_op= [rand_remove, worst_remove, dist_remove, tw_remove, shaw_remove] # remove operators
    vec_remove = collect(1:length(remove_op))
    selected_rm_op = remove_op[rand(vec_remove)]
    RMOP(selected_rm_op, soluT, darp, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
    insert_unserved(soluT, recorder_routes, recorder_lscost_routes, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    # RMOP(selected_rm_op, soluT, darp, lgraph, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
    # insert_unserved(solution, soluT,  recorder_routes, recorder_lscost_routes, fast_chg, global_parameter, instance, darp, lgraph, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, T, type_veh, cap_type_veh) 
    update_sol(solution, soluT) 
end


# update temporal solution
function update_tmp_sol(soluT::Solution, user, route_of_requests, r_i, vec_routes, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
    
    n_route = soluT.n_route
    n_cus = darp.n_cus
    soluT.pre[user] = soluT.succ[user] = soluT.pre[user+ n_cus] = soluT.succ[user+ n_cus] = -1

    if length(vec_routes[r_i])==2 # remove this route
       soluT.n_route -= 1
        for i in r_i+1:n_route
            soluT.RI[i-1,:] = soluT.RI[i,:]
            # soluT.TT_walk_routes[i-1]   = soluT.TT_walk_routes[i]
            soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
        end 
        for r in collect(r_i+1:n_route)      
            for node in vec_routes[r]
                route_of_requests[node] = r-1
            end  
        end
        deleteat!(vec_routes, r_i)
    else           
        route_auxiliary_update(soluT, darp, vec_routes[r_i], r_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)   
    end    
end

# insert the unserved customers according to their cost, based on Ropke and Pisinger (2016)
function greedy_insert_unserved(soluT::Solution, fast_chg, global_parameter, recorder_routes, recorder_lscost_routes, instance, darp, unserved_users, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back)
    
    # @show("greedy_insert_unserved call!")
    n_cus, start_depot, end_depot, qi=  darp.n_cus, darp.start_depot, darp.end_depot, darp.qi 
    dist_all, dist_orig, set_physical  = darp.dist_all, darp.dist_orig, darp.set_physical
    parameter_energy = instance.parameter_energy
    penalty = global_parameter.penalty
    T_max, K = instance.T_max, instance.K
    best_pos, best_route_temp= Int32[], Int32[]
    
    flag_init_sol = false; chg_succ = false 
    n_route = soluT.n_route # num of routes of the current solution
    routes= collect(1:n_route)
    
    route_r = 0; idx_min_cost = 0
    vec_cost_routes = zeros(n_route)  
    # @show(unserved_users, n_route)
    for user in unserved_users
        for (idx, _r) in enumerate(routes)
            vec_cost_routes[idx], best_pos, best_route_temp = cost_greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, flag_init_sol) 
            recorder_routes[idx].r, recorder_routes[idx].best_pos, recorder_routes[idx].new_route = _r, best_pos, best_route_temp
        end
        # @show(user, n_route, vec_cost_routes)
        vec_idx_sorted_cost= sortperm(vec_cost_routes) # sorted index set
        idx_min_cost = vec_idx_sorted_cost[1]          
        recorder_lscost_routes[user].insert_cost = vec_cost_routes[idx_min_cost] # store the best position and route
        recorder_lscost_routes[user].r           = recorder_routes[idx_min_cost].r 
        recorder_lscost_routes[user].best_pos    = recorder_routes[idx_min_cost].best_pos 
        recorder_lscost_routes[user].new_route   = recorder_routes[idx_min_cost].new_route    
    end
    
    vec_insert_cost= [recorder_lscost_routes[user].insert_cost for user in unserved_users]
    idx_vec_insert_cost = sortperm(vec_insert_cost)
    sorted_unserved_users = unserved_users[idx_vec_insert_cost]    
    # @show(unserved_users, sorted_unserved_users, length(unserved_users), length(sorted_unserved_users))
    flag_route_changed = zeros(Bool, n_route)
    copy_unserved = Int32[]
    for user in sorted_unserved_users
        route_r = recorder_lscost_routes[user].r  
        veh_id = soluT.RI[route_r,5]
        best_pos, best_route_temp = recorder_lscost_routes[user].best_pos, recorder_lscost_routes[user].new_route
        # @show(user,route_r, veh_id, best_pos, best_route_temp, flag_route_changed[route_r] )
        if length(best_pos) == 0 # no feasible route
            push!(copy_unserved, user)
        else
            if ! flag_route_changed[route_r] 
                if parameter_energy.is_electric[veh_id] == true # EV 
                    compute_dt_r(best_route_temp, dt_r, dist_orig)
                    ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
                    if ec_route > parameter_energy.max_ec[veh_id]
                        if fast_chg.n_fast_chg_installed > 0
                            if insert_charging_route(soluT, darp, route_r, best_route_temp, dist_all, ec_route, fast_chg, parameter_energy, dt_r, dt_r_back, flag_init_sol)[1]
                                insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                                update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem 
                                flag_route_changed[route_r] = true
                            else
                                push!(copy_unserved, user)
                            end 
                        else
                            push!(copy_unserved, user)
                        end
                    else    
                        soluT.vec_chg_event[route_r,:] .*= 0           
                        insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                        update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem  
                        flag_route_changed[route_r] = true
                    end
                else
                    soluT.vec_chg_event[route_r,:] .*= 0 
                    insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                    update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
                    flag_route_changed[route_r] = true
                end  
            else  # route has been changed due to inserting precedent unserved customers, so we insert the customer to the best pos of this route
                if !greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, route_r, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back, flag_init_sol)
                    push!(copy_unserved, user)
                end
            end
        end
    end
    
    # try again 
    # @show("try again : ", copy_unserved)
    remaining_unserved = Int32[]
    for user in copy_unserved
        success = greedy_insert_request(soluT, fast_chg, global_parameter, recorder_routes, instance, darp, user, 0, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back)
        if ! success
            push!(remaining_unserved, user)
        end
    end
    # @show(remaining_unserved)
    soluT.unserved_users = remaining_unserved 
    if length(remaining_unserved) > 0
        # create a new route 
        used_vehs = soluT.RI[1:n_route, end]
        unused_vehs = setdiff!(collect(1:K), used_vehs)
        set_type_veh = parameter_energy.is_electric[used_vehs]
        if global_parameter.co2_reduc_target < 1 && (sum(set_type_veh) == n_route) # all EVs
           veh_id = findfirst(x->x==0, parameter_energy.is_electric[unused_vehs]) # use a disel veh
           vec_tmp =  [parameter_energy.cap_passenger[veh_id], veh_id] 
        else
           vec_tmp =  soluT.RI[end, 4:5]
        end

        soluT.RI[n_route+2:end, :] = soluT.RI[n_route+1:end-1,: ] # shift a row
        soluT.RI[n_route+1,1:3] = [0,0,0]; soluT.RI[n_route+1,4:5] = vec_tmp
        soluT.n_route += 1
        route_i = soluT.n_route ; user = remaining_unserved[1]  
        soluT.vec_chg_event[route_i, :] .*= 0
        route = [start_depot, user, user+n_cus, end_depot]
        insert_consecutive(soluT, route_i, start_depot, user, end_depot, n_cus, qi) 
        update_auxiliary_var(soluT, darp, route_i, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back) # need to revise for flexbus problem
        setdiff!(soluT.unserved_users, user)
        for user in remaining_unserved[2:end]        
            inserted = greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, route_i, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back, flag_init_sol)
            inserted  && setdiff!(soluT.unserved_users, user)
        end 
        # update charging
        # route = get_route(soluT, route_i, start_depot, end_depot)
        # n_charger_installed, success, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
        # if ! success && (n_charger_installed > 0)                            
        #     repair_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, pos_insert_chg, forward_time_slack, vec_t_start)
        # end
        # length(soluT.unserved_users)>0 &&  @show("greedy_insert_unserved : ", length(soluT.unserved_users), soluT.n_route, soluT.RI[1:soluT.n_route,end])
    end
    
end


# insert the request to the least cost position of the other routes, route_i is the route of the request "user", route_i = 0 if it is the unserved pool
function greedy_insert_request(soluT::Solution, fast_chg, global_parameter, recorder_routes, instance, darp, user, route_i, e_r, cap_r, l_r, maxcap_r,dt_r,dt_r_back)
    
    start_depot, end_depot, qi= darp.start_depot, darp.end_depot, darp.qi 
    dist_all, dist_orig  = darp.dist_all, darp.dist_orig
    parameter_energy =  instance.parameter_energy  
    bigM =  global_parameter.bigM 
    
    best_pos, best_route_temp= Int32[], Int32[]
  
    flag_init_sol = false 
    n_route = soluT.n_route # num of candidate routes
    routes= collect(1:n_route)
    if route_i > 0 
        setdiff!(routes, route_i)  
        n_route -= 1
    end 
    route_r = 0; idx_min_cost = 0
    vec_cost_routes = zeros(n_route) # cost_increase, route_i, best_position
     
    for (idx, _r) in enumerate(routes)
        vec_cost_routes[idx], best_pos, best_route_temp = cost_greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, flag_init_sol) 
        recorder_routes[idx].r, recorder_routes[idx].best_pos, recorder_routes[idx].new_route = _r, best_pos, best_route_temp
    end
    vec_idx_sorted_cost= sortperm(vec_cost_routes) # sorted index set
    for idx_min_cost in vec_idx_sorted_cost 
        if vec_cost_routes[idx_min_cost] < bigM   
            route_r = recorder_routes[idx_min_cost].r  # route to insert the rejected request
            veh_id = soluT.RI[route_r,5]
            best_pos, best_route_temp = recorder_routes[idx_min_cost].best_pos, recorder_routes[idx_min_cost].new_route
            if parameter_energy.is_electric[veh_id] == true # EV 
                compute_dt_r(best_route_temp, dt_r, dist_orig)
                ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
                if ec_route > parameter_energy.max_ec[veh_id]
                    if fast_chg.n_fast_chg_installed > 0
                        if insert_charging_route(soluT, darp, route_r, best_route_temp, dist_all, ec_route, fast_chg, parameter_energy, dt_r, dt_r_back, flag_init_sol)[1]
                            insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                            update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem 
                            return true 
                        end 
                    end
                else    
                    soluT.vec_chg_event[route_r,:] .*= 0           
                    insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                    update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem  
                    return true
                end
            else
                soluT.vec_chg_event[route_r,:] .*= 0 
                insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
                return true
            end  
        end  
    end 
    return false 
end

# "relocate operator removes a user request from its current route and tries to reinsert it in the least cost position of the other routes"
# revised on 1.9.2023
function relocate(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)

    n_route = solution.n_route
    if n_route < 2 
        return false
    end
    
    n_cus, start_depot, end_depot, qi= darp.n_cus, darp.start_depot, darp.end_depot, darp.qi 
    dist_all, set_physical, dist_orig  = darp.dist_all, darp.set_physical, darp.dist_orig
    K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy  
    bigM =  global_parameter.bigM 
    
    best_pos, best_route_temp= Int32[], Int32[]
    
    flag_init_sol = false
    chg_succ = false
    update_sol(soluT,solution) #copy solution as temporary solution soluT
    route_i=0; routes=Int32[]; users=Int32[]    
    
    route_i = rand(1:n_route)                     # selected route to eject a request
    routes = setdiff!(collect(1:n_route),route_i) # other routes
    
    route = get_route(soluT, route_i, start_depot, end_depot)
    users = route[findall(x->1<x<n_cus+2, route)]
    shuffle!(users) 
    # @show(users)
    user = 0;  route_r = 0
    vec_cost_routes = zeros(n_route-1) # cost_increase, route_i, best_position
    done =false
    for user in users  
        for (idx, _r) in enumerate(routes)
            vec_cost_routes[idx], best_pos, best_route_temp = cost_greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, flag_init_sol) 
            recorder_routes[idx].r, recorder_routes[idx].best_pos, recorder_routes[idx].new_route = _r, best_pos, best_route_temp
        end
        idx_min_cost= sortperm(vec_cost_routes)[1]
        
        if vec_cost_routes[idx_min_cost] < bigM
            done=true
            route_r = recorder_routes[idx_min_cost].r  # route to insert the rejected request
            veh_id = solution.RI[route_r,5]
            best_pos, best_route_temp = recorder_routes[idx_min_cost].best_pos, recorder_routes[idx_min_cost].new_route
            if parameter_energy.is_electric[veh_id] == true # EV 
                compute_dt_r(best_route_temp, dt_r, dist_orig)
                ec_route = dt_r[end_depot] * parameter_energy.beta[veh_id] # energy consumption of the route
                if ec_route > parameter_energy.max_ec[veh_id]
                    # @show(veh_id, ec_route, parameter_energy.max_ec[veh_id], get_route(soluT, route_r, start_depot, end_depot))
                    if fast_chg.n_fast_chg_installed == 0
                        done = false
                    else
                        if insert_charging_route(soluT, darp, route_r, best_route_temp, dist_all, ec_route, fast_chg, parameter_energy, dt_r, dt_r_back, flag_init_sol)[1]
                            insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                            update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem 
                        else
                            done = false
                        end 
                    end
                else    
                    soluT.vec_chg_event[route_r,:] .*= 0           
                    insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                    update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem  
                end             
            else
                soluT.vec_chg_event[route_r,:] .*= 0 
                insert_user_route(soluT, route_r, best_pos, start_depot, end_depot, qi)
                update_auxiliary_var(soluT, darp, route_r, best_pos, e_r,cap_r,l_r,maxcap_r,dt_r,dt_r_back) # need to revise for flexbus problem
            end       
        end  
        if done              
            if length(route)==4 # remove this route 
                soluT.RI[end,4:5] =  soluT.RI[route_i,4:5] # move this unused veh (capacity and veh_id) to the end of the list
                soluT.n_route -= 1
                for i in route_i+1:n_route
                    soluT.RI[i-1,:] = soluT.RI[i,:] 
                    soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
                end      
                update_sol(solution, soluT) 
                # verify_solution(solution, instance, darp, fast_chg)              
                return true 
            else #update this route          
                filter!(x->x∉[user, user + n_cus], route) 
                route_auxiliary_update(soluT, darp, route , route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)  
                if soluT.vec_chg_event[route_i, 1] > 0  
                    _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
                    if  chg_succ
                        update_sol(solution, soluT)                               
                    else
                        repair_charging(soluT, darp, route_i, route, fast_chg, parameter_energy, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start) && update_sol(solution, soluT)
                    end
                end  
                return true 
            end        
        end  
    end
    return false  # means that this operator has been done and count for performance score update     
end

#repair an infeasible route in a greed way
function repair(solution::Solution,soluT::Solution, demand_data, fast_chg, global_parameter, instance, darp, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back,  Q, Li)
    
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    dist_all, set_physical = darp.dist_all, darp.set_physical 
    K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy

    flag_init_sol= false 
    update_sol(soluT,solution) #copy solution as temporary solution soluT
    n_route = soluT.n_route        
    route = get_route(soluT, route_i, start_depot, end_depot) 
    users= route[findall(x->1<x<n_cus+2, route)]
    routes = setdiff!(collect(1:n_route), route_i)
    user = shuffle!(users)
    for user in users
        shuffle!(routes)    
        done = false 
        for _r in routes
            if greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back,flag_init_sol) 
                done = true
                break
            end
        end
        if done 
            if length(route)==4  # remove this route
                soluT.RI[end,4:5] =  soluT.RI[route_i,4:5] # move this unused veh (capacity and veh_id) to the end of the list
                soluT.n_route -= 1
                for i in route_i+1:n_route
                    soluT.RI[i-1,:] = soluT.RI[i,:] 
                    soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
                end   
                update_sol(solution, soluT)               
                return true 
            else #update this route           
                filter!(x->x∉[user, user + n_cus], route)
                # soluT.TT_walk_routes[route_i] -= walk_tt_bus_stop[user-1]
                route_auxiliary_update(soluT, darp, route , route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)  
                # if eight_step(route, route_i, ei, li, s_t, dist_orig, Q, Li, n_cus, TH)
                    if soluT.vec_chg_event[route_i, 1] > 0 
                        _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
                        if  chg_succ
                            update_sol(solution, soluT)  
                        else
                            repair_charging(soluT, darp, route_i, route, fast_chg, parameter_energy, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start) && update_sol(solution, soluT) 
                        end
                    end
                    return true
                # end
            end
        end 
    end   
    return false  
end


# #swap two segments of two different routes
function swap_seg(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back,  Q, T, Li)
    
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    ei, li, s_t, qi, dist, dist_orig=  darp.ei, darp.li, darp.s_t, darp.qi, darp.dist, darp.dist_orig
    dist_all, set_physical  = darp.dist_all, darp.set_physical 
    T_max, parameter_energy, K=  instance.T_max, instance.parameter_energy, instance.K
    TH = darp.TH

    update_sol(soluT, solution) #copy solution as temporary solution soluT
    chg_succ_i=chg_succ_j=false
    n_route= soluT.n_route
    if n_route < 2
        return false
    end
    si=Int32[];sj=Int32[]
    nodes_route1=Int32[];nodes_route2=Int32[]
    ri_new = zeros(Int32,2*n_cus); rj_new = zeros(Int32,2*n_cus)
    
    r_i ,r_j = shuffle!(collect(1:n_route))[1:2]    
    
    route_i  = get_route(soluT,r_i, start_depot, end_depot)
    reset_cap_r(soluT, r_i, cap_r, darp)
    cap_r_ri = [cap_r[vi] for vi in route_i] 
    idx_ri   = findall(x->x==0, cap_r_ri) # find location of zero load    
    route_j  = get_route(soluT, r_j, start_depot, end_depot)
    reset_cap_r(soluT, r_j,  cap_r, darp)
    cap_r_rj = [cap_r[vi] for vi in route_j] 
    idx_rj   = findall(x->x==0, cap_r_rj)
    pre_idx1 = idx_ri[1]
    n_ri_new = 0; n_rj_new = 0
    eh_vi=0; eh_vj=0
    for idx1 in idx_ri[2:end-1]
        # @show(pre_idx1,idx1)
        si = route_i[pre_idx1+1:idx1]
        pre_si = route_i[pre_idx1]
        suc_si = route_i[idx1+1]
        n_si = length(si)
        pre_idx2 = idx_rj[1]
        for idx2 in idx_rj[2:end-1]
            # @show(pre_idx2,idx2) 
            sj = route_j[pre_idx2+1:idx2]
            # if  is_compatible(si[1], sj[1], layer_nodes, lyrs_compatible, end_depot)
                pre_sj = route_j[pre_idx2]
                suc_sj = route_j[idx2+1]
                n_sj = length(sj)
                dist_saving = dist[pre_si,si[1]]+ dist[si[end],suc_si]+dist[pre_sj,sj[1]]+ dist[sj[end],suc_sj]
                - (dist[pre_si,sj[1]]+ dist[sj[end],suc_si] + dist[pre_sj,si[1]] + dist[si[end],suc_sj])
                if dist_saving > 0 
                    eh_si = max(ei[si[1]], e_r[pre_sj]+s_t[pre_sj]+dist[pre_sj,si[1]])
                    eh_sj = max(ei[sj[1]], e_r[pre_si]+s_t[pre_si]+dist[pre_si,sj[1]])
                    if eh_si <= li[si[1]] && eh_sj <= li[sj[1]] 
                        eh_pre = eh_si ; pre_v = si[1]
                        feasible = true
                        for v in si[2:end]  
                            eh_vi = max(ei[v], eh_pre+s_t[pre_v]+dist[pre_v,v])
                            if eh_vi > li[v]  
                                feasible = false
                                break
                            end
                            eh_pre = eh_vi ; pre_v = v
                        end  
                        if feasible
                            eh_pre = eh_sj ; pre_v = sj[1]
                            for v in sj[2:end]  
                                eh_vj = max(ei[v], eh_pre+s_t[pre_v]+dist[pre_v,v])
                                if eh_vj > li[v]  
                                    feasible = false
                                    break
                                end
                                eh_pre = eh_vj ; pre_v = v
                            end                        
                            if feasible
                                eh_suc_si = max(ei[suc_si], eh_vj+s_t[sj[end]]+dist[sj[end],suc_si])
                                eh_suc_sj = max(ei[suc_sj], eh_vi+s_t[si[end]]+dist[si[end],suc_sj])
                                if eh_suc_si <= li[suc_si] && eh_suc_sj <= li[suc_sj]
                                    n_ri_new = pre_idx1 + n_sj + length(route_i) - (pre_idx1+n_si)
                                    n_rj_new = pre_idx2 + n_si + length(route_j) - (pre_idx2+n_sj)
                                    ri_new[1:pre_idx1], ri_new[pre_idx1+1: pre_idx1+n_sj], ri_new[pre_idx1+n_sj+1:n_ri_new] = route_i[1:pre_idx1], sj, route_i[idx1+1:end]
                                    rj_new[1:pre_idx2], rj_new[pre_idx2+1: pre_idx2+n_si], rj_new[pre_idx2+n_si+1:n_rj_new] = route_j[1:pre_idx2], si, route_j[idx2+1:end]
                                    # @show(ri_new[1:n_ri_new],rj_new[1:n_rj_new])
                                    # nodes_route1 = nodes[ri_new[1:n_ri_new]];nodes_route2 = nodes[rj_new[1:n_rj_new]]
                                    # if  no_double_tour(nodes_route1, D′) &&  no_double_tour(nodes_route2, D′)
                                        if eight_step(ri_new[1:n_ri_new],r_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi) && eight_step(rj_new[1:n_rj_new],r_j, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)
                                            # @show("eight_step successfully !!")
                                            route_update(soluT, ri_new[1:n_ri_new], r_i, n_cus)
                                            route_update(soluT, rj_new[1:n_rj_new], r_j, n_cus)    
                                            update_auxiliary_var(soluT, darp, r_i, [pre_si, sj[1], sj[2], sj[end-1], sj[end], suc_si],e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back) 
                                            update_auxiliary_var(soluT, darp, r_j, [pre_sj, si[1], si[2], si[end-1], si[end], suc_sj],e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
                                            n_charger_installed, chg_succ_i, pos_insert_chg_i, forward_time_slack_i, vec_t_start_i = update_charging(soluT, darp, r_i, ri_new[1:n_ri_new], dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
                                            n_charger_installed, chg_succ_j, pos_insert_chg_j, forward_time_slack_j, vec_t_start_j = update_charging(soluT, darp, r_j, rj_new[1:n_rj_new], dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
                                            if  chg_succ_i && chg_succ_j
                                                update_sol(solution, soluT)  
                                            else     
                                                if  n_charger_installed > 0                            
                                                    ! chg_succ_j && (chg_succ_j = repair_charging(soluT, darp, r_j, rj_new[1:n_rj_new], fast_chg, parameter_energy, dt_r, pos_insert_chg_j, forward_time_slack_j, vec_t_start_j))   
                                                    compute_dt_r(ri_new[1:n_ri_new], dt_r, dist_orig)   
                                                    ! chg_succ_i && (chg_succ_i = repair_charging(soluT, darp, r_i, ri_new[1:n_ri_new], fast_chg, parameter_energy, dt_r, pos_insert_chg_i, forward_time_slack_i, vec_t_start_i))                                  
                                                    if chg_succ_i && chg_succ_j
                                                        update_sol(solution, soluT)   
                                                    end
                                                end 
                                            end
                                            return true
                                        end
                                    # end
                                end
                            end
                        end
                    end
                end
            pre_idx2 = idx2
        end
        pre_idx1 = idx1
    end
    return false  # means that this operator has been done and count for performance score update
end

# exchange two segments
# "only arcs for which the vehicle is empty are considered to be removed
function two_opt_star(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    TH, ei, li, s_t, dist, dist_orig, qi =  darp.TH, darp.ei, darp.li, darp.s_t, darp.dist, darp.dist_orig, darp.qi
    dist_all, set_physical  = darp.dist_all, darp.set_physical 
    K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy  

    update_sol(soluT, solution) #copy solution as temporary solution soluT
    chg_succ_i = chg_succ_j = false
    n_route= soluT.n_route
    if n_route < 2
        return false
    end
    r_i = 0;r_j=0
    route_i = Int32[] ;route_j =Int32[];cap_r_ri=Int32[];cap_r_rj=Int32[]
    idx_ri  = Int32[] ;idx_rj = Int32[];seg_i1=Int32[];n_segi1=0;n_segi2=0;n_segj1=0;n_segj2=0
    # nodes_route1=Int32[] ;nodes_route2=Int32[] 
    ri_new = zeros(Int32,2*n_cus); rj_new = zeros(Int32,2*n_cus)  
    
    rand_start_ri = rand_start_rj=2
    r_i ,r_j = shuffle!(collect(1:n_route))[1:2]    
    route_i  = get_route(soluT,r_i, start_depot, end_depot)
    reset_cap_r(soluT, r_i, cap_r, darp)
    cap_r_ri = [cap_r[vi] for vi in route_i] 
    idx_ri   = findall(x->x==0, cap_r_ri) # find location of zero load   
    n_ri = length(idx_ri) 
    n_ri>3 && (rand_start_ri = rand(collect(2:n_ri-1)))
    
    route_j  = get_route(soluT, r_j, start_depot, end_depot)
    reset_cap_r(soluT, r_j,  cap_r, darp)
    cap_r_rj = [cap_r[vi] for vi in route_j] 
    idx_rj   = findall(x->x==0, cap_r_rj)
    n_rj = length(idx_rj) 
    n_rj>3 && (rand_start_rj = rand(collect(2:n_rj-1)))
    
    for idx1 in idx_ri[rand_start_ri:end-1]
        vi = route_i[idx1]
        pre_vi = soluT.pre[vi]
        suc_vi = soluT.succ[vi]
        seg_i1=route_i[1:idx1]  
        seg_i2=route_i[idx1+1:end]  
        n_segi1=idx1; n_segi2= length(seg_i2)  
        # @show(pre_vi,vi,suc_vi,seg_i1,seg_i2)
        for idx2 in idx_rj[rand_start_rj:end-1]
            vj =route_j[idx2]
            #check dist savings
            pre_vj = soluT.pre[vj]
            suc_vj = soluT.succ[vj]
            # @show(pre_vj,vj,suc_vj)
            if  dist[vi,suc_vi]+dist[vj,suc_vj] - (dist[vi,suc_vj] + dist[vj,suc_vi]) > T #T is the temparature of the SA 
                #check feasibility
                eh_suc_vi = max(ei[suc_vi], e_r[vj]+s_t[vj]+dist[vj,suc_vi])
                eh_suc_vj = max(ei[suc_vj], e_r[vi]+s_t[vi]+dist[vi,suc_vj])
                if eh_suc_vi <= li[suc_vi] && eh_suc_vj <= li[suc_vj] 
                    seg_j1= route_j[1:idx2] 
                    seg_j2= route_j[idx2+1:end]  
                    n_segj1=idx2; n_segj2 = length(seg_j2) 
                    n_ri_new=n_segi1+n_segj2 ;n_rj_new= n_segj1+n_segi2
                    # @show(seg_j1,seg_j2) #avoid using vcat(start_depot,seg_i1,seg_j2, end_depot): it is much slower
                    ri_new[1:n_segi1], ri_new[n_segi1+1:n_segi1+n_segj2] = seg_i1, seg_j2
                    rj_new[1:n_segj1], rj_new[n_segj1+1:n_segj1+n_segi2] = seg_j1, seg_i2
                    # @show(ri_new[1:n_ri_new],rj_new[1:n_rj_new])
                    # nodes_route1 = ri_new[1:n_ri_new];nodes_route2 = rj_new[1:n_rj_new]
                    # if  no_double_tour(nodes_route1, D′) &&  no_double_tour(nodes_route2, D′)  
                        if eight_step(ri_new[1:n_ri_new],r_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi) && eight_step(rj_new[1:n_rj_new],r_j, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)
                            # @show("eight_step ok")
                            route_update(soluT, ri_new[1:n_ri_new], r_i, n_cus)
                            route_update(soluT, rj_new[1:n_rj_new],r_j, n_cus)  
                            vec_tmp_i=[pre_vi, vi, suc_vj, pre_vi, vi, suc_vj]
                            vec_tmp_j=[pre_vj, vj, suc_vi, pre_vj, vj, suc_vi]
                            update_auxiliary_var(soluT, darp, r_i, vec_tmp_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back) 
                            update_auxiliary_var(soluT, darp, r_j, vec_tmp_j,e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
                            n_charger_installed, chg_succ_i, pos_insert_chg_i, forward_time_slack_i, vec_t_start_i = update_charging(soluT, darp, r_i, ri_new[1:n_ri_new], dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
                            n_charger_installed, chg_succ_j, pos_insert_chg_j, forward_time_slack_j, vec_t_start_j = update_charging(soluT, darp, r_j, rj_new[1:n_rj_new], dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)                               
                            if  chg_succ_i && chg_succ_j
                                update_sol(solution, soluT)  
                            else                               
                                if  n_charger_installed > 0 
                                    ! chg_succ_j && (chg_succ_j = repair_charging(soluT, darp, r_j, rj_new[1:n_rj_new], fast_chg, parameter_energy, dt_r, pos_insert_chg_j, forward_time_slack_j, vec_t_start_j))                                   
                                    compute_dt_r(ri_new[1:n_ri_new], dt_r, dist_orig)  
                                    ! chg_succ_i && (chg_succ_i = repair_charging(soluT, darp, r_i, ri_new[1:n_ri_new], fast_chg, parameter_energy, dt_r, pos_insert_chg_i, forward_time_slack_i, vec_t_start_i))                                   
                                    if chg_succ_i && chg_succ_j
                                        update_sol(solution, soluT)   
                                    end
                                end                                 
                            end
                            return true
                        end
                    # end
                end
            end
        end
    end    
    return false  # means that this operator has been done and count for performance score update 
end

 
function swap_two_users(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    TH, ei, li, s_t, dist, dist_orig, qi =  darp.TH, darp.ei, darp.li, darp.s_t, darp.dist, darp.dist_orig, darp.qi
    dist_all, set_physical  = darp.dist_all, darp.set_physical 
    T_max, parameter_energy, K = instance.T_max, instance.parameter_energy, instance.K 
    
    flag_init_sol = false 
    update_sol(soluT, solution) #copy solution as temporary solution soluT
    # insert from route_i to route j
    init_route_info(e_r,cap_r,l_r,maxcap_r, li, start_depot, end_depot, dt_r, dt_r_back)
    n_route = soluT.n_route
    if n_route < 2
        return false
    end
    r_i ,r_j = shuffle!(collect(1:n_route))[1:2]
    
    seg_i1=Int32[];seg_i2=Int32[];seg_j1=Int32[];seg_j2=Int32[];copy_seg_j2 =Int32[];copy_seg_i2 =Int32[]
    route_j_new = zeros(Int32,2*n_cus+2); route_i_new = zeros(Int32,2*n_cus+2); new_route = zeros(Int32,4)
    nodes_route = Int32[]; route_j_tmp = Int32[]
    idx3=0 ;  dist_saving = 0; route =Int32[]
    if soluT.RI[r_i,3]>3
        route_i = get_route(soluT, r_i, start_depot, end_depot)
        route_j = get_route(soluT, r_j, start_depot, end_depot)
        # insert from route_i to route j
        idx_users_ri = findall(x->x<n_cus+2, route_i)
        idx_users_rj = findall(x->x<n_cus+2, route_j)
        for idx in 2:length(idx_users_ri)
            idx1 = idx_users_ri[idx]
            vi   = route_i[idx1]
            for _idx in 2:length(idx_users_rj)
                idx2 = idx_users_rj[_idx]
                vj = route_j[idx2]
                # if lyrs_compatible[lyrs_nodes[vi],lyrs_nodes[vj]]
                pre_vi = route_i[idx1-1]
                suc_vi = route_i[idx1+1]
                pre_vj = route_j[idx2-1]
                suc_vj = route_j[idx2+1]
                v_ni = vi+n_cus
                v_nj = vj+n_cus
                seg_i1=route_i[1:idx1-1];seg_i2=route_i[idx1+1:end]
                seg_j1=route_j[1:idx2-1];seg_j2=route_j[idx2+1:end]  
                pre_vj=soluT.pre[vj]  
                pre_v_ni = soluT.pre[v_ni]
                suc_v_ni = soluT.succ[v_ni] 
                pre_v_nj = soluT.pre[v_nj]
                suc_v_nj = soluT.succ[v_nj]         
                dist_saving = dist[pre_vi,vi]+dist[vi,suc_vi]+dist[pre_vj,vj]+dist[vj,suc_vj]
                + dist[pre_v_ni,v_ni]+dist[v_ni,suc_v_ni]+dist[pre_v_nj,v_nj]+dist[v_nj,suc_v_nj]
                - (dist[pre_vj,vi]+dist[vi,suc_vj]+dist[pre_v_nj,v_ni]+dist[v_ni,suc_v_nj] 
                + dist[pre_vi,vj]+dist[vj,suc_vi]+dist[pre_v_ni,v_nj]+dist[v_nj,suc_v_ni])
                if dist_saving > T
                    # check tw feasibility if inserting vi at the position of vj
                    eh_vi = max(ei[vi], e_r[pre_vj]+s_t[pre_vj]+dist[pre_vj,vi])
                    if eh_vi <= li[vi]  
                        # @show(route_j,seg_j1,vj,seg_j2)
                        idx3 = findfirst(x->x==v_nj, seg_j2)
                        feasible = true
                        pre_v = vi 
                        eh_pre_v = eh_vi
                        for i in 1:idx3-1
                            next = seg_j2[i]
                            eh_pre_v = max(ei[next], eh_pre_v + s_t[pre_v]+dist[pre_v,next])
                            if eh_pre_v > li[next]  
                                feasible = false
                                break
                            end
                            pre_v = next
                        end
                        if feasible 
                            eh_v_ni    = max(ei[v_ni], eh_pre_v + s_t[pre_v]+dist[pre_v,v_ni])
                            eh_suc_v_nj= max(ei[suc_v_nj], eh_v_ni + s_t[v_ni]+dist[v_ni,suc_v_nj]) 
                            if eh_v_ni<=li[v_ni] && eh_suc_v_nj <= li[suc_v_nj]
                                if eh_v_ni - li[vi] - s_t[vi] <= Li[vi] 
                                    copy_seg_j2 =copy(seg_j2) 
                                    copy_seg_j2[idx3] = v_ni
                                    n1=length(seg_j1); n2=length(copy_seg_j2)
                                    route_j_new[1:n1],route_j_new[n1+1] = seg_j1,vi
                                    route_j_new[n1+2:n1+n2+1] = copy_seg_j2 
                                    route_j_tmp = route_j_new[1:n1+n2+1]
                                    if eight_step(route_j_tmp, r_j, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)                                           
                                        route_auxiliary_update(soluT, darp, route_j_tmp, r_j, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
                                        n_charger_installed, chg_succ_j, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, r_j, route_j_tmp, dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
                                        if  ! chg_succ_j && (n_charger_installed > 0)
                                            chg_succ_j = repair_charging(soluT, darp, r_j, route_j_tmp, fast_chg, parameter_energy, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start)                                  
                                        end
                                        if chg_succ_j
                                            copy_seg_i2 = copy(seg_i2)
                                            deleteat!(copy_seg_i2,  findfirst(x->x==v_ni, copy_seg_i2));
                                            n1=length(seg_i1) ; n2 =length(copy_seg_i2)
                                            route_i_new[1:n1], route_i_new[n1+1:n1+n2]=  seg_i1,copy_seg_i2
                                            route_i_tmp =  route_i_new[1:n1+n2]
                                            route_auxiliary_update(soluT, darp, route_i_tmp, r_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
                                            if greedy_insert(soluT, fast_chg, global_parameter, instance, darp, vj, r_i, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, flag_init_sol)
                                                update_sol(solution, soluT)
                                                return true
                                            else  
                                                routes = setdiff!(collect(1:n_route), [r_i,r_j])                              
                                                if length(routes)>0
                                                    shuffle!(routes) 
                                                    for _r in routes
                                                        if greedy_insert(soluT, fast_chg, global_parameter, instance, darp, vj, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, flag_init_sol)                        
                                                            route = get_route(soluT, r_i, start_depot, end_depot)
                                                            if eight_step(route, r_i, ei, li, s_t, dist_orig, Q, Li, n_cus, TH, qi)
                                                                _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, r_i, route, dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
                                                                if  chg_succ                                                                 
                                                                    update_sol(solution, soluT)  
                                                                else
                                                                    if  repair_charging(soluT, darp, r_i, route, fast_chg, parameter_energy, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start) 
                                                                        update_sol(solution, soluT)    
                                                                    end                   
                                                                end
                                                            else
                                                                if repair(soluT, soluT, demand_data, fast_chg, global_parameter, instance, darp, r_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, Li) 
                                                                    update_sol(solution, soluT)  
                                                                end
                                                            end
                                                            return true
                                                        end
                                                    end                                                 
                                                    return true                    
                                                else     
                                                    return true
                                                end
                                            end
                                        end
                                    end
                                end                        
                            end
                        end
                    end
                end
                # end
            end
        end
    end 
    return false
end


function four_opt_intra(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)

    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    TH, ei, li, s_t, dist, qi=  darp.TH, darp.ei, darp.li, darp.s_t, darp.dist, darp.qi

    init_route_info(e_r,cap_r,l_r,maxcap_r, li, start_depot, end_depot, dt_r, dt_r_back)
    n_route = solution.n_route
    if n_route < 1
        return false
    end
    route_i = shuffle(collect(1:n_route))[1]
    n_nodes = solution.RI[route_i,3]
    combis =  collect(permutations(collect(1:3)))[2:end]
    idx=0; best_route = Int32[];seg =Int32[]; nodes_route = Int32[]; nodes_tmp =Int32[]
    if n_nodes > 3 
        route = get_route(solution, route_i, start_depot, end_depot)
        # @show(route_i,route)
        idx =2 # 1 is the depot
        pre_v = route[1]
        n_iter = floor(Int32, n_nodes/3)
        for iter in 1:n_iter
            done = false
            eh_prev = e_r[pre_v]
            nodes_tmp = route[idx:idx+2]
            suc_v = route[idx+3]
            dist_0 = dist[pre_v,nodes_tmp[1]]+dist[nodes_tmp[1],nodes_tmp[2]]+dist[nodes_tmp[2],nodes_tmp[3]]+dist[nodes_tmp[3],suc_v]
            best_dist_saving = 0
            for comb in combis
                seg = copy(nodes_tmp[comb])
                dist_saving = dist_0 - (dist[pre_v,seg[1]]+dist[seg[1],seg[2]]+dist[seg[2],seg[3]]+dist[seg[3],suc_v])
                if dist_saving > 0
                    # check order feasibility   
                    if check_pre_seg(seg, n_cus)
                        #check tw
                        eh_vi_1  = max(ei[seg[1]], eh_prev+s_t[pre_v]+dist[pre_v,seg[1]])
                        eh_vi_2  = max(ei[seg[2]], eh_vi_1+s_t[seg[1]]+dist[seg[1],seg[2]])
                        eh_vi_3  = max(ei[seg[3]], eh_vi_2+s_t[seg[2]]+dist[seg[2],seg[3]])
                        eh_suc_v = max(ei[suc_v],  eh_vi_3+s_t[seg[3]]+dist[seg[3],suc_v])
                        check =[eh_vi_1,eh_vi_2,eh_vi_3,eh_suc_v]-[li[seg[1]],li[seg[2]],li[seg[3]],li[suc_v]]
                        if maximum(check) <= 0 # tw feasible 
                            route_temp = vcat(route[1:idx-1],seg,route[idx+3:end])
                            # nodes_route = nodes[route_temp[2:end-1]]
                            if dist_saving > best_dist_saving
                                # if no_double_tour(nodes_route, D′) && eight_step(route_temp,route_i, ei, li, s_t, dist, Q, Li, n_cus, TH) 
                                if eight_step(route_temp,route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi) 
                                    done = true
                                    best_route = copy(route_temp)
                                    best_dist_saving = dist_saving
                                end
                            end
                        end
                    end
                end
            end
            # update solution
            done && route_auxiliary_update(solution, darp, best_route, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)   
            # verify_solution(solution)            
            idx += 3; pre_v = nodes_tmp[end]   
        end
    end
    return false  # means that this operator has been done and count for performance score update  
end



# destroy a random route and some part of other routes and re-insert them into the solution, a set of remove and insert heuristics are randomly paired
function remove_route(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)                                                         
    
    update_sol(soluT, solution) #copy solution as temporary solution soluT  
    n_cus  = darp.n_cus ; n_route = soluT.n_route
    max_degree_destruction = 0.275 # based on the tuned value from Roman Lutz (2014) Adaptive Large Neighborhood Search: A heuristic for the Rich Pickup
                                   # and Delivery Problem with Time Windows. Bachelor thesis at Ulm University.
    unserved_users =  soluT.unserved_users
    route_i = rand(1:n_route)     
    route   = get_route_no_depot(soluT, route_i)
    users   = route[findall(x->x<n_cus+2, route)] 
    soluT.vec_chg_event[route_i,:] .*= 0
    # put the users of the removed route into the pool
    unserved_users =  soluT.unserved_users
    union!(unserved_users, users) 
    candidate_users_remove = setdiff!(collect(2:n_cus+1), unserved_users) # pay attention here, user is from 2:n_cus+1
    length(candidate_users_remove) == 0 && return 
    tmp = floor(Int, length(candidate_users_remove) * max_degree_destruction) 
    degree_destruct = min(60, tmp)
    degree_destruct == 0 && (degree_destruct += 1)
    n_remove = rand(collect(1:degree_destruct)) # num of requests to remove from the other routes
 
    soluT.RI[end,4:5] =  soluT.RI[route_i,4:5] # move this unused veh (capacity and veh_id) to the end of the list
    soluT.vec_chg_event[end,:] .*= 0
    soluT.n_route -= 1    
    for i in route_i+1:n_route
        soluT.RI[i-1,:] = soluT.RI[i,:] 
        soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
    end 

    remove_op= [rand_remove, worst_remove, dist_remove, tw_remove, shaw_remove] # remove operators
    vec_remove = collect(1:length(remove_op))
    selected_rm_op = remove_op[rand(vec_remove)]
    RMOP(selected_rm_op, soluT, darp, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
    insert_unserved(soluT, recorder_routes, recorder_lscost_routes, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    # RMOP(selected_rm_op, soluT, darp, lgraph, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
    # insert_unserved(solution, soluT,  recorder_routes, recorder_lscost_routes, fast_chg, global_parameter, instance, darp, lgraph, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, T, type_veh, cap_type_veh) 
    update_sol(solution, soluT) 
end

# destroy a random route and some part of other routes and re-insert them into the solution, a set of remove and insert heuristics are randomly paired
# function remove_route(solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    
#     # @show("remove_route call !")
#     if solution.n_route < 2 
#         return 
#     end 
#     update_sol(soluT, solution) #copy solution as temporary solution soluT  
#     n_cus  = darp.n_cus ; n_route = soluT.n_route
#     max_degree_destruction = 0.275 # based on the tuned value from Roman Lutz (2014) Adaptive Large Neighborhood Search: A heuristic for the Rich Pickup
#                                    # and Delivery Problem with Time Windows. Bachelor thesis at Ulm University.                   
#     route_i = rand(1:n_route)     
#     route   = get_route_no_depot(soluT, route_i)
#     users   = route[findall(x->x<n_cus+2, route)] 
#     soluT.vec_chg_event[route_i,:] .*= 0
#     # put the users of the removed route into the pool
#     unserved_users =  soluT.unserved_users
#     union!(unserved_users, users) 

#     # compute the num of users to be removed from the other routes
#     candidate_users_remove = setdiff!(collect(2:n_cus+1), unserved_users) # users to be removed from the other routes,  user is from 2:n_cus+1
#     tmp = floor(Int, length(candidate_users_remove) * max_degree_destruction) 
#     degree_destruct = min(60, tmp)
#     degree_destruct == 0 && (degree_destruct += 1)
#     n_remove = rand(collect(1:degree_destruct)) # num of requests to remove from the other routes

#     soluT.RI[end,4:5] =  soluT.RI[route_i,4:5] # move this unused veh (capacity and veh_id) to the end of the list
#     soluT.vec_chg_event[end,:] .*= 0
#     soluT.n_route -= 1    
#     for i in route_i+1:n_route
#         soluT.RI[i-1,:] = soluT.RI[i,:] 
#         soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
#     end 
 
#     remove_op= [rand_remove, worst_remove, dist_remove, tw_remove, shaw_remove] # remove operators
#     vec_remove = collect(1:length(remove_op))
#     # select remove operator 
#     selected_rm_op = remove_op[rand(vec_remove)]
#     RMOP(selected_rm_op, soluT, darp, fast_chg, instance, n_remove, unserved_users, candidate_users_remove, relatedness, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
#     # @show("n_route after remove: ", selected_rm_op, soluT.n_route, n_remove, unserved_users)
#     # insert all unserved users using greedy or regret insertion
#     insert_unserved(soluT, recorder_routes, recorder_lscost_routes, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
#     update_sol(solution, soluT) 
# end

# # relocate rand a user and insert it to a rand route when first feasible  
# function relocate_rand(solution::Solution, soluT::Solution, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    
#     n_route = solution.n_route
#     if n_route < 2
#         return false
#     end
#     n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
#     qi =  darp.qi
#     dist_all, set_physical  = darp.dist_all, darp.set_physical 
#     K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy
#     penalty  = global_parameter.penalty   
#     chg_succ = false
#     update_sol(soluT,solution) #copy solution as temporary solution soluT
#     route_i=0; routes=Int32[]; users=Int32[]  

#     route_i = rand(1:n_route) 
#     routes = setdiff!(collect(1:n_route), route_i)

#     route = get_route(soluT, route_i, start_depot, end_depot)
#     users = route[findall(x->1<x<n_cus+2, route)]
   
#     shuffle!(users)
#     for user in users
#         shuffle!(routes)    
#         done =false 
#         for _r in routes
#             if  rand_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, start_depot, end_depot, Li) 
#                 done=true
#                 break
#             end
#         end        
#         if done             
#             if length(route)==4  # remove this route 
#                 soluT.RI[end,4:5] =  soluT.RI[route_i,4:5] # move this unused veh (capacity and veh_id) to the end of the list
#                 soluT.n_route -= 1
#                 for i in route_i+1:n_route
#                     soluT.RI[i-1,:] = soluT.RI[i,:] 
#                     soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
#                 end        
#                 update_sol(solution, soluT)              
#             else #update this route           
#                 filter!(x->x∉[user, user + n_cus], route)
#                 # soluT.TT_walk_routes[route_i] -= walk_tt_bus_stop[user-1]
#                 route_auxiliary_update(soluT, darp, route , route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)  
#                 # if eight_step(route, route_i, ei, li, s_t, dist_orig, Q, Li, n_cus, TH)
#                 if soluT.vec_chg_event[route_i, 1] > 0 
#                     _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
#                     if  chg_succ 
#                         update_sol(solution, soluT)   
#                     else
#                         repair_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start) && update_sol(solution, soluT)
#                     end
#                 end                     
#                 return true   
#             end 
#         end 
#     end
#     return false  # means that this operator has been done and count for performance score update    
# end 



# # update temporal solution
# function update_tmp_sol(soluT::Solution, r_i, vec_routes, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
    
#     n_route = soluT.n_route
#     if length(vec_routes[r_i])==2 # remove this route
#         soluT.vec_chg_event[r_i,:] .*= 0
#         soluT.RI[end,4:5] =  soluT.RI[r_i, 4:5] # move this unused veh (capacity and veh_id) to the end of the list
#         soluT.n_route -= 1     
#         for i in r_i+1:n_route
#             soluT.RI[i-1,:] = soluT.RI[i,:] 
#             soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
#         end
#     else           
#         route_auxiliary_update(soluT, darp, vec_routes[r_i], r_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)  
#     end
    
# end

#remove rand veh 
# or when co2 emission is grater than the threshold
# function remove_route(solution::Solution, soluT::Solution, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
#     # @show("remove_route call ")
#     n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
#     dist_all, set_physical  = darp.dist_all, darp.set_physical 
#     K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy   
#     flag_init_sol = false
#     n_route = solution.n_route
#     # @show(n_route) 
#     if n_route >1        
#         update_sol(soluT, solution) #copy solution as temporary solution soluT    
#         route_i = rand(1:n_route)     
#         route = get_route(soluT, route_i, start_depot, end_depot)
#         users = route[findall(x->1<x<n_cus+2, route)] 
#         route_copy= copy(route)
#         routes = setdiff!(collect(1:n_route), route_i)
#         shuffle!(users)
#         for user in users
#             shuffle!(routes)    
#             for _r in routes
#                 if rand() > 0.5  
#                     if greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, flag_init_sol) 
#                         filter!(x->x∉[user, user + n_cus], route_copy)
#                         break
#                     end
#                 else
#                     if rand_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, start_depot, end_depot, Li) 
#                         filter!(x->x∉[user, user + n_cus], route_copy)
#                         break
#                     end
#                 end
#             end
#         end
#         if length(route_copy)==2 # remove this route
#             soluT.n_route -= 1
#             for i in route_i+1:n_route
#                 soluT.RI[i-1,:] = soluT.RI[i,:] 
#                 soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
#             end
#             update_sol(solution, soluT)               
#             return true 
#         else #update this route           
#             route_auxiliary_update(soluT, darp, route_copy , route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back) 
#             # if  eight_step(route, route_i, ei, li, s_t, dist_orig, Q, Li, n_cus, TH)
#                 if soluT.vec_chg_event[route_i, 1] > 0 
#                     _ , chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, route_i, route_copy, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
#                     if  chg_succ
#                         update_sol(solution, soluT)                               
#                     else
#                         repair_charging(soluT, darp, route_i, route_copy, dist_all, fast_chg, parameter_energy, set_physical, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start) && update_sol(solution, soluT)
#                     end
#                 end
#             # else
#             #     repair_rand(soluT, soluT, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, Li, layer_nodes, lyrs_compatible)   
#             #     update_sol(solution, soluT)  
#             # end   
#             return true 
#         end
#     end 
# end


# # not used
# function create_route(solution::Solution, soluT::Solution, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    
#     n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
#     qi=  darp.qi
#     dist_all, set_physical  = darp.dist_all, darp.set_physical 
#     K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy  

#     n_route = solution.n_route 

#     if n_route < K && length(solution.unserved_users) > 0
#     # if  length(solution.unserved_users) > 0
#         route_i = n_route+1 
#         veh_used = solution.RI[1:n_route, 5]
#         vehs = setdiff(collect(1:K),veh_used)
#         veh_id = rand(vehs)
#         solution.RI[route_i, 4] = parameter_energy.cap_passenger[veh_id]
#         solution.RI[route_i, 5] = veh_id         
#         solution.n_route += 1 
#         user = rand(collect(solution.unserved_users)) # user denotes a bus stop
#         solution.RI[route_i,3]= 0 # make it clean
#         route = [start_depot, user, user+n_cus, end_depot]
#         insert_consecutive(solution, route_i, start_depot, user, end_depot, n_cus, qi) 
#         update_auxiliary_var(solution, darp, route_i, [start_depot, user, user+n_cus, user, user+n_cus, end_depot], e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back) # need to revise for flexbus problem
#         # @show(route_i, user, get_route(solution, route_i, start_depot, end_depot), "unserved_users (BE): ", solution.unserved_users)
#         setdiff!(solution.unserved_users, user)
#         solution.penalty_unserved -= qi[user]* penalty
#         if fast_chg.n_fast_chg_installed >0
#             update_charging(solution, darp, route_i, route , dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
#         end
#         return true
#     end
#     return false 
       
# end


# randomly run one of the three relocation LS operators
# function relocate_ensemble(solution::Solution, soluT::Solution, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)

#     ls_relocation= [relocate, relocate_worst, relocate_rand, insert_unserved]
#     ls_idx= rand(collect(1:length(ls_relocation)))
#     ls = ls_relocation[ls_idx]    
#     return LS(ls, solution, soluT, demand_data, fast_chg, global_parameter, instance, darp, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, T, Li)

# end


# function relocate(solution::Solution, soluT::Solution, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    
#     n_route = solution.n_route
#     if n_route < 2 
#         return false
#     end

#     n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
#     qi =  darp.qi
#     dist_all, set_physical  = darp.dist_all, darp.set_physical 
#     K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy  
#     # penalty, max_n_removed = global_parameter.penalty, global_parameter.max_n_request_removed
#     # n_request_removed = rand(1:max_n_removed) # added on 30.8.2023
#     penalty = global_parameter.penalty 

#     flag_init_sol = false
#     chg_succ = false
#     update_sol(soluT,solution) #copy solution as temporary solution soluT
#     route_i=0; routes=Int32[]; users=Int32[]    
#     # if length(soluT.unserved_users) > 0
#     #     route_i = rand(1:n_route+1) # n_route +1 denotes unserved_users
#     #     routes = setdiff!(collect(1:n_route+1),route_i)
#     # else
#         route_i = rand(1:n_route) 
#         routes = setdiff!(collect(1:n_route),route_i)
#     # end
    
#     # if route_i > n_route
#     #     users = collect(soluT.unserved_users)
#     # else
#         route = get_route(soluT, route_i, start_depot, end_depot)
#         users = route[findall(x->1<x<n_cus+2, route)]
#     # end
    
#     shuffle!(users)
#     # n_request_removed = min(n_request_removed, length(users))
#     # removed_users = users[1:n_request_removed] 
#     for user in users
#         shuffle!(routes)    
#         done =false 
#         for _r in routes
#             if greedy_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, flag_init_sol) 
#                 done=true
#                 break
#             end
#         end
#         if done 
#             # if route_i > n_route
#             #     setdiff!(soluT.unserved_users, user)
#             #     soluT.penalty_unserved -= qi[user] * penalty     
#             #     update_sol(solution,soluT)  
#             #     return true
#             # else
#                 if length(route)==4 # remove this route
#                     # @show("relocate: route_removed")
#                     soluT.n_route -= 1
#                     for i in route_i+1:n_route
#                         soluT.RI[i-1,:] = soluT.RI[i,:]
#                         # soluT.TT_walk_routes[i-1] = soluT.TT_walk_routes[i]
#                         soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
#                     end      
#                     update_sol(solution, soluT)               
#                     return true 
#                 else #update this route           
#                     filter!(x->x∉[user, user + n_cus], route)
#                     # soluT.TT_walk_routes[route_i] -= walk_tt_bus_stop[user-1]
#                     route_auxiliary_update(soluT, darp, route , route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back) 
#                     # if  eight_step(route, route_i, ei, li, s_t, dist_orig, Q, Li, n_cus, TH)
#                         if soluT.vec_chg_event[route_i, 1] > 0 
#                             _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
#                             if  chg_succ
#                                 update_sol(solution, soluT)                               
#                             else
#                                 repair_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start) && update_sol(solution, soluT)
#                             end
#                         end 
#                     return true 
#                 end
#             # end
#         end 
#     end
#     return false  # means that this operator has been done and count for performance score update     
# end

# relocate a segment of a route (zero split, Parragh et al.,2010) and insert these users one by one to another rand route
# time consuiming and not very effective, not used finally
# function relocate_seg(solution::Solution, soluT::Solution, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back,  Q, T, Li, layer_nodes, lyrs_compatible, type_veh, cap_type_veh)
    
#     start_depot=1
#     end_depot = 2*n_cus+2 # end depot node
#     update_sol(soluT,solution) #copy solution as temporary solution soluT
#     n_route = soluT.n_route
#     ri = rand(1:n_route)
#     route_i  = get_route(soluT, ri, start_depot, end_depot)
#     reset_cap_r(soluT, ri, cap_r)
#     cap_r_ri = [cap_r[vi] for vi in route_i] 
#     idx_ri   = collect(findall(x->x==0, cap_r_ri)) # find location of zero load    
#     routes   = setdiff!(collect(1:n_route),ri)
 
#     si =Int32[]
#     idx = rand(2:length(idx_ri)-1)
#     pre_idx = idx_ri[idx-1]
    
#     si = route_i[pre_idx+1:idx_ri[idx]]
#     route_i_new = zeros(Int32,2*n_cus)
#     users = si[findall(x->1<x<n_cus+2, si)]
#     n_route_i_new = length(route_i)-length(users)*2
#     length(users) > 1 && shuffle!(users)
#     seg_temp = Int32[]
#     # @show( ri, route_i, si, n_route_i_new, users, routes)
#     for user in users
#          shuffle!(routes)
#         if ! greedy_insert(soluT, user, routes[1], e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, ei, li, start_depot, end_depot, Li, layer_nodes, lyrs_compatible) 
#              push!(seg_temp, user, user+n_cus)
#         else
#             # soluT.RI[route_i,5] -= qi[user]
#             soluT.TT_walk_routes[route_i] -= walk_tt_bus_stop[user-1]
#         end
#     end
#     # @show(seg_temp)
#     if length(seg_temp) == 0
#         route_i_new[1:pre_idx], route_i_new[pre_idx+1:n_route_i_new] = route_i[1:pre_idx], route_i[idx_ri[idx]+1:end]
#         route_i_new = route_i_new[1:n_route_i_new]
#     else
#         n_seg = length(seg_temp)
#         n_route_i_new += n_seg
#         # @show( route_i[1:pre_idx],seg_temp, route_i[idx_ri[idx]+1:end] )
#         route_i_new[1:pre_idx], route_i_new[pre_idx+1 : pre_idx + n_seg] = route_i[1:pre_idx], seg_temp
#         route_i_new[pre_idx+n_seg+1:n_route_i_new] = route_i[idx_ri[idx]+1:end]
#         route_i_new = route_i_new[1:n_route_i_new]
#     end
#     ################
#     # if use this LS, check the following @show()
#     ################
#     @show(route_i_new, get_route(soluT, ri,start_depot,end_depot))  
   
#     if length(route_i_new)==2  # remove this route
#         soluT.n_route -= 1
#         !homogenous_fleet && push, soluT.RI[route_i,4])
#         for i in ri+1:n_route
#             soluT.RI[i-1,:] = soluT.RI[i,:]
#             soluT.TT_walk_routes[i-1] = soluT.TT_walk_routes[i]
#             soluT.vec_chg_event[i-1, :] = soluT.vec_chg_event[i, :]
#         end
#         update_sol(solution,soluT) 
#         return true 
#     else #update this route      
#         @show(route_i_new, get_route(soluT, ri,start_depot,end_depot))  
#         @error("sdsdsd")   
#         route_auxiliary_update(soluT, route_i_new , ri, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)  
#         if verify_eight_step(soluT, ri, Q, Li) 
#                update_charging( soluT, ri, route_i_new, dist_all, fast_chg, parameter_energy, set_physical, nodes, dt_r, dt_r_back, T_max)          
#         else
#               repair(soluT, soluT, ri, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back,  Q, Li, layer_nodes, lyrs_compatible)    
#         end   
#         update_sol(solution,soluT)  
#         return true      
#     end    
#     return false  # means that this operator has been done and count for performance score update     
# end


# function insert_unserved(solution::Solution, soluT::Solution, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    
#     if length(soluT.unserved_users) > 0

#         start_depot, end_depot = darp.start_depot, darp.end_depot
#         qi = darp.qi
#         penalty= global_parameter.penalty
#         update_sol(soluT, solution) #copy solution as temporary solution soluT
#         n_route = soluT.n_route
#         users = shuffle!(collect(soluT.unserved_users))
#         routes= shuffle!(collect(1:n_route))  
#         for user in users
#             done =false 
#             for _r in routes
#                 if  rand_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r, cap_r, l_r, maxcap_r,  dt_r, dt_r_back, Q, start_depot, end_depot, Li)
#                     done=true
#                     break
#                 end
#             end
#             if done 
#                 setdiff!(soluT.unserved_users, user)
#                 soluT.penalty_unserved -= qi[user] * penalty  
#                 update_sol(solution, soluT)  
#             end 
#         end
#     end
# end


# # relocate a worst user of a rand route and insert it to a rand route when first feasible 
# function relocate_worst(solution::Solution, soluT::Solution, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
    
#     n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
#     qi = darp.qi
#     dist_all, set_physical, dist  = darp.dist_all, darp.set_physical, darp.dist
#     K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy  
#     penalty = global_parameter.penalty
#     flag_init_sol = false
   
#     update_sol(soluT, solution) #copy solution as temporary solution soluT
#     route_i = 0; n_route = soluT.n_route; route = Int32[]; routes = Int32[]
#     # if length(soluT.unserved_users) > 0
#     #     route_i = rand(1:n_route+1) # n_route +1 denotes unserved_users
#     #     routes = setdiff!(collect(1:n_route+1), route_i)
#     # else
#         route_i = rand(1:n_route) 
#         routes = setdiff!(collect(1:n_route), route_i)
#     # end
#     # route_i < n_route+1 && (route = get_route(soluT, route_i, start_depot, end_depot))
#     route = get_route(soluT, route_i, start_depot, end_depot)
#     user = worst_user(soluT, qi, route_i, start_depot, end_depot, n_cus, dist)
#     # @show(user,route_i, soluT.n_route, routes)
#     if user > 0
#         shuffle!(routes)    
#         done =false 
#         for _r in routes            
#             if  rand_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r, cap_r, l_r, maxcap_r,  dt_r, dt_r_back, Q, start_depot, end_depot, Li)
#                 done=true
#                 break
#             end
#         end
#         if done 
#             # if route_i > n_route
#             #     setdiff!(soluT.unserved_users, user)
#             #     soluT.penalty_unserved -= qi[user] * penalty  
#             #     update_sol(solution, soluT)  
#             #     # @show("test 32: ") 
#             #     # verify_solution(solution, Q)
#             #     return true            
#             # else
#                 if length(route)==4  # remove this route
#                     soluT.n_route -= 1
#                     for i in route_i+1:n_route
#                         soluT.RI[i-1,:] = soluT.RI[i,:]
#                         soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
#                     end     
#                     update_sol(solution, soluT)            
#                     return true  
#                 else #update this route         
#                     filter!(x->x∉[user, user + n_cus], route) 
#                     route_auxiliary_update(soluT, darp, route , route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)  
#                     # if eight_step(route, route_i, ei, li, s_t, dist_orig, Q, Li, n_cus, TH)
#                         if soluT.vec_chg_event[route_i, 1] > 0 
#                             _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
#                             if  chg_succ
#                                 update_sol(solution, soluT)  
#                             else
#                                 repair_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start) && update_sol(solution, soluT)
#                             end
#                         end
#                     # else
#                     #     repair_rand(soluT, soluT, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, Li, layer_nodes, lyrs_compatible)   
#                     #     update_sol(solution, soluT)
#                     #     # @show("test 35: ")  ;    verify_solution(solution, Q)  
#                     # end   
#                     return true 
#                 end
#             # end
#         end 
#     end
#     return false  # means that this operator has been done and count for performance score update
# end 


# relocate a worst user on a rand route and insert it to a rand route in greedy mode
# not effective and tend to be trapped in local optimal and not complementary for the other LS operators
# not used, need to modify if used
# function relocate_worst_greedy(solution::Solution, soluT::Solution, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li, layer_nodes, lyrs_compatible, type_veh, cap_type_veh)
#     start_depot=1; end_depot = 2*n_cus+2 # end depot node
#     update_sol(soluT, solution) #copy solution as temporary solution soluT
#     route_i = 0; n_route = soluT.n_route; route = Int32[]; routes = Int32[]
#     if length(soluT.unserved_users) > 0
#         route_i = rand(1:n_route+1) # n_route +1 denotes unserved_users
#         routes = setdiff!(collect(1:n_route+1), route_i)
#     else
#         route_i = rand(1:n_route) 
#         routes = setdiff!(collect(1:n_route), route_i)
#     end
#     route_i < n_route+1 && (route = get_route(soluT, route_i, start_depot, end_depot))
#     user = worst_user(soluT, route_i, start_depot, end_depot)
#     # @show(user,route_i, soluT.n_route, routes)
#     if user > 0
#         shuffle!(routes)    
#         done =false 
#         for _r in routes
#             if greedy_insert(soluT, user, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q,  ei, li, start_depot, end_depot, Li, layer_nodes, lyrs_compatible)
#                done=true
#                 # @show(user, _r,solution.unserved_users)
#                 break
#             end 
#         end
#         if done 
#             if route_i > n_route
#                 setdiff!(soluT.unserved_users, user)
#                 soluT.penalty_unserved -= qi[user] * penalty  
#                 update_sol(solution, soluT)   
#                 return true            
#             else
#                 if length(route)==4  # remove this route
#                     soluT.n_route -= 1
#                     !homogenous_fleet && push, soluT.RI[route_i,4])
#                     for i in route_i+1:n_route
#                         soluT.RI[i-1,:] = soluT.RI[i,:]
#                         soluT.TT_walk_routes[i-1] = soluT.TT_walk_routes[i]
#                         soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
#                     end
#                     # no need to update chg state occupancy 
#                     # verify_solution(soluT) 
#                     update_sol(solution,soluT) 
#                     return true 
#                 else #update this route           
#                     filter!(x->x∉[user, user + n_cus], route)
#                     # soluT.RI[route_i,5] -= qi[user]
#                     soluT.TT_walk_routes[route_i] -= walk_tt_bus_stop[user-1]
#                     route_auxiliary_update(soluT, route , route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)  
#                     if verify_eight_step(soluT, route_i, Q, Li) 
#                          update_charging(soluT, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, nodes, dt_r, dt_r_back, T_max)
#                     else
#                          repair_rand(soluT, soluT, route_i,e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, Li, layer_nodes, lyrs_compatible)   
#                     end    
#                     update_sol(solution,soluT)  
#                     return true                           
#                 end
#             end
#         end 
#     end
#     return false  # means that this operator has been done and count for performance score update
    
# end 


# "two user requests are swapped, the pickup (delivery) vertex of the first request
# # may only be inserted at the same compatible layers position of pickup (delivery) vertex of the second user request in the second route, while
# # the pickup and delivery vertices of the second user request may be inserted at any position in the first route."
# function swap_two_users(solution::Solution, soluT::Solution, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li, layer_nodes, lyrs_compatible, type_veh, cap_type_veh)
    
#     start_depot=1; end_depot = 2*n_cus+2 # end depot node
#     update_sol(soluT, solution) #copy solution as temporary solution soluT
#     # insert from route_i to route j
#     init_route_info(e_r,cap_r,l_r,maxcap_r, li, start_depot, end_depot, dt_r, dt_r_back)
#     n_route = soluT.n_route
#     if n_route < 2
#         return false
#     end
#     route_1 ,route_2 = shuffle!(collect(1:n_route))[1:2]
#     if soluT.RI[route_1,3] < soluT.RI[route_2,3]
#         r_i = route_1; r_j=route_2
#     else
#         r_i = route_2; r_j=route_1
#     end 
#     seg_i1=Int32[];seg_i2=Int32[];seg_j1=Int32[];seg_j2=Int32[];copy_seg_j2 =Int32[];copy_seg_i2 =Int32[]
#     route_j_new = zeros(Int32,2*n_cus); route_i_new = zeros(Int32,2*n_cus); new_route = zeros(Int32,4)
#     nodes_route = Int32[]; route_j_tmp = Int32[]
#     idx3=0 ;  dist_saving = 0; route =Int32[]
#     if soluT.RI[r_i,3]>3
#         route_i = get_route(soluT, r_i, start_depot, end_depot)
#         route_j = get_route(soluT, r_j, start_depot, end_depot)
#         # insert from route_i to route j
#         idx_users_ri   = collect(findall(x->x<n_cus+2, route_i))
#         idx_users_rj   = collect(findall(x->x<n_cus+2, route_j))
#         # @show(route_i,idx_users_ri,route_j,idx_users_rj)
#         for idx in 2:length(idx_users_ri)
#             idx1 = idx_users_ri[idx]
#             idx2 = idx_users_rj[idx]
#             vi   = route_i[idx1]
#             vj   = route_j[idx2]
#             # if  is_compatible(vi, vj, layer_nodes, lyrs_compatible, end_depot)
#             pre_vi = route_i[idx1-1]
#             suc_vi = route_i[idx1+1]
#             pre_vj = route_j[idx2-1]
#             suc_vj = route_j[idx2+1]
#             v_ni = vi+n_cus
#             v_nj = vj+n_cus
#             seg_i1=route_i[1:idx1-1];seg_i2=route_i[idx1+1:end]
#             seg_j1=route_j[1:idx2-1];seg_j2=route_j[idx2+1:end]  
#             pre_vj=soluT.pre[vj]  
#             pre_v_ni = soluT.pre[v_ni]
#             suc_v_ni = soluT.succ[v_ni] 
#             pre_v_nj = soluT.pre[v_nj]
#             suc_v_nj = soluT.succ[v_nj]         
#             dist_saving = dist[pre_vi,vi]+dist[vi,suc_vi]+dist[pre_vj,vj]+dist[vj,suc_vj]
#             + dist[pre_v_ni,v_ni]+dist[v_ni,suc_v_ni]+dist[pre_v_nj,v_nj]+dist[v_nj,suc_v_nj]
#             - (dist[pre_vj,vi]+dist[vi,suc_vj]+dist[pre_v_nj,v_ni]+dist[v_ni,suc_v_nj] 
#             + dist[pre_vi,vj]+dist[vj,suc_vi]+dist[pre_v_ni,v_nj]+dist[v_nj,suc_v_ni])
#             if dist_saving > T
#                 # check tw feasibility if inserting vi at the position of vj
#                 eh_vi = max(ei[vi], e_r[pre_vj]+s_t[pre_vj]+dist[pre_vj,vi])
#                 if eh_vi <= li[vi]  
#                     # @show(route_j,seg_j1,vj,seg_j2)
#                     idx3 = findfirst(x->x==v_nj, seg_j2)
#                     feasible = true
#                     pre_v = vi 
#                     eh_pre_v = eh_vi
#                     for i in 1:idx3-1
#                         next = seg_j2[i]
#                         eh_pre_v = max(ei[next], eh_pre_v + s_t[pre_v]+dist[pre_v,next])
#                         if eh_pre_v > li[next]  
#                             feasible = false
#                             break
#                         end
#                         pre_v = next
#                     end
#                     if feasible 
#                         eh_v_ni    = max(ei[v_ni], eh_pre_v + s_t[pre_v]+dist[pre_v,v_ni])
#                         eh_suc_v_nj= max(ei[suc_v_nj], eh_v_ni + s_t[v_ni]+dist[v_ni,suc_v_nj]) 
#                         if eh_v_ni<=li[v_ni] && eh_suc_v_nj <= li[suc_v_nj]
#                             if eh_v_ni - li[vi] - s_t[vi] <= Li[vi] 
#                                 copy_seg_j2 =copy(seg_j2) 
#                                 copy_seg_j2[idx3] = v_ni
#                                 n1=length(seg_j1); n2=length(copy_seg_j2)
#                                 route_j_new[1:n1],route_j_new[n1+1] = seg_j1,vi
#                                 route_j_new[n1+2:n1+n2+1] = copy_seg_j2
#                                 nodes_route = nodes[route_j_new[2:n1+n2]]
#                                 route_j_tmp = route_j_new[1:n1+n2+1]
#                                 if no_double_tour(nodes_route, D′) && eight_step(route_j_tmp, r_j, ei, li, s_t, dist, Q, Li, n_cus, TH)   
#                                     route_auxiliary_update(soluT, route_j_tmp, r_j, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
#                                     chg_succ_j, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, r_j, route_j_new[1:n1+n2+1], dist_all, fast_chg, parameter_energy, set_physical, nodes, dt_r, dt_r_back, T_max)
#                                     if  ! chg_succ_j
#                                         chg_succ_j = repair_charging(soluT, r_j, route_j_tmp, dist_all, fast_chg, parameter_energy, set_physical, nodes, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start)                                  
#                                     end
#                                     if chg_succ_j
#                                         copy_seg_i2 = copy(seg_i2)
#                                         deleteat!(copy_seg_i2,  findfirst(x->x==v_ni, copy_seg_i2));
#                                         n1=length(seg_i1) ; n2 =length(copy_seg_i2)
#                                         route_i_new[1:n1], route_i_new[n1+1:n1+n2]=  seg_i1,copy_seg_i2
#                                         route_i_tmp =  route_i_new[1:n1+n2]
#                                         route_auxiliary_update(soluT, route_i_tmp, r_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)
#                                         if greedy_insert(soluT, vj, r_i, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, ei, li, start_depot, end_depot, Li, layer_nodes, lyrs_compatible)
#                                             # @show("test 2: ")
#                                             # @show(vj, r_i, get_route(soluT, r_i, start_depot, end_depot), soluT.TT_walk_routes)
#                                             soluT.TT_walk_routes[r_j] += walk_tt_bus_stop[vi-1]
#                                             soluT.TT_walk_routes[r_i] -= walk_tt_bus_stop[vi-1]      
#                                             soluT.TT_walk_routes[r_j] -= walk_tt_bus_stop[vj-1] 
#                                             update_sol(solution, soluT)
#                                             return true
#                                         else  
#                                             routes = setdiff!(collect(1:n_route), [r_i,r_j])                              
#                                             if length(routes)>0
#                                                 shuffle!(routes) 
#                                                 for _r in routes
#                                                     if greedy_insert(soluT, vj, _r, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, ei, li, start_depot, end_depot, Li, layer_nodes, lyrs_compatible)                        
#                                                         if verify_eight_step(soluT, r_i, Q, Li) 
#                                                             route = get_route(soluT, r_i, start_depot, end_depot)
#                                                             chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, r_i, route, dist_all, fast_chg, parameter_energy, set_physical, nodes, dt_r, dt_r_back, T_max)
#                                                             if  chg_succ
#                                                                 soluT.TT_walk_routes[r_j] += walk_tt_bus_stop[vi-1]
#                                                                 soluT.TT_walk_routes[r_i] -= walk_tt_bus_stop[vi-1]      
#                                                                 soluT.TT_walk_routes[r_j] -= walk_tt_bus_stop[vj-1] 
#                                                                 update_sol(solution, soluT)  
#                                                                 # @show("test 3 : ", r_i, r_j)
#                                                                 # @show([get_route(soluT, i, start_depot,end_depot) for i in 1:soluT.n_route]) 
#                                                                 # @show(soluT.vec_chg_event[1:K, :])
#                                                                 # verify_solution(soluT, Q) 
#                                                             else
#                                                                 if repair_charging(soluT, r_i, route, dist_all, fast_chg, parameter_energy, set_physical, nodes, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start) 
#                                                                     soluT.TT_walk_routes[r_j] += walk_tt_bus_stop[vi-1]
#                                                                     soluT.TT_walk_routes[r_i] -= walk_tt_bus_stop[vi-1]      
#                                                                     soluT.TT_walk_routes[r_j] -= walk_tt_bus_stop[vj-1] 
#                                                                     update_sol(solution, soluT)    
#                                                                 end  
#                                                                 # @show("test 4: ")
#                                                                 # verify_solution(soluT, Q)                  
#                                                             end
#                                                         else
#                                                             if repair(soluT, soluT, r_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, Li, layer_nodes, lyrs_compatible) 
#                                                                 soluT.TT_walk_routes[r_j] += walk_tt_bus_stop[vi-1]
#                                                                 soluT.TT_walk_routes[r_i] -= walk_tt_bus_stop[vi-1]      
#                                                                 soluT.TT_walk_routes[r_j] -= walk_tt_bus_stop[vj-1] 
#                                                                 update_sol(solution, soluT)  
#                                                             end
#                                                         end
#                                                         return true
#                                                     end
#                                                 end                                                 
#                                                 return true                    
#                                             else        
#                                                 return true
#                                             end
#                                         end
#                                     end
#                                 end
#                             end                        
#                         end
#                     end
#                 end
#             end
#             # end
#         end 
#     end
#     return false
# end


#repair an infeasible route in a first feasible way
# function repair_rand(solution::Solution,soluT::Solution, demand_data, fast_chg, global_parameter, instance, darp, route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, Li, layer_nodes, lyrs_compatible)
    
#     n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
#     dist_all, set_physical  = darp.dist_all, darp.set_physical 
#     K, T_max, parameter_energy = instance.K, instance.T_max, instance.parameter_energy  

#     chg_succ = false
#     update_sol(soluT,solution) #copy solution as temporary solution soluT
#     n_route = soluT.n_route        
#     route = get_route(soluT, route_i, start_depot, end_depot) 
#     users= route[findall(x->1<x<n_cus+2, route)]
#     routes = setdiff!(collect(1:n_route),route_i)
#     user = shuffle!(users)
#     for user in users
#         shuffle!(routes)    
#         done = false 
#         for _r in routes
#             if rand_insert(soluT, fast_chg, global_parameter, instance, darp, user, _r, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, start_depot, end_depot, Li ) 
#                 done = true
#                 break
#             end
#         end
#         if done 
#             if length(route)==4  # remove this route
#                 soluT.n_route -= 1      
#                 for i in route_i+1:n_route
#                     soluT.RI[i-1,:] = soluT.RI[i,:]
#                     # soluT.TT_walk_routes[i-1]   = soluT.TT_walk_routes[i]                    
#                     soluT.vec_chg_event[i-1,:]  = soluT.vec_chg_event[i,:]
#                 end 
#                 update_sol(solution, soluT)                 
#                 return true 
#             else #update this route           
#                 filter!(x->x∉[user, user + n_cus], route)
#                 # soluT.TT_walk_routes[route_i] -= walk_tt_bus_stop[user-1]
#                 route_auxiliary_update(soluT, darp, route , route_i, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back)  
#                 # if eight_step(route, route_i, ei, li, s_t, dist_orig, Q, Li, n_cus, TH)
#                     if soluT.vec_chg_event[route_i, 1] > 0 
#                         _, chg_succ, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back, T_max, K)
#                         if  chg_succ
#                             update_sol(solution, soluT)  
#                         else
#                             repair_charging(soluT, darp, route_i, route, dist_all, fast_chg, parameter_energy, set_physical, dt_r, pos_insert_chg, forward_time_slack, vec_t_start) && update_sol(solution, soluT)
#                         end
#                     end
#                     return true
#                 # end         
#             end
#         end 
#     end     
# end