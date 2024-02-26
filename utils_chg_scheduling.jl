


function status_chg_occupancy(solution::Solution, set_unchanged_route, parameter_energy, fast_chg, T_max)

    discrte_t_interval = 1/6 # 10 second
    n_route = solution.n_route  
    n_t_max = floor(Int32, (T_max+20)/discrte_t_interval) # 840 discrete charging occupancy state over T_max+20 minutes
    n_charger = fast_chg.n_fast_chg
    occ_state_chg = zeros(Int32, n_charger, n_t_max) 
    vec_chg_event = solution.vec_chg_event[1:n_route,:] 
    for r in set_unchanged_route 
        for i in 1:floor(Int32, vec_chg_event[r, 2])
            idx_chg, t1, t2 = vec_chg_event[r, 2+(i-1)*4+2], vec_chg_event[r, 2+(i-1)*4+3],vec_chg_event[r, 2+(i-1)*4+4]
            idx_chg = floor(Int32,idx_chg)
            h_t1, h_t2 = floor(Int32, t1/discrte_t_interval)+1, floor(Int32, t2/discrte_t_interval)+1
            occ_state_chg[idx_chg, h_t1: h_t2] .+=  1 
        end
    end
    return occ_state_chg
end


function update_chg_state(idx_chg, t1, t2, occ_state_chg, discrte_t_interval)    
    
    idx_chg = floor(Int32,idx_chg)
    h_t1, h_t2 = floor(Int32, t1/discrte_t_interval)+1, floor(Int32, t2/discrte_t_interval)+1
    occ_state_chg[idx_chg, h_t1: h_t2] .=  1
    
end


# schedule charging events for the new routes in a greedy way without conflicts with other vehicles
# v1 in vec_chg_event is node id in V_DF
function update_chg_new_sol(solution::Solution, lgraph, parameter_energy, routes_new, changed_routes, dist_all, set_physical, fast_chg, V_DF, T_max, nodes)                 
    
    discrte_t_interval = 1/6
    n_route = solution.n_route
    vec_chg_event = copy(solution.vec_chg_event[1:n_route,:])
    changed_routes= sort(collect(changed_routes))
    set_unchanged_route = setdiff(collect(1:n_route),changed_routes) 
    occ_state_chg= status_chg_occupancy(solution, set_unchanged_route, parameter_energy, fast_chg, T_max) 
    #  # change node index for the unchanged routes
    # for r in set_unchanged_route
    #     n_chg_event = floor(Int32, vec_chg_event[r,2])
    #     for i in 1:n_chg_event
    #         idx = 2+(i-1)*4 + 1 
    #         v1= floor(Int32, vec_chg_event[r,idx])
    #         vec_chg_event[r,idx] = nodes[v1] 
    #     end
    # end
    for (idx, r) in enumerate(changed_routes)
        route= routes_new[r];veh_id = solution.RI[r,5] 
        success, info_chg = insert_charging_cus_re_assign(solution, lgraph, r, route, fast_chg, parameter_energy, occ_state_chg, V_DF)
        if success  
            vec_chg_event[r,:] = info_chg
            if info_chg[2]>0
                idx_chg, t1, t2 = info_chg[4], info_chg[5], info_chg[6]
                update_chg_state(idx_chg, t1, t2, occ_state_chg, discrte_t_interval)    
            end
        else
             return false, vec_chg_event
        end
    end   
    return true, vec_chg_event  
end
