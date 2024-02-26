
function preprocess(n_cus, li, dist_orig, s_t)
 
    li[1]=li[2*n_cus+2] = maximum([li[i]+s_t[i] + dist_orig[i, 2*n_cus+2] for i in 2:2*n_cus+1])
 
end


function parameter_settings(set_seed_number, percentage_full, step_size_rm_route)

    # parameters depending on user-defined problems
    step_size_remove_route = step_size_rm_route # default for num of requests <=50, for larger instances, can set it as 250
    step_size_create_route = 500 
    max_num_chg = 4 # of for c500 case, E_init=100%, reduce 90% CO2 emissions
    n_run_darp = 3   
    is_gasoline_fleet = false # initialization
    co2_reduc_target = 0  
    n_type_veh = 2 # updated in the read_data file gasoline and EB-1, EB-2, will be replaced based on the data input
    detour_factor =1.5
    duplcate_charger = 3 
    bigM = 1000; 
    penalty = 60#40 # 100 # 50 calbrated penalty for one unserved customer (fleetsize), 200 (old value)
    penalty_per_veh_used = 200#150 #calbrated penalty (fleetsize), 100 old value
    big_penalty_co2  = 100*1000
    degree_rand_chg_policy =1 # don't change it
    n_init_sol = 1000
    flag_output_to_file = true
    N_ITER = 100 *1000             # calibrated value of fleetsize
    max_stagnant_multiplier = 150#200  # calibrated value 200 
    T_red = 200; t_max = 2.1;  n_imp = 100 # calibrated value of our revised TRE paper, ref
    tw_width = 10 # 
    VoWT = 1 # not used
    parameter_sa = Parameter_SA(N_ITER, step_size_remove_route, max_stagnant_multiplier, T_red, t_max, n_imp,penalty,penalty_per_veh_used)
    n_scenario=1
    flag_initialization = true  
    unit_charger_cost = 10000 # dummy init
    vec_cs_site_cost  = Float32[]
    vec_max_elec_supply=Float32[]
    co2_gasoline_fleet=0
    rho = 0.9 # rho-quantile to determine fleet size
    co2_threshold = 1000*1000 # initialization
    co2_threshold_init = 1000*1000 # initialization
    max_passenger_request = 4 # max num of passengers per ride request
    parameter_chg_infra = Parameter_chg_infra(unit_charger_cost, vec_cs_site_cost, vec_max_elec_supply)
    solve_veh_exchange_t_limit = 2 * 60  # sec.
    # max_n_request_removed = 3
    #instantiate 
 
    return Global_parameter(is_gasoline_fleet, n_run_darp, set_seed_number, penalty, n_type_veh, 
                detour_factor, duplcate_charger, bigM,n_init_sol, flag_output_to_file, tw_width, VoWT, parameter_sa,
                n_scenario, flag_initialization, parameter_chg_infra, co2_reduc_target,
                co2_gasoline_fleet, co2_threshold, co2_threshold_init, rho, big_penalty_co2, max_passenger_request, 
                percentage_full, penalty_per_veh_used, max_num_chg, degree_rand_chg_policy, solve_veh_exchange_t_limit,
                step_size_create_route) 
end

# check whether the dropoff nodes are visited after the pickup nodes of customers
function check_pre_seg(segment::Vector{Int32}, n_cus)

    N_NODE = 2*n_cus+2
    visited_order =ones(Int32,N_NODE) * N_NODE
    idx=0
    for i in 1:length(segment)
        idx +=1
        visited_order[segment[i]] = idx
    end
    users = segment[findall(x->x<n_cus+2,  segment)]
    for user in users
        if visited_order[user+n_cus] - visited_order[user] < 0
            return false 
        end
    end
    return true
end

function check_route(route::Vector{Int32})
    idxs = findfirst(x->x<1, route)
    if !isnothing(idxs)
        return false
    else
        return true
    end
end

function check_precedence(route::Vector{Int32}, n_cus)

    visited_order =zeros(Int32,2*n_cus+2)
    idx=0; feasible = true; failed_user=0
    for i in 1:length(route)-1
        idx +=1
        visited_order[route[i]] = idx
        if route[i]<1 || route[i]>2*n_cus+2
            feasible = false
            break
        end
    end
    users = route[findall(x->1<x<n_cus+2,  route)]
    for user in users
        if visited_order[user+n_cus] - visited_order[user] < 0
            feasible = false; failed_user = user
             break
        end       
    end
    if ! feasible
        @show("check_precedence failed ! route : ", route, failed_user)
    end
    return feasible

end


function check_tw(route, t0, ei, li, s_t, dist)
 
    N = length(route)
    A =zeros(Float32, N)  #Arrival time at vertex i 
    B =zeros(Float32, N)  #Beginning of service at vertex i 
    W =zeros(Float32, N)  #Wait time before starting service at vertex i 
    D =zeros(Float32, N)  #Departure time at vertex i  
    A[1] = t0 ;  B[1] = max(t0, ei[route[1]])
    D[1] = B[1] + s_t[route[1]];  W[1] = B[1] - A[1]
    for idx in 2:N  
        A[idx] = D[idx-1] + dist[route[idx-1], route[idx]]
        B[idx] = max(A[idx], ei[route[idx]])
        W[idx] = B[idx] - A[idx]
        D[idx] = B[idx] + s_t[route[idx]]    
        if round(B[idx] - li[route[idx]], digits= DIGITS) > 0 
            return false 
        end
    end
    return true
    
end

function check_tw_ride_time(route, t0, ei, li, s_t, dist)
 
    N = length(route)
    A =zeros(Float32, N)  #Arrival time at vertex i 
    B =zeros(Float32, N)  #Beginning of service at vertex i 
    W =zeros(Float32, N)  #Wait time before starting service at vertex i 
    D =zeros(Float32, N)  #Departure time at vertex i  
    RT =zeros(Float32, N) 
    A[1] = t0 ;  B[1] = max(t0, ei[route[1]])
    D[1] = B[1] + s_t[route[1]];  W[1] = B[1] - A[1]
    for idx in 2:N-1
        A[idx] = D[idx-1] + dist[route[idx-1], route[idx]]
        B[idx] = max(A[idx], ei[route[idx]])
        W[idx] = B[idx] - A[idx]
        D[idx] = B[idx] + s_t[route[idx]]    
        if round(B[idx] - li[route[idx]], digits= DIGITS) > 0 
            return false 
        end
        if  route[idx] > n_cus+1
            idx2 = findfirst(x->x==route[idx]-n_cus, route)
            # RT[idx] = A[idx] - D[idx2]
            RT[idx] = B[idx] - D[idx2]
            if  round(RT[idx] - Li[route[idx]], digits= DIGITS) > 0
                return false
            end
        end
    end
    return true
    
end

# check tw constraint
function verify_tw_ride_time(fast_chg, darp, instance, route, r, ei_route, li_route, s_t_route, dist_orig, qi_route, Li_route,set_physical, dist_all)
 
    n_cus, TH, s_t = darp.n_cus, darp.TH, darp.s_t
    N_NODE = 2*n_cus+2
    Q = instance.Q_max
    N=length(route)
    t_start = ei_route[1]

    A, B, W, D, Load, RT = calculate_ads_chg(fast_chg, darp, route, t_start, ei_route, li_route, s_t_route, dist_orig, n_cus, qi_route,set_physical, dist_all)
    if ! check_tw_cap_chg(route, B, li_route, Load, r, Q) 
        @show("check_tw_cap(route, B ,Load) failed!! ")
        @show(route, B, li_route, Load, Q)
         return false, -1,-1,-1,-1,-1,-1
    end

    F0 = calc_Fi_chg(1, route, W, B, RT, Li_route, li_route)
    Wp = sum(W[2:N-1])
    t_start_new = ei_route[1] + min(F0, Wp)
    A, B, W, D, Load, RT = calculate_ads_chg(fast_chg, darp, route, t_start_new, ei_route, li_route, s_t_route, dist_orig, n_cus, qi_route, set_physical, dist_all)

    if check_ridetime_chg(route, RT, Li_route) && check_routeduration(route, B, TH) 
        return true, A, B, W, D, Load, RT 
    end
 
    for i in 2:N-1
        if  1 < route[i] <= n_cus+1 
            Fi = calc_Fi_chg(i, route, W, B, RT, Li_route, li_route)
            Wp = sum(W[i+1:N-1])              
            W[i] += min(Fi, Wp)
            B[i] = A[i] + W[i]
            D[i] = B[i] + s_t[route[i]]
            for j in i+1:N
                if  route[j-1] > N_NODE 
                    v1,v2 = set_physical[route[j-1]], set_physical[route[j]]
                    A[j] = D[j-1] + dist_all[v1, v2] 
                elseif route[j] > N_NODE
                    v1,v2 = set_physical[route[j-1]] , set_physical[route[j]] 
                    A[j] = D[j-1] + dist_all[v1, v2] 
                else
                    A[j] = D[j-1] + dist_orig[route[j-1], route[j]]
                end
                B[j] = max(A[j], ei_route[j])
                W[j] = B[j] - A[j]
                D[j] = B[j] + s_t_route[j]
                if  n_cus+1 < route[j] < N_NODE
                    idx2 = findfirst(x->x==route[j]-n_cus, route)
                    # RT[j] = A[j] - D[idx2]
                    RT[j] = B[j] - D[idx2]
                end
            end
        end
    end
    #check whether the drive time is respected for all dropoff vertices after i, if so the route is feasible 
    if check_routeduration(route, B, TH) && check_ridetime_chg(route, RT, Li_route)
        return true, A, B, W, D, Load, RT 
    else   
        ! check_routeduration(route, B, TH) && (@show("check_routeduration failed : ", check_routeduration(route, B, TH)))
        if ! check_ridetime_chg(route, RT, Li_route)
            @show("check_ridetime_chg failed : see detail at tw_ride_time_detail_r.xls", check_ridetime_chg(route, RT, Li_route))   
            tw_detail = DataFrame(node = route, A=A, B=B, W=W, D=D, RT=RT, Load=Load, ei = ei_route, li = li_route, st = s_t_route, Li_route=Li_route)
            XLSX.writetable("tw_ride_time_detail_r" * string(r) * ".xlsx", overwrite=true, tw_detail)   
            error("check_ridetime_chg failed !")
        end
        return false, -1,-1,-1,-1,-1,-1
    end

end


# function to check the final solution if the charging operations of the solution do not 
# violate the TW and ride time constraints of customers
function eight_step_chg(solution::Solution, r, chg_operation_r, fast_chg::Fast_chargers, darp)

    ei, li, qi, s_t, Li, dist_orig = darp.ei, darp.li, darp.qi, darp.s_t, darp.Li, darp.dist_orig
    dist_all, set_physical, n_cus =  darp.dist_all,  darp.set_physical,  darp.n_cus 

    start_depot =1; end_depot = 2* n_cus+2
    route_0 = get_route(solution, r, start_depot, end_depot)
    route = copy(route_0)
    ei_route  = ei[route]
    li_route  = li[route]
    s_t_route = s_t[route]
    qi_route  = qi[route]
    Li_route  = Li[route] 
    for i in 1:floor(Int32, chg_operation_r[2])
        idx0 = 2+(i-1)*4
        v1    = floor(Int32, chg_operation_r[idx0+1]) # node precede the charging station (charger) node 
        v_chg= fast_chg.v_chg[floor(Int32, chg_operation_r[idx0+2])]     
        t_chg_start = chg_operation_r[idx0+3]
        dur_chg     = chg_operation_r[idx0+4] - chg_operation_r[idx0+3]
        idx_insert = findfirst(x->x==v1, route) 
        insert!(route, idx_insert+1, v_chg); insert!(ei_route, idx_insert+1, t_chg_start ) 
        insert!(li_route, idx_insert+1, li[route[1]]); insert!(s_t_route, idx_insert+1,dur_chg )
        insert!(qi_route,idx_insert+1, 0 ); insert!(Li_route, idx_insert+1,li[end])
 
    end
    #check TW and ride time
    # @show(route, chg_operation_r)                            
    success, A, B, W, D, Load, RT = verify_tw_ride_time(fast_chg, darp, instance, route, r, ei_route, li_route, s_t_route, dist_orig, qi_route, Li_route, set_physical, dist_all)
    tw_detail = DataFrame(node = route, A=A, B=B, W=W, D=D, RT=RT, Load=Load, ei = ei_route, li = li_route, st = s_t_route, Li_route=Li_route )
    XLSX.writetable("tw_ride_time_detail_r" * string(r) * ".xlsx",overwrite=true, tw_detail)
    if !success 
        @show("verify_tw_ride_time failed! : ", r, route_0, route, chg_operation_r)
        error("verify_tw_ride_time failed! : ")
        return false
    else
        return true
    end
end

function verify_tw_chg_operation(solution::Solution, fast_chg::Fast_chargers, darp, instance)
    
    # n_cus, ei, li, qi, s_t, dist_orig, dist_all, set_physical = darp.n_cus, darp.ei, darp.li, darp.qi, darp.s_t, darp.dist_orig, darp.dist_all, darp.set_physical
    for r in 1:solution.n_route
        if  solution.vec_chg_event[r,2] > 0  
            chg_operation_r = solution.vec_chg_event[r,:]          
            eight_step_chg(solution, r, chg_operation_r, fast_chg, darp)
            # check_tw_RT_chg(solution, r, chg_operation_r, fast_chg, ei, li, s_t, dist_orig, dist_all)
        end
    end
    return true
end

function export_route_detail(path_result, solution::Solution, instance, darp, fast_chg)
    
    # @show(" export_route_detail call !", path_result)
    n_cus, start_depot, end_depot= darp.n_cus, darp.start_depot, darp.end_depot
    ei, li, s_t, dist, dist_orig, Li, qi, TH = darp.ei, darp.li, darp.s_t, darp.dist, darp.dist_orig, darp.Li, darp.qi, darp.TH
    Q= instance.Q_max
    for r in 1:solution.n_route
        if  solution.vec_chg_event[r,2] > 0  
            chg_operation_r = solution.vec_chg_event[r,:]   
            disp_route_chg(path_result, solution, darp, instance, r, chg_operation_r, fast_chg, ei, li, Li, s_t, dist_orig)
        else
            route = get_route(solution, r, start_depot, end_depot )
            disp_route(path_result, route, r, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)
        end
    end   
     
end

# for tw and ride time check when involving charging operations on a route
function calculate_arr_t(route::Vector{Int32}, t_start, ei_route, li_route, s_t_route, dist_orig, n_cus, qi_route, fast_chg)
    # 1 and 2n+2 are the depots, 2...n+1 (pickup), n+2...2n+1 dropoff
    n_node_darp = 2 * n_cus + 2
    N = length(route)
    A = zeros(Float32, N)  #Arrival time at vertex i 
    B = zeros(Float32, N)  #Beginning of service at vertex i 
    W = zeros(Float32, N)  #Wait time before starting service at vertex i 
    D = zeros(Float32, N)  #Departure time at vertex i 
    Load = zeros(Int32, N) #cummulative laod after leaving vertex i 
    RT = zeros(Float32, N) #ride time, only specified for pickup vertices else 0 
    A[1] = t_start
    B[1] = max(t_start, ei_route[1])
    D[1] = B[1] + s_t[route[1]]
    W[1] = B[1] - A[1]
    for idx in 2:N
        if route[idx-1] > n_node_darp
            v1, v2 = fast_chg.v_physical_chg[route[idx-1]-10000], set_physical[route[idx]]
            A[idx] = D[idx-1] + dist_all[v1, v2]
            # @show((route[idx-1], route[idx]),(v1, v2), dist_all[v1, v2])
        elseif route[idx] > n_node_darp
            v1, v2 = set_physical[route[idx-1]], fast_chg.v_physical_chg[route[idx]-10000]
            A[idx] = D[idx-1] + dist_all[v1, v2]
            # @show((route[idx-1], route[idx]),(v1, v2), dist_all[v1, v2])
        else
            A[idx] = D[idx-1] + dist_orig[route[idx-1], route[idx]]
        end
        B[idx] = max(A[idx], ei_route[idx])
        W[idx] = B[idx] - A[idx]
        D[idx] = B[idx] + s_t_route[idx]
        Load[idx] = Load[idx-1] + qi_route[idx]
    end
    for idx in 2:N-1  # pickup or dropoff
        if route[idx] > n_cus + 1 && (route[idx] < 2 * n_cus + 2)
            idx2 = findfirst(x -> x == route[idx] - n_cus, route)
            # RT[idx] = A[idx] - D[idx2]
            RT[idx] = B[idx] - D[idx2]
        end
    end

    return A, B, W, D, Load, RT
end

# Get the total waiting times when arriving earlier than ei associated with the transit stations of the vehicles
# Given the fact that TW are associated at transit stations, whenever there are earlier arrivals at transit stations,
# the beginning time of service at the pickup location could always be delayed as much as the waiting time, so the 
# buses can always arrive within the desired TW at transit stations.
# Consequently, this function is used for the final best solution check
function get_wait_time_sol(solution::Solution, qi, ei, li, s_t, dist_orig, Li, fast_chg)
    start_depot = 1; end_depot = 2*n_cus+2
    wait_time_sol = 0
    for r in 1:solution.n_route
        route = get_route(solution, r, start_depot, end_depot)
        if  solution.vec_chg_event[r,2] > 0  
            chg_operation_r = solution.vec_chg_event[r,:]   
            wait_time_sol +=  get_wait_time_chg(route, fast_chg, chg_operation_r, dist_orig, dist_all, ei, li, s_t, Li, n_cus) 
        else
            wait_time_sol +=  get_wait_time_no_chg(route, dist_orig, ei, li, s_t, Li, n_cus, qi)
        end
    end    
    return wait_time_sol

end

# check tw constraint
function get_wait_time_chg(route_0::Vector{Int32}, fast_chg::Fast_chargers, chg_operation_r, dist_orig, dist_all, ei, li, s_t, Li, n_cus)
 
    route = copy(route_0)
    ei_route  = ei[route]
    li_route  = li[route]
    s_t_route = s_t[route] 
    Li_route  = Li[route] 
    for i in 1:floor(Int32, chg_operation_r[2])
        idx0 = 2+(i-1)*4
        v1    = floor(Int32, chg_operation_r[idx0+1]) # node precede the charging station (charger) node 
        v_chg= fast_chg.v_chg[floor(Int32, chg_operation_r[idx0+2])]     
        t_chg_start = chg_operation_r[idx0+3]
        dur_chg     = chg_operation_r[idx0+4] - chg_operation_r[idx0+3]
        idx_insert  = findfirst(x->x==v1, route) 
        insert!(route,     idx_insert+1, v_chg)
        insert!(ei_route,  idx_insert+1, t_chg_start) 
        insert!(li_route,  idx_insert+1, li[route[1]]) # planning horizon
        insert!(s_t_route, idx_insert+1, dur_chg)
        insert!(Li_route,  idx_insert+1, li[end])
    end
 
    N=length(route)
    t_start = ei_route[1]
    A, B, W, D, RT = calculate_ads_t_wait_chg(route, t_start, ei_route, s_t_route, dist_orig, dist_all, n_cus, fast_chg) 
    F0 = calc_Fi_chg(1, route, W, B, RT, Li_route, li_route)
    Wp = sum(W[2:N-1])
    t_start_new = ei_route[1] + min(F0, Wp)
    t_wait = calculate_t_wait_chg(route, t_start_new, ei_route, li_route, s_t_route, dist_orig, n_cus, fast_chg)
    if t_wait > 0
        @show("t_wait > 0 in get_wait_time_chg !! ", route, t_wait) 
    end
    return t_wait
end

function get_wait_time_no_chg(route::Vector{Int32}, dist_orig, ei, li, s_t, Li, n_cus, qi)

    N=length(route)
    t_start = ei[route[1]]
    A, B, W, D, Load, RT = calculate_ads(route, t_start, ei, li, s_t, dist_orig, n_cus, qi)
    F0 = calc_Fi(1, route, W, B, RT, Li, li)
    Wp = sum(W[2:N-1])
    t_start_new = ei[route[1]] + min(F0, Wp)
    t_wait = calculate_t_wait_no_chg(route, t_start_new, ei, li, s_t, dist_orig, n_cus)
    if t_wait > 0
        @show("t_wait > 0 !! in get_wait_time_no_chg", route, t_wait) 
    end
    return t_wait  
    
end


#function to check whether the charged energy is sufficient for the route
function check_energy(solution::Solution, darp, route, r, parameter_energy, fast_chg::Fast_chargers)
 
    n_route =solution.n_route
    set_physical, dist_all = darp.set_physical, darp.dist_all
    add_dist = 0 
    for i in 1:floor(Int32, solution.vec_chg_event[r,2])
        v1    = floor(Int32, solution.vec_chg_event[r, 2+(i-1)*4+1])
        idx = findfirst(x->x==v1, route)
        if isnothing(idx)
            @show("check_energy error in findfirst(x->x==v1, route)! ", v1, r, route, solution.vec_chg_event[1:n_route,:] )
        end
        v2= route[idx+1]
        v1, v2   = set_physical[v1], set_physical[v2]
        chg_id   = floor(Int32,solution.vec_chg_event[r, 2+(i-1)*4+2]) # chg id on the list of used fast chargers
        v_physi_chg       =  fast_chg.v_physical_chg[chg_id]
        t_access_1 = dist_all[v1, v_physi_chg] 
        t_access_2 = dist_all[v_physi_chg, v2]
        dist_v1_v2 = dist_all[v1, v2]
        add_dist += (t_access_1 + t_access_2 - dist_v1_v2)
    end
    veh_id = solution.RI[r,5]
    e_consum = parameter_energy.beta[veh_id] * (length_route_energy_check(route, darp) + add_dist) # pay attention here, don't use length_route(route)
    _ , eng_chged = get_route_chg_time_energy(solution, r, fast_chg)
    #  @show("check energy :", r, route,  parameter_energy.E_init[veh_id], e_consum, eng_chged, parameter_energy.E_init[veh_id] + eng_chged - e_consum - parameter_energy.E_min[veh_id])
    if round(e_consum - eng_chged - parameter_energy.max_ec[veh_id], digits = DIGITS) > 0
        @show(r, route, veh_id, parameter_energy.E_init[veh_id], eng_chged, e_consum, parameter_energy.E_min[veh_id], parameter_energy.max_ec[veh_id], solution.vec_chg_event[r,:] )
        # @show(parameter_energy.E_init[veh_id] + eng_chged - e_consum - parameter_energy.E_min[veh_id], e_consum - eng_chged - parameter_energy.max_ec[veh_id])
        error("check energy failed !")
        return false
    else
        # @show("check_energy successful !" )
        return true
    end
end


#function to check whether the charged energy is sufficient for the route
function check_energy_silent(solution::Solution, darp, route, r, parameter_energy, fast_chg::Fast_chargers)
 
    set_physical, dist_all = darp.set_physical, darp.dist_all
    add_dist = 0 
    for i in 1:floor(Int32, solution.vec_chg_event[r,2])
        v1    = floor(Int32, solution.vec_chg_event[r, 2+(i-1)*4+1])
        idx = findfirst(x->x==v1, route)
        if isnothing(idx)
            @show("check_energy error in findfirst(x->x==v1, route)! ", v1, r, route, solution.vec_chg_event[1:K,:] )
        end
        v2= route[idx+1]
        v1, v2   = set_physical[v1], set_physical[v2]
        chg_id   = floor(Int32,solution.vec_chg_event[r, 2+(i-1)*4+2]) # chg id on the list of used fast chargers
        v_physi_chg       =  fast_chg.v_physical_chg[chg_id]
        t_access_1 = dist_all[v1, v_physi_chg] 
        t_access_2 = dist_all[v_physi_chg, v2]
        dist_v1_v2 = dist_all[v1, v2]
        add_dist += (t_access_1 + t_access_2 - dist_v1_v2)
    end
    veh_id = solution.RI[r,5]
    e_consum = parameter_energy.beta[veh_id] * (length_route_energy_check(route, darp) + add_dist) # pay attention here, don't use length_route(route)
    _ , eng_chged = get_route_chg_time_energy(solution, r, fast_chg)
    if round(e_consum - eng_chged - parameter_energy.max_ec[veh_id], digits = DIGITS) > 0
        return false
    else 
        return true
    end
end



function verify_solution(solution::Solution, instance, darp, fast_chg, global_parameter)

    # solution = best_solu
    K, Q = instance.K, instance.Q_max 
    n_cus, start_depot, end_depot, dist = darp.n_cus, darp.start_depot, darp.end_depot, darp.dist
    parameter_energy = instance.parameter_energy
    ei, li, qi, Li, s_t, dist_orig, TH = darp.ei, darp.li, darp.qi, darp.Li, darp.s_t, darp.dist_orig, darp.TH
    count_sol= zeros(Int32,2*n_cus+1)
    infeasible_precedence = 0; infeasible_tw_cap_ride_time = 0; infeasible_energy=0
    n_route= solution.n_route
 
    for r in 1:n_route
        route = get_route(solution, r, start_depot, end_depot)  
        # @show(r, route)
        infeasible_energy += !check_energy(solution, darp, route, r, parameter_energy, fast_chg)
        infeasible_precedence += !check_precedence(route, n_cus)
        infeasible_tw_cap_ride_time += !eight_step(route, r, ei, li, s_t, dist_orig, Q, Li, n_cus, TH, qi)
        for i in 2:length(route)-1
            count_sol[route[i]] +=1 
        end 
    end
    unserved_user = solution.unserved_users
    n_userved_request = length(unserved_user)
    for user in unserved_user
        count_sol[user] +=1 ;  count_sol[user+n_cus] +=1
    end
    solution.co2_emission > global_parameter.co2_threshold ? infeasible_co2emission = true : infeasible_co2emission = false
    if infeasible_co2emission || sum(count_sol) != 2*n_cus || infeasible_precedence > 0 || infeasible_tw_cap_ride_time > 0  || infeasible_energy > 0 || n_userved_request > 0
        @show("verify_solution failed: ",count_sol, infeasible_precedence, infeasible_tw_cap_ride_time, infeasible_energy)
        if infeasible_precedence > 0
            @show("check precedence failed !!", infeasible_precedence)
            for r in 1:n_route
                route = get_route(solution, r, start_depot, end_depot) 
                @show(r, route)
            end
        end
        if n_userved_request >0
            @show("check all requests  served failed !!", n_userved_request, solution.unserved_users)
        end

        infeasible_co2emission && @show("infeasible_co2emission!!", infeasible_co2emission)

        if infeasible_tw_cap_ride_time > 0
            for r in 1:solution.n_route 
                if !verify_eight_step(solution, darp, instance,  r)           
                    route = get_route(solution, r, start_depot, end_depot)
                    A, B, W, D, Load, RT = calculate_ads(route, 0, ei, li, s_t, dist_orig, n_cus, qi)
                    @show("tw check failed : ", r,  route ) 
                    @show(ei[route], li[route], A, B, li[route]-B, W, D, Load, RT, Li[route])                    
                    # route_check = DataFrame(node = route, ei = ei[route], li = li[route], t_arr=A, B=B, li_minus_B = li[route]-B, W=W
                    #                  , D=D, Load=Load, RT=RT, Li=Li[route] )
                    # str1 = path * "route_check_haru.xlsx" 
                    # XLSX.writetable(str1, overwrite=true, route_check)
                end 
            end
        end
        infeasible_energy > 0 && @show("infeasible_energy error : ")
        error("verify_solution failed")
    else
        @show("valid solution")
        # @show(count_sol, infeasible_precedence, infeasible_tw_cap_ride_time, infeasible_energy)
    end
end

function verify_solution_silent(solution::Solution, instance, darp, fast_chg)

    K, Q = instance.K, instance.Q_max 
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
    parameter_energy = instance.parameter_energy
    ei, li, qi, Li, s_t, dist_orig, TH = darp.ei, darp.li, darp.qi, darp.Li, darp.s_t, darp.dist_orig, darp.TH
    count_sol= zeros(Int32,2*n_cus+1)
    infeasible_precedence = 0; infeasible_tw_cap_ride_time = 0; infeasible_energy=0
    n_route= solution.n_route
 
    for r in 1:n_route
        route = get_route(solution, r, start_depot, end_depot) 
        infeasible_energy += !check_energy(solution, darp, route, r, parameter_energy, fast_chg)
        infeasible_precedence += !check_precedence(route, n_cus)
        infeasible_tw_cap_ride_time += !eight_step(route, r, ei, li, s_t, dist_orig, Q, Li, n_cus, TH, qi)
        for i in 2:length(route)-1
            count_sol[route[i]] +=1 
        end
    end
    unserved_user = solution.unserved_users
    for user in unserved_user
        count_sol[user] +=1 ;  count_sol[user+n_cus] +=1
    end
    if sum(count_sol) != 2*n_cus || infeasible_precedence > 0 || infeasible_tw_cap_ride_time > 0  || infeasible_energy > 0
        return false
    else
        return true
    end
end

# eight step fesibility check of DARP (Parragh et al., 2010)
function verify_eight_step(solution::Solution, darp, instance, route_i)

    Q =  instance.Q_max 
    n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot 
    ei, li, qi, Li, s_t, dist, TH = darp.ei, darp.li, darp.qi, darp.Li, darp.s_t, darp.dist, darp.TH
    route = get_route(solution, route_i, start_depot, end_depot)
    return eight_step(route, route_i, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)
  
end
# verify whether the route info are correctly calculated for a solution
function verify_route_info(solution::Solution)

    e_r   = copy(ei)
    cap_r = zeros(Int32,N_NODE)
    l_r   = copy(li)
    maxcap_r = zeros(Int32,N_NODE)
    for r in solution.n_route
        v_i = start_depot
        v_j = solution.RI[r,1]
        while v_j > 0
            e_r[v_j] = max(ei[v_j], e_r[v_i]+s_t[v_i]+dist[v_i,v_j])
            cap_r[v_j] =  cap_r[v_i] + qi[v_j]
            v_i,v_j = v_j, solution.succ[v_j]
        end
        # backward updates
        v_j = end_depot# start from the next node of the drop off node
        v_i = solution.RI[r,2]
        while v_i > 0 
            l_r[v_i] = min(li[v_i], l_r[v_j] - dist[v_j,v_i]-s_t[v_i])
            maxcap_r[v_i] = max(0, maxcap_r[v_j] + qi[v_i])
            v_i,v_j = solution.pre[v_i], v_i
        end
    end
   @show(e_r,cap_r,l_r,maxcap_r )
end
###########################################################################
# wrapper function for call different local search operator using ls_name
########################################################################
function LS(ls_name, solution::Solution, soluT::Solution, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, T, Li) 
    ls_name(solution, soluT, recorder_routes, recorder_lscost_routes, relatedness, demand_data,  fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
          
end

# copy the solution to sol 
function update_sol(sol::Solution,solution::Solution)

    sol.pre  = copy(solution.pre)
    sol.succ = copy(solution.succ)
    sol.RI   = copy(solution.RI)
    sol.unserved_users = copy(solution.unserved_users)
    sol.vec_chg_event = copy(solution.vec_chg_event)
    sol.n_route = solution.n_route
    sol.total_cost = solution.total_cost
    sol.total_cost_with_penalty = solution.total_cost_with_penalty
    sol.penalty_unserved = solution.penalty_unserved
    sol.total_chg_time =  solution.total_chg_time
    sol.co2_emission = solution.co2_emission
    
end

function copy_route_info(e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, _er, _capr, _lr, _maxcapr, _dt_r, _dt_r_back)
    e_r= copy(_er); cap_r = copy(_capr)
    l_r= copy(_lr) ; maxcap_r = copy(_maxcapr)
    dt_r = copy(_dt_r) ; dt_r_back = copy(_dt_r_back)
end

function init_chargers(parameter_energy, n_chg_station, max_n_charger, v_chg_physical, idxs_chg_staion,
                       max_n_chg_by_site, list_chgs_idx_by_site, cost_rapid_chg)
   
    idxs_rapid_chg =collect(1:max_n_charger)
    pw = parameter_energy.pw_charger[idxs_rapid_chg]   #physical nodes for all fast chargers 
    n_fast_chg_installed = 0; num_chg_installed_by_site = zeros(Int32,n_chg_station)
    set_id_chg_installed = Int32[]
    return Fast_chargers(pw, max_n_charger, idxs_rapid_chg, idxs_chg_staion, v_chg_physical, 
              collect(10000+1:10000+max_n_charger), max_n_chg_by_site, list_chgs_idx_by_site,
              n_fast_chg_installed, num_chg_installed_by_site, set_id_chg_installed, cost_rapid_chg, zeros(Int32,2,2))
end 

# compute access time to all the fast chargers  
function compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg::Fast_chargers)
   
    v1, v2     = set_physical[v[1]], set_physical[v[2]]
    chg        = fast_chg.v_physical_chg[idx_chg]
    t_access_1 = dist_all[v1, chg] 
    t_access_2 = dist_all[chg, v2]
    t_v_v1     = dist_all[v1, v2]
    return  t_access_1, t_access_2, t_v_v1
    
end



# get to the nearest DC fast charger 
function get_chg_greedy(v, dist_all, set_physical, fast_chg::Fast_chargers, veh_id, parameter_energy, ec_route)
   
    idxs=  collect(1:fast_chg.n_fast_chg_installed) 
   
    v1, v2      = set_physical[v[1]], set_physical[v[2]]
    chgs        = fast_chg.v_physical_chg[idxs]
    vec_dist_tmp = (dist_all[v1, chgs] + dist_all[chgs, v2]) .- dist_all[v1, v2] 
    vec_delta_e  =  parameter_energy.E_min[veh_id] .- (parameter_energy.E_init[veh_id] .- (ec_route .+ parameter_energy.beta[veh_id] .* vec_dist_tmp))
    vec_t_chg    = vec_delta_e ./ fast_chg.pw[idxs]
    vec_add_chg_op_time= vec_dist_tmp .+ vec_t_chg 
    idx_chg_selected = argmin(vec_add_chg_op_time) 
    # @show( v1, chgs, vec_dist_tmp, vec_delta_e, vec_t_chg, vec_add_chg_op_time, idx_chg_selected)

    return  idx_chg_selected 
end

#schedule a charging event (start time, end time, duration and the energy to be recharged at next charging location of the route)
function schedule_chg(info_chg, darp, parameter_energy, e_charged_route, tt_acc_chg_total, dt_r, remain_e_chg, r, route, veh_id, forward_time_slack, vec_t_start, idx_pre_chg_pos, t_arr_node_after_chg, fast_chg::Fast_chargers, pos_insert_chg, idx_first_loc, count_chg, v1, idx_chg, t_access_1, t_access_2, t_v_v1)
    
    e_charged_tmp = t_chg_constr= 0
    add_access_chg_time = t_access_1+ t_access_2 - t_v_v1
    if count_chg > 1
        remain_e_chg += (add_access_chg_time * parameter_energy.beta[veh_id])
    end
    idx_loc = idx_first_loc - 1 + count_chg
    idx_pos = pos_insert_chg[idx_loc,3]
    e_actual = parameter_energy.E_init[veh_id] - (dt_r[v1] + t_access_1) * parameter_energy.beta[veh_id]
    e_actual += (e_charged_route - tt_acc_chg_total * parameter_energy.beta[veh_id])
    t_chg_max = (parameter_energy.E_max[veh_id] - e_actual) / fast_chg.pw[idx_chg]

    if count_chg > 1 # update forward_time_slack calculation
        W, B = calculate_ads_chg_schedule(route[idx_pre_chg_pos+1:end], t_arr_node_after_chg, darp)  # re_calculate A,B,W,D after the precedent chg event
        idx_Fi = idx_pos - idx_pre_chg_pos
        forward_time_slack_new, t_start_fw_slack = get_forward_time_slack_revise(route[idx_pre_chg_pos+1:end], idx_Fi, W, B, darp.li)
        tmp          = max(0, forward_time_slack_new - add_access_chg_time)
        t_chg_constr = min(tmp, t_chg_max)
    else
        tmp = max(0, forward_time_slack[idx_loc] - add_access_chg_time,0)
        t_chg_constr = min(tmp, t_chg_max)
        t_start_fw_slack = vec_t_start[idx_loc] 
    end
    
    t_chg_desired = remain_e_chg / fast_chg.pw[idx_chg]
    t_chg= min(t_chg_constr, t_chg_desired)
    (t_chg <0.01) && (t_chg =0)
    additional_time = add_access_chg_time +t_chg 
    info_chg[1] += additional_time;  info_chg[2] +=1
    idxs = collect(2+(count_chg-1)*4+1 : 2+count_chg*4)
    # @show(forward_time_slack[count_chg] , add_access_chg_time,t_chg_constr,t_chg_desired)
    # @show(count_chg, t_chg_constr, t_chg_desired, t_chg )
    if  t_chg_desired > t_chg_constr
        t_chg_start = t_start_fw_slack + t_access_1
        t_chg_end = t_chg_start + t_chg
        info_chg[idxs] = [v1, idx_chg, t_chg_start, t_chg_end]
        remain_e_chg -= (t_chg * fast_chg.pw[idx_chg])
        t_arr_node_after_chg = t_chg_end + t_access_2 
        if t_chg > 0 
            e_charged_tmp        = t_chg * fast_chg.pw[idx_chg]
            return true, remain_e_chg, t_arr_node_after_chg, t_chg, e_charged_tmp, add_access_chg_time   #continue to recharge
        else
            return false, 0, 0 , t_chg, 0, 0 # no feasible , t_chg=0
        end
    else
        # @show(remain_e_chg=0)
        t_chg_start = t_start_fw_slack + t_access_1 + max(0, t_chg_constr - t_chg) * rand() # random start time
        t_chg_end = t_chg_start + t_chg
        info_chg[idxs] = [v1, idx_chg, t_chg_start, t_chg_end]
        e_charged_tmp  = t_chg * fast_chg.pw[idx_chg]
        return false, 0, 0, t_chg, e_charged_tmp, add_access_chg_time   # charging operation finished
    end 
end




# insert charging events on a route. First, identify possible charging positions of the route, then insert randomly the charging event 
# on a random position and a random (random policy)/best(greedy policy) charger 
function insert_charging_route(solution::Solution, darp, r, route::Vector{Int32}, dist_all, ec_route, fast_chg::Fast_chargers, parameter_energy, dt_r, dt_r_back, flag_init_sol)
                            
    degree_rand_chg_policy = parameter_energy.degree_rand_chg_policy
    v=Int32[]; idx_chg = t_access_1 = t_access_2 = t_v_v1 = vec_delta_e = t_chg = 0 # some are dummy declarations
    dist_orig, n_cus, end_depot = darp.dist_orig, darp.n_cus, darp.end_depot 
    ei, li, s_t, dist, qi = darp.ei, darp.li, darp.s_t, darp.dist, darp.qi
    set_physical, dist_all = darp.set_physical, darp.dist_all
 
    # walk_tt_bus_stop = darp.walk_tt_bus_stop
    max_num_chg = parameter_energy.max_num_chg # max num of charging events limits for a route 
    success_chg= false; remain_e_chg =0; idx_loc =0; charging_time=0
    n_fast_chg = fast_chg.n_fast_chg_installed
    pos_insert_chg = zeros(Int32, 30, 3) # store the pair nodes (precedent and successive nodes) of the location of the charger
    n_pos_insert_chg = 1; t_v_v1 = 0; t_arr_node_after_chg = 0; idx_pre_chg_pos =0
    # find the positions (nodes) where the SOC (state of charge) of the vehicle is insufficient  
    pos_insert_chg[n_pos_insert_chg,:] = [route[1], route[2], 1]
    veh_id = solution.RI[r,5]
    max_ec = parameter_energy.max_ec[veh_id]
    i_ec_violate= 0 # will alwayse find a positive position 
    dt_r_eng = dt_r[route] .* parameter_energy.beta[veh_id]
    i_ec_violate = findfirst(x->x-max_ec>0, dt_r_eng)
    e_charged_route = tt_acc_chg_total = 0 #get the charged energy and total access times for the chg operations of a route
    # @show("insert_charging_route: ", r, veh_id, max_ec, i_ec_violate)
    # need to check whether all customers on board have beed drop-offed to be a legitimate candidate charging position
    for (idx, v) in enumerate(route[1:i_ec_violate-1])
        if  v > n_cus+1 && (route[idx+1] < n_cus+2 || route[idx+1] == end_depot)
            if (idx-1)%2 ==0 # check if all customers on board have beed drop-offed
                n_pos_insert_chg += 1
                pos_insert_chg[n_pos_insert_chg, :] = [v, route[idx+1], idx]
            end
        end
    end  
    #randomly select a fast charger
    idx_first_loc = 1
    if rand() > 0.75 && (n_pos_insert_chg > 1) # 25% start from random position
       idx_first_loc =  rand(collect(2:n_pos_insert_chg))
    end
    # chg is the fast charger index in the list of all fast chargers (vec_idx_chg)
    v= pos_insert_chg[idx_first_loc,1:2]
    idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg, veh_id, parameter_energy, ec_route)
    t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
    # veh has no sufficient energy to arrive the charger's location, so it is infeasible
    if (dt_r[v[1]] + t_access_1) * parameter_energy.beta[veh_id] > parameter_energy.max_ec[veh_id] && idx_first_loc > 1
        idx_first_loc -= 1
        v= pos_insert_chg[idx_first_loc,1:2]
        t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
    end 
    vec_delta_e  =  parameter_energy.E_min[veh_id] - (parameter_energy.E_init[veh_id] - (ec_route + parameter_energy.beta[veh_id] * (t_access_1+ t_access_2 - t_v_v1)))
    t_chg = vec_delta_e / fast_chg.pw[idx_chg]
    # veh cannot charge to higher than its max
    t_chg_max = (parameter_energy.E_max[veh_id] - parameter_energy.E_init[veh_id] ) / fast_chg.pw[idx_chg] # 
    t_chg > t_chg_max && (idx_first_loc = 1) # if is the case, start charging from the first chg candidate location

    idxs_chg_loc = pos_insert_chg[1:n_pos_insert_chg, 3]
    # @show(idx_loc, v, t_chg, idxs_chg_loc)
    forward_time_slack, vec_t_start = get_forward_time_slack(route, ei, li, s_t, dist, n_cus, qi, idxs_chg_loc, dist_orig, darp)
    # @show(forward_time_slack)
    count_chg = 1; next_recharge = true; 
    n_candicate_loc =  n_pos_insert_chg - idx_first_loc +1
    # @show(pos_insert_chg, idxs_chg_loc)
    remain_e_chg = vec_delta_e # desired charged energy at the beginning
    
    info_chg = zeros(Float32, 2+max_num_chg*4) # store at most 4 charging operations  
    if t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[idx_first_loc] && (t_chg_max >= t_chg )
        # @show("t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[1] !")
        additional_time = t_access_1 + t_access_2 - t_v_v1 +  t_chg 
        # t_start_chg = vec_t_start[1] + t_access_1 + max(0, (forward_time_slack[idx_first_loc]-t_chg))*rand() # random start time
        delta_t = rand() *(forward_time_slack[idx_first_loc] - ( t_chg + t_access_1 + t_access_2 - t_v_v1))
        t_start_chg = vec_t_start[idx_first_loc] + t_access_1 +  delta_t # random start time
        info_chg[1] += additional_time;  info_chg[2] +=1
        info_chg[3:6] = [v[1], idx_chg, t_start_chg, t_start_chg+t_chg]
        # e_charged_route  += (t_chg * fast_chg.pw[idx_chg])
        # tt_acc_chg_total += (t_access_1 + t_access_2 - t_v_v1)
        # soc_af_chg = parameter_energy.E_init[veh_id] + e_charged_route
        # eng_consum = ec_route + parameter_energy.beta[veh_id] * tt_acc_chg_total 
        # if soc_af_chg - eng_consum >= parameter_energy.E_min[veh_id]
            success_chg = true 
        # else
            # success_chg = false
        # end
    else        
        while next_recharge  
            if count_chg >= n_candicate_loc  || count_chg >= max_num_chg
                success_chg = false
                break
            else
                # @show(count_chg,v[1])
                next_recharge, remain_e_chg, t_arr_node_after_chg, charging_time, e_charged_tmp, add_access_chg_time =  schedule_chg(info_chg, darp, parameter_energy, e_charged_route, tt_acc_chg_total, dt_r, remain_e_chg, r, route, veh_id, forward_time_slack, vec_t_start, idx_pre_chg_pos, t_arr_node_after_chg, fast_chg, pos_insert_chg, idx_first_loc, count_chg, v[1], idx_chg, t_access_1, t_access_2, t_v_v1)
                e_charged_route += e_charged_tmp; tt_acc_chg_total += add_access_chg_time
                idx_loc_pre = idx_first_loc - 1 + count_chg
                idx_pre_chg_pos = pos_insert_chg[idx_loc_pre, 3]
                # next location
                count_chg += 1
                if count_chg < n_candicate_loc +1
                    v = pos_insert_chg[idx_first_loc-1+count_chg,1:2]   
                    idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg, veh_id, parameter_energy, ec_route) 
                    t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
                end
            end
        end 
        # ! next_recharge && (success_chg = true)
        if !next_recharge && charging_time>0        
            # soc_af_chg = parameter_energy.E_init[veh_id] + e_charged_route
            # eng_consum = ec_route + parameter_energy.beta[veh_id] * tt_acc_chg_total 
            # if soc_af_chg - eng_consum >= parameter_energy.E_min[veh_id]
                success_chg = true
            # else
            #     success_chg = false
            # end
        else
            success_chg = false
        end
    end
    if success_chg              
        solution.vec_chg_event[r,:] =  info_chg
        # check_energy(solution, darp, route, r, parameter_energy, fast_chg)
        return true,  pos_insert_chg, forward_time_slack, vec_t_start
    else     
        return false, pos_insert_chg, forward_time_slack, vec_t_start
    end
end


# swap electric to gasoline if have cost saving (based on the input setting)
function veh_exchange(solution::Solution, darp, fast_chg, parameter_energy, instance, global_parameter)
    # @show("veh_exchange call!")
    start_depot = darp.start_depot 
    n_route = solution.n_route
    co2_threshold = global_parameter.co2_threshold 
    cost_saving_best = 0; route_id_exchanged = 0
    veh_ids_types    = parameter_energy.veh_ids_types
    for r in 1:n_route 
        veh_id= solution.RI[r,5]  # type 1 is gazo, type 2 is electric
        if parameter_energy.is_electric[veh_id]   
            if solution.co2_emission + co2_route_gv(solution, r, start_depot, darp, parameter_energy) <= co2_threshold
                cost_saving = cost_saving_exchange(solution, r, darp, parameter_energy, fast_chg)
                if cost_saving > cost_saving_best
                    cost_saving_best = cost_saving
                    route_id_exchanged = r 
                end
            end
        end 
    end
    if route_id_exchanged > 0 
        veh_used = solution.RI[1:n_route, 5]
        vehs = setdiff(veh_ids_types[1], veh_used) # set of unused agso vehs
        veh_id = vehs[1]  
        solution.RI[route_id_exchanged,4] = parameter_energy.cap_passenger[veh_id]
        solution.RI[route_id_exchanged,5] = veh_id  
        solution.vec_chg_event[route_id_exchanged,:] .*= 0  
        cost_solution(solution, global_parameter, darp, parameter_energy, fast_chg, instance) 
    end
    # # exchange from gasoline to ev
    # if solution.co2_emission > co2_threshold 
    #     for r in 1:n_route 
    #         veh_id= solution.RI[r,5]  # type 1 is gazo, type 2 is electric
    #         if !parameter_energy.is_electric[veh_id]  
    #             solution.co2_emission - co2_route_gv(solution, r, start_depot, darp, parameter_energy) <= co2_threshold
    #             veh_used = solution.RI[1:n_route, 5]
    #             vehs = setdiff(veh_ids_types[1], veh_used) # set of unused agso vehs
    #             veh_id = vehs[1]  
    #             solution.RI[r,4] = parameter_energy.cap_passenger[veh_id]
    #             solution.RI[r,5] = veh_id  
    #             solution.vec_chg_event[route_id_exchanged,:] .*= 0  
    #             cost_solution(solution, global_parameter, darp, parameter_energy, fast_chg, instance) 
    #             break
    #         end
    #     end
    # end
end

# function insert_chg_af_exchange()

#     route = get_route(soluT, r_i, start_depot, end_depot)
#     n_charger_installed, success, pos_insert_chg, forward_time_slack, vec_t_start =update_charging(soluT, darp, r_i, route, dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
#     if ! success && (n_charger_installed > 0)                            
#         repair_charging(soluT, darp, r_i, route, fast_chg, parameter_energy, dt_r, pos_insert_chg, forward_time_slack, vec_t_start)
#     end
      
# end

# insert charging operation at a random location
function update_charging(solution::Solution, darp, route_i, route::Vector{Int32}, dist_all, fast_chg::Fast_chargers, parameter_energy, dt_r, dt_r_back)
    
    flag_init_sol = false
    n_cus, dist_orig = darp.n_cus, darp.dist_orig
    end_depot = darp.end_depot
    veh_id = solution.RI[route_i,5]  
    n_charger_installed = fast_chg.n_fast_chg_installed
    if parameter_energy.is_electric[veh_id] == true 
        compute_dt_r(route, dt_r, dist_orig)  
        ec_route =  dt_r[end_depot] * parameter_energy.beta[veh_id]
        if  ec_route > parameter_energy.max_ec[veh_id] # avoid using the end_depot to check it as the dt_r[end_depot] might be changed by other routes
            if n_charger_installed > 0
                a, b, c, d = insert_charging_route(solution, darp, route_i, route, dist_all, ec_route, fast_chg, parameter_energy, dt_r, dt_r_back, flag_init_sol)
                return n_charger_installed, a, b, c, d
            else
                @warn("n_charger_installed = 0 !",veh_id, ec_route, parameter_energy.max_ec[veh_id])
                return n_charger_installed, false, 0,0,0
            end
        else    
            solution.vec_chg_event[route_i,:] .*= 0
            return n_charger_installed, true, 0, 0, 0
        end
    else
        solution.vec_chg_event[route_i,:] .*= 0
        return n_charger_installed, true, 0, 0, 0
    end
    
end

# repair the charging operations by inserting charging operations starting from the depot
function repair_charging(solution::Solution, darp, r, route::Vector{Int32}, fast_chg::Fast_chargers, parameter_energy, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start)
    
    degree_rand_chg_policy = parameter_energy.degree_rand_chg_policy
    v=Int32[]; idx_chg=0
    success_chg= false;charging_time=0
    set_physical, dist_all = darp.set_physical, darp.dist_all
    max_num_chg =  parameter_energy.max_num_chg # max num of charging events on a route
    n_fast_chg = fast_chg.n_fast_chg_installed
    n_pos_insert_chg = length(forward_time_slack); t_v_v1 = 0; t_arr_node_after_chg =0; idx_pre_chg_pos = 0
    info_chg = zeros(Float32, 2+max_num_chg*4) # store at most 4 charging operations  
    veh_id = solution.RI[r,5]
    ec_route =  dt_r[darp.end_depot] * parameter_energy.beta[veh_id]
    e_charged_route = tt_acc_chg_total = 0 #get the charged energy and total access times for the chg operations of a route
    compute_dt_r(route, dt_r, darp.dist_orig)       
    # rand() > degree_rand_chg_policy ? rand_policy =true : rand_policy = false
    v= pos_insert_chg[1,1:2]; idx_loc =1 # start charging operation from the depot

    idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg, veh_id, parameter_energy, ec_route)
  
    t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
    
    vec_delta_e  = parameter_energy.E_min[veh_id] - (parameter_energy.E_init[veh_id] - (ec_route + parameter_energy.beta[veh_id] * (t_access_1+ t_access_2 - t_v_v1)))
    # vec_delta_e <= 0 && error("vec_delta_e <=0 error!!")
    t_chg = vec_delta_e / fast_chg.pw[idx_chg]
    t_chg_max = (parameter_energy.E_max[veh_id] - parameter_energy.E_init[veh_id] ) / fast_chg.pw[idx_chg] # 
    t_chg > t_chg_max && (idx_first_loc = 1) # if is the case, start charging from the first chg candidate location

    count_chg = 1; next_recharge = true; idx_first_loc =1 
    # @show(pos_insert_chg, idxs_chg_loc)
    remain_e_chg = vec_delta_e
    if t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[1] && (t_chg_max >= t_chg )
        additional_time = t_access_1 + t_access_2 - t_v_v1 +  t_chg 
        # t_start_chg = vec_t_start[1] + t_access_1 + max(0, (forward_time_slack[1] - t_chg))*rand() # random start time
        delta_t = rand() *(forward_time_slack[1] - ( t_chg + t_access_1 + t_access_2 - t_v_v1))
        t_start_chg = vec_t_start[1] + t_access_1 + delta_t # random start time
        info_chg[1] += additional_time;  info_chg[2] +=1
        info_chg[3:6] = [v[1], idx_chg, t_start_chg, t_start_chg+t_chg]
        # e_charged_route  += (t_chg * fast_chg.pw[idx_chg])
        # tt_acc_chg_total += (t_access_1 + t_access_2 - t_v_v1)
        # soc_af_chg = parameter_energy.E_init[veh_id] + e_charged_route
        # eng_consum = ec_route + parameter_energy.beta[veh_id] * tt_acc_chg_total 
        # if soc_af_chg - eng_consum >= parameter_energy.E_min[veh_id]
           success_chg = true
        # else
        #     success_chg = false
        # end
    else        
        while next_recharge  
            if count_chg >= n_pos_insert_chg || count_chg >= max_num_chg
                success_chg = false
                break
            else
                next_recharge, remain_e_chg, t_arr_node_after_chg, charging_time, e_charged_tmp, add_access_chg_time =  schedule_chg(info_chg, darp, parameter_energy, e_charged_route, tt_acc_chg_total, dt_r, remain_e_chg, r, route, veh_id, forward_time_slack, vec_t_start, idx_pre_chg_pos, t_arr_node_after_chg, fast_chg, pos_insert_chg, idx_first_loc, count_chg, v[1], idx_chg, t_access_1, t_access_2, t_v_v1)
                e_charged_route += e_charged_tmp; tt_acc_chg_total += add_access_chg_time
                idx_loc_pre = idx_first_loc - 1 + count_chg
                idx_pre_chg_pos = pos_insert_chg[idx_loc_pre, 3]
                count_chg += 1
                if count_chg < n_pos_insert_chg +1
                    v = pos_insert_chg[count_chg,1:2] 
                    idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg, veh_id, parameter_energy, ec_route)
                    t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
                end
            end
        end
        if !next_recharge && charging_time>0
            # soc_af_chg = parameter_energy.E_init[veh_id] + e_charged_route
            # eng_consum = ec_route + parameter_energy.beta[veh_id] * tt_acc_chg_total 
            # if soc_af_chg - eng_consum >= parameter_energy.E_min[veh_id]
                success_chg = true
            # else
            #     success_chg = false
            # end
        else
            success_chg = false
        end
    end
    if success_chg              
        solution.vec_chg_event[r,:] = info_chg
        # @show(" test 24: ", success_chg, r, route, info_chg ) 
        # check_energy(solution, darp, route, r, parameter_energy, fast_chg)              
    end
    return success_chg
end


# function to compute costs of the routes of the solution for each type of vehicle
function compute_matrix_cost_v_type(solution::Solution, darp, parameter_energy)

    n_route    = solution.n_route
    veh_info   = parameter_energy.veh_info
    n_type_veh = parameter_energy.veh_info.n_type_veh
    if n_route == 1        
        matrix_cost_v_type = zeros(n_type_veh)
        matrix_co2_v_type  = zeros(n_type_veh) # the first type is gaso 
        tt_route = length_route(solution, 1, darp)       
        veh_id= solution.RI[1,5] 
        for j in 1:n_type_veh
            matrix_cost_v_type[j] = veh_info.cost_energy_type[j] * tt_route + veh_info.daily_purchase_cost_type[j]
            matrix_co2_v_type[j] =  veh_info.CO2_type[j] * tt_route
            if parameter_energy.veh_type[veh_id] != j       
                check_cap_r(solution, darp, 1, parameter_energy, j, matrix_cost_v_type) 
            end
        end
    else
        matrix_cost_v_type = zeros(n_route, n_type_veh)
        matrix_co2_v_type  = zeros(n_route, n_type_veh) # the first type is gaso 
        for r in 1:n_route
            tt_route = length_route(solution, r, darp)       
            veh_id= solution.RI[r,5] 
            for j in 1:n_type_veh
                matrix_cost_v_type[r,j] = veh_info.cost_energy_type[j] * tt_route + veh_info.daily_purchase_cost_type[j]
                matrix_co2_v_type[r,j] =  veh_info.CO2_type[j] * tt_route
                if parameter_energy.veh_type[veh_id] != j       
                    check_cap_r(solution, darp, r, parameter_energy, j, matrix_cost_v_type) 
                end
            end
        end 
    end
    return matrix_cost_v_type, matrix_co2_v_type
end



function reset_veh_ids(soluT::Solution, parameter_energy, x_opt, Ni, Nj)

    n_route         = soluT.n_route
    veh_ids_types   = parameter_energy.veh_ids_types
    veh_ids_used_be = soluT.RI[1:n_route, 5]
    for j in Nj # Nj is the types of veh
        veh_ids = veh_ids_types[j] # set of veh ids of the type j
        count_tmp = 0
        for r in Ni
            if  x_opt[r,j] == 1 
                count_tmp +=1 
                veh_id = veh_ids[count_tmp]
                soluT.RI[r, 5] = veh_id
                soluT.RI[r, 4] =  parameter_energy.cap_passenger[veh_id] 
            end
        end
    end   
    K_max =size(soluT.RI)[1]
    vehs_unused = setdiff(collect(1:K_max), soluT.RI[1:n_route,5])
    soluT.RI[n_route+1:end, 5] = vehs_unused
    soluT.RI[n_route+1:end, 4] = parameter_energy.cap_passenger[vehs_unused]
    veh_ids_used_af = soluT.RI[1:n_route, 5]
    vec_tmp = setdiff(veh_ids_used_be, veh_ids_used_af)
    if isnothing(vec_tmp)
        return false
    else
        return true # successfully exchanged
    end
end

function update_chaging_exchange(solution::Solution, soluT::Solution, global_parameter, instance, darp, fast_chg, parameter_energy, dt_r, dt_r_back)

    n_route = soluT.n_route
    dist_all = darp.dist_all
    start_depot, end_depot = darp.start_depot, darp.end_depot 
    success_chg_sol = true
    for r_i in collect(1:n_route)
        route = get_route(soluT, r_i, start_depot, end_depot)
        n_charger_installed, success, pos_insert_chg, forward_time_slack, vec_t_start = update_charging(soluT, darp, r_i, route, dist_all, fast_chg, parameter_energy, dt_r, dt_r_back)
        if ! success    
            if n_charger_installed == 0 
                success_chg_sol = false; break 
            else        
                success = repair_charging(soluT, darp, r_i, route, fast_chg, parameter_energy, dt_r, pos_insert_chg, forward_time_slack, vec_t_start)
                if ! success
                    success_chg_sol = false; break
                end 
            end
        end
    end
    if success_chg_sol 
        cost_solution(soluT, global_parameter, darp, parameter_energy, fast_chg, instance)
        update_sol(solution, soluT) 
        # @show("succ exchanged :", solution.total_cost, solution.total_cost_with_penalty, solution.RI[1:n_route,5])
    end     
end
 

# swap electric to gasoline if have cost saving (based on the input setting)
function veh_exchange_milp(solution::Solution, soluT::Solution, instance, darp, fast_chg, parameter_energy, global_parameter, dt_r, dt_r_back)
    
    n_route = solution.n_route
    update_sol(soluT, solution)
    n_type_veh                 = parameter_energy.veh_info.n_type_veh
    solve_veh_exchange_t_limit = global_parameter.solve_veh_exchange_t_limit
    co2_threshold              = global_parameter.co2_threshold 
    veh_type                   = parameter_energy.veh_type
    matrix_cost_v_type, matrix_co2_v_type = compute_matrix_cost_v_type(soluT, darp, parameter_energy) # veh_ids * n_vec_type
    if n_route == 1
        idxs= findall(x->x<=co2_threshold, matrix_co2_v_type)
        idx_min = argmin(matrix_cost_v_type[idxs])
        id_type_f= idxs[idx_min]
        if id_type_f != veh_type[soluT.RI[1,5]] #veh_type is changed
            veh_id = veh_ids_types[id_type_f][1] # set of veh ids of the type j
            soluT.RI[r, 5] = veh_id
            soluT.RI[r, 4] =  parameter_energy.cap_passenger[veh_id] 
            K_max =size(soluT.RI)[1]
            vehs_unused = setdiff(collect(1:K_max), soluT.RI[1,5])
            soluT.RI[n_route+1:end, 5] = vehs_unused
            soluT.RI[n_route+1:end, 4] = parameter_energy.cap_passenger[vehs_unused]
            update_chaging_exchange(solution, soluT, global_parameter, instance, darp, fast_chg, parameter_energy, dt_r, dt_r_back)
          
        end
    else    # veh_exchange milp formulation
        Ni =collect(1:n_route); Nj = collect(1:n_type_veh)
        
        model = direct_model(Gurobi.Optimizer(GRB_ENV))
        set_optimizer_attribute(model, "OutputFlag", 0)
        set_optimizer_attribute(model, "Timelimit", solve_veh_exchange_t_limit)
        set_silent(model)
        
        @variable(model, x[1:n_route, 1:n_type_veh], Bin) 
        @objective(model, Min, sum(matrix_cost_v_type[i, j] * x[i, j] for i in Ni for j in Nj))        
        @constraint(model, [i in Ni], sum(x[i, j] for j in Nj) == 1) 
        @constraint(model, sum(matrix_co2_v_type[i,j] * x[i, j] for i in Ni for j in Nj)  <= co2_threshold)
        
        optimize!(model)
 
        if termination_status(model) == MOI.OPTIMAL  
            x_opt = round.(Int, value.(x)) 
            flag_exchanged =  reset_veh_ids(soluT , parameter_energy, x_opt, Ni, Nj) # reset the solution for the type and capacity
            if flag_exchanged
                update_chaging_exchange(solution, soluT, global_parameter, instance, darp, fast_chg, parameter_energy, dt_r, dt_r_back)
            end
        else
            @show(termination_status(model)) 
            @error("veh_exchange_milp canoot find solution, check termination_status(model) !!") 
        end
    end
end


# # exchange two routes(vehices) with charging operations to save charging times 
# # Given the E_init of vehicles are different, there might have saving to change vehicles to save charging times
# function route_exchange(solution::Solution, instance, cap_r, dist_all, Q, fast_chg, parameter_energy, set_physical, nodes, dt_r, dt_r_back, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval)
    
#     start_depot = 1; end_depot = 2*n_cus+2
#     n_route = solution.n_route
#     if solution.n_route > 1
#         r_e_init_sorted = sortperm(parameter_energy.E_init[1:n_route], rev=true)# sort routes in decreasing order of E_init
#         r_chg_sorted    = sortperm(solution.vec_chg_event[1:n_route,1], rev = true) # sort routes in decreasing order of charging times
#         if solution.vec_chg_event[r_chg_sorted[1],1] > 0 && sum(abs.(r_e_init_sorted.-r_chg_sorted)) > 0
#             tabu_list=falses(K)
#             for (idx_1, r1) in enumerate(r_chg_sorted)
#                 tabu_list[r1]= true
#                 route_1= get_route(solution, r1, start_depot, end_depot)
#                 if solution.vec_chg_event[r1,1] > 0  
#                     for (idx_2, r2) in enumerate(r_e_init_sorted)
#                         if  parameter_energy.E_init[r2] > parameter_energy.E_init[r1] 
#                             if ! tabu_list[r2]  
#                                 route_2= get_route(solution, r2, start_depot, end_depot)
#                                 success, chg_time_saving =  exchange_feasible(solution, instance, r1, r2, route_1, route_2, cap_r,dist_all, Q, fast_chg, parameter_energy, set_physical, nodes, dt_r, dt_r_back, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval)
#                                 if success# update route info
#                                     tabu_list[r2] = true
#                                     solution.RI[r1,1:3], solution.RI[r2,1:3] = solution.RI[r2,1:3], solution.RI[r1,1:3]
#                                     # solution.TT_walk_routes[r1],solution.TT_walk_routes[r2] = solution.TT_walk_routes[r2], solution.TT_walk_routes[r1]
#                                     solution.total_cost -= chg_time_saving
#                                     break
#                                 end
#                             end
#                         else
#                             break
#                         end
#                     end 
#                 else
#                     break
#                 end
#             end
#         end
#     end
# end

# using route_i to serve customers on the exchanged route_j (route below)
function update_charging_route_exchange(solution::Solution, instance, darp, route_i, route::Vector{Int32}, dist_all, fast_chg::Fast_chargers, parameter_energy, set_physical, dt_r, dt_r_back)
    
    flag_init_sol = false
    n_cus, dist_orig = darp.n_cus, darp.dist_orig
    end_depot = 2*n_cus+2    
    n_chg_installed = fast_chg.n_fast_chg_installed
    compute_dt_r(route, dt_r, dist_orig)  
    veh_id = solution.RI[route_i,5]  
    ec_route =  dt_r[end_depot] * parameter_energy.beta[veh_id]
    if  ec_route > parameter_energy.max_ec[veh_id] # avoid using the end_depot to check it as the dt_r[end_depot] might be changed by other routes
        return n_chg_installed, insert_charging_route(solution, darp, route_i, route, dist_all, ec_route, fast_chg, parameter_energy, dt_r, dt_r_back, flag_init_sol)
    else    
        solution.vec_chg_event[route_i,:] .*= 0
        return n_chg_installed, true, 0, 0, 0
    end
end

function exchange_feasible(solution::Solution, instance, r1, r2, route_1, route_2, cap_r,dist_all, Q, fast_chg, parameter_energy, set_physical, nodes, dt_r, dt_r_back, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval)
    
    cap_tmp = 0; route=Int32[] 
    solution.RI[r1,4] == solution.RI[r2,4] ? feasible_capacity = true : feasible_capacity = false 
    dur_chg_init= sum(solution.vec_chg_event[1:K,1]) 
    if ! feasible_capacity 
        if solution.RI[r1,4] > solution.RI[r2,4]
            route = copy(route_1)
            cap_tmp = solution.RI[r2,4] 
        else
            route = copy(route_2)
            cap_tmp = solution.RI[r1,4] 
        end
        idxs = findfirst(x->x>cap_tmp, cap_r[route[2:end-1]])
        if isnothing(idxs)
            feasible_capacity = true
        end         
    end
    # re-schedule charging operations
    success_r1 = success_r2 = false
    if feasible_capacity
        vec_chg_event = copy(solution.vec_chg_event) 
        # update charging operations on the exchanged route r1 
        n_chg_installed, success_r1, pos_insert_chg, forward_time_slack, vec_t_start = update_charging_route_exchange(solution, instance, darp, r1, route_2, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back)
        if !success_r1 && n_chg_installed > 0
            success_r1= repair_charging(solution, darp, r1, route_2, fast_chg, parameter_energy, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start)                               
        end
         # update charging operations on the exchanged route r2 
         _, success_r2, pos_insert_chg, forward_time_slack, vec_t_start = update_charging_route_exchange(solution, instance, darp, r2, route_1, dist_all, fast_chg, parameter_energy, set_physical, dt_r, dt_r_back)
        if !success_r2 && n_chg_installed > 0
            success_r2= repair_charging(solution, darp, r2, route_1, fast_chg, parameter_energy, dt_r, pos_insert_chg, forward_time_slack, vec_t_start)                               
        end
        if !(success_r1 * success_r2)  
            solution.vec_chg_event = vec_chg_event # restore charging operations information
            return false, 0
        else 
            chg_time_saving  = dur_chg_init - sum(solution.vec_chg_event[1:K,1])
            if chg_time_saving > 0 && check_chg_occ_constr(solution, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval) 
                return true, chg_time_saving
            else
                solution.vec_chg_event = vec_chg_event # restore charging operations information
                return false, 0
            end
        end
    else
        return false, 0 
    end
end


function SA(solution::Solution, soluT::Solution, solu_init::Solution, demand_data, instance, fast_chg, darp, global_parameter, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, avg_dist_between_nodes, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval)
  
    # solution = solu
    step_size_create_route = global_parameter.step_size_create_route
    penalty, step_size_remove_route = global_parameter.penalty, global_parameter.parameter_sa.step_size_remove_route
    K, parameter_energy, Q = instance.K, instance.parameter_energy, instance.Q_max
    co2_threshold = global_parameter.co2_threshold
    start_depot, Li, n_cus= darp.start_depot, darp.Li, darp.n_cus
    # init a route container for greedy insert etc.
    size_container = 100 # 30 veh, if not enough, increase it
    recorder_routes =  Vector{Container_route}(undef, size_container) # added on 1.9.2023
    container_route = Container_route(0, zeros(2), 0, Int32[], Int32[])
    for i in collect(1:size_container)
        container_route.r = i
        recorder_routes[i]=deepcopy(container_route)
    end
    # to store the least cost, regret, route_id and position to insert a request for all unserved customers
    recorder_lscost_routes =  Vector{Container_route}(undef, n_cus+1) #user id is from 2:n_cus+1
    for i in collect(1:n_cus+1) 
        recorder_lscost_routes[i]=deepcopy(container_route)
    end
    #compute relatedness:  the first is the depot, customer is from 2:n_cus+1
    rel_dist =zeros(Float32, n_cus+1, n_cus+1); rel_tw = zeros(Float32, n_cus+1, n_cus+1) 
    rel_shaw = zeros(Float32, n_cus+1, n_cus+1) # use li
    relatedness = Relatedness(rel_dist, rel_tw, rel_shaw)
    compute_relatedness(darp, global_parameter, relatedness)

    # SA parameter
    T_red  = global_parameter.parameter_sa.T_red
    t_max  = global_parameter.parameter_sa.t_max
    n_imp  = global_parameter.parameter_sa.n_imp
    N_ITER = global_parameter.parameter_sa.N_ITER
    max_stagnant_multiplier = global_parameter.parameter_sa.max_stagnant_multiplier

    sol_best = deepcopy(solution)  # x_best
    sol      = deepcopy(solution)  # x, x' => solution 
    copy_e_r   = copy(e_r);     copy_cap_r = copy(cap_r)
    copy_l_r   = copy(l_r);     copy_maxcap_r = copy(maxcap_r)
    copy_dt_r  = copy(dt_r);    copy_dt_r_back = copy(dt_r_back)
    best_e_r   = copy(e_r);     best_cap_r = copy(cap_r)
    best_l_r   = copy(l_r);     best_maxcap_r = copy(maxcap_r)
    best_dt_r = copy(dt_r);     best_dt_r_back = copy(dt_r_back)
    
    T_max = t_max *avg_dist_between_nodes
    T    = T_max
      
    i_imp=0 ; 

    ls_name= [relocate_ensemble, destroy_repair, two_opt, four_opt_intra, two_opt_star, swap_two_users, swap_seg] # destroy_repair, swap_seg
   
    n_op = length(ls_name)
    idx_op =collect(1:n_op)
    step_size = 100; total_cost=-1; pre_best_cost=100000
  
    count_stagnant = 0
    #rand weight
    weights_rand = ones(Float32, length(ls_name)) .* 1/length(ls_name)
    weights = Weights(weights_rand)
    for idx in 1:N_ITER
        if idx % step_size == 0 
            if  abs(pre_best_cost - sol_best.total_cost_with_penalty) < 0.0005 * pre_best_cost 
                count_stagnant += 1
            else
                count_stagnant = 0
            end
            pre_best_cost = sol_best.total_cost_with_penalty
        end
        
        if count_stagnant ==  max_stagnant_multiplier 
            return sol_best
        end
        
        i_imp += 1
        # Random selecting LS operator
        ls_idx = sample(idx_op, weights)
        ls = ls_name[ls_idx]    
        #  @show( ls)
        LS(ls, solution, soluT, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, T, Li)
        cost_solution(solution, global_parameter, darp, parameter_energy, fast_chg, instance)
        
        if idx % step_size_remove_route == 0  
            remove_route(solution, soluT, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
            cost_solution(solution, global_parameter, darp, parameter_energy, fast_chg, instance)
        end
        if idx % step_size_create_route == 0  
           create_route(solution, soluT, recorder_routes, recorder_lscost_routes, relatedness, demand_data, fast_chg, global_parameter, instance, darp, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, Q, T, Li)
           cost_solution(solution, global_parameter, darp, parameter_energy, fast_chg, instance)
        end

        if (solution.total_cost_with_penalty < sol.total_cost_with_penalty + T && check_chg_occ_constr(solution, occ_state_chg_ri, occ_state_chg_rj, discrte_t_interval)) || (solution.n_route < sol.n_route)               
            
            update_sol(sol, solution)
            copy_route_info(copy_e_r, copy_cap_r, copy_l_r, copy_maxcap_r, copy_dt_r, copy_dt_r_back, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
            if isempty(sol.unserved_users) && (sol.total_cost_with_penalty < sol_best.total_cost_with_penalty || (sol.n_route < sol_best.n_route && (sol.co2_emission <= co2_threshold))) 
             
                veh_exchange(sol, darp, fast_chg, parameter_energy, instance, global_parameter)
                update_sol(sol_best, sol) 
                copy_route_info(best_e_r, best_cap_r, best_l_r, best_maxcap_r, best_dt_r, best_dt_r_back, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back) 
                i_imp = 0
            end             
        else
            update_sol(solution, sol) 
            copy_route_info(e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, copy_e_r, copy_cap_r, copy_l_r, copy_maxcap_r, copy_dt_r, copy_dt_r_back)
        end
        
        if i_imp > 0            
            T -= T_max / T_red
            if T<0
                T = rand()*T_max
                if i_imp > n_imp * sol_best.n_route
                    update_sol(sol, sol_best)
                    update_sol(solution, sol_best)
                    copy_route_info(e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, best_e_r, best_cap_r, best_l_r, best_maxcap_r, best_dt_r, best_dt_r_back)
                    i_imp = 0
                end
            end            
        end
    end    
    return sol_best
end




function set_co2_threshold(id_dataset, n_cus, co2_reduc_target, is_gasoline_fleet, n_run, global_parameter) 

    id_dataset ∉ collect(1:3) && @warn("id_dataset ∉ collect(1:3) !!",id_dataset )
    vec_co2_max_milp = zeros(Float32,14) # c10,c20,...,c60,70,80,90,100,200,300,400,500 
    if id_dataset == 1
        vec_co2_max_milp[1:10]=[26.3277762667357, 39.139795180004, 58.2477553979516 , 72.9762773663471, 85.571741114915, # max co2 emission for c10 to c50 instance based on the milp solutions, used for testing the darp solver with co2 constraints
                                 102.2032318, 123.4625702, 137.2053223, 151.567337, 162.2761993]
        vec_co2_max_milp[6:14]= [104.209, 127.889, 142.713, 156.182 , 165.025,  293.256, 449.850, 525.444, 679.207] # obtained from SA using all gv
     
    elseif id_dataset == 2
        vec_co2_max_milp[1:10]=[21.89301682, 36.55890274, 62.94337082, 71.57552338,	94.22901917,
         101.9602051, 120.0026627, 132.0667877, 157.1488953, 173.6827087] # max co2 emission for c10 to c50 instance based on the milp solutions, used for testing the darp solver with co2 constraints                

    elseif id_dataset == 3      
        vec_co2_max_milp[1:10]=[22.68206024, 41.37122345, 54.78794479, 79.08572388, 94.9900589, 111.3295822, 
                              117.338501, 119.8001556, 164.4554901, 166.1737366]
        # max co2 emission for c10 to c50 instance based on the milp solutions, used for testing the darp solver with co2 constraints
    end

    idx=0
    if  n_cus > 100
        idx= 9 + floor(Int, n_cus/100)
    else
        idx= floor(Int, n_cus/10)
    end
    idx > length(vec_co2_max_milp) && error("idx > length(vec_co2_max_milp) error!!, instance and vec_co2_max_milp not match!!")
    
    global_parameter.is_gasoline_fleet  = is_gasoline_fleet
    if !is_gasoline_fleet
        global_parameter.co2_gasoline_fleet = vec_co2_max_milp[idx]
        global_parameter.co2_threshold_init = vec_co2_max_milp[idx]
        global_parameter.co2_reduc_target   = co2_reduc_target # % of co2 reduction, user-defined parameter 
        global_parameter.co2_threshold      = global_parameter.co2_threshold_init * (1.0 - co2_reduc_target)
    end
    global_parameter.n_run              = n_run
 
end


function eight_step_test(route::Vector{Int32}, ri, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)  
                                   
    N_NODE = 2*n_cus+2
    N=length(route)
    t_start = ei[route[1]]
    A, B, W, D, Load, RT = calculate_ads(route, t_start, ei, li, s_t, dist, n_cus, qi)
    if ! check_tw_cap(route, B, li, Load, ri, Q) 
         @show("check_tw_cap(route, B ,Load) failed!! ")
         return false,-1,-1,-1,-1,-1,-1
    end
    F0 = calc_Fi(1, route, W, B, RT, Li, li)
    Wp = sum(W[2:N-1])
    t_start_new = ei[route[1]] + min(F0, Wp)
    A, B, W, D, Load, RT = calculate_ads(route, t_start_new, ei, li, s_t, dist, n_cus, qi)
    if check_ridetime(route, RT, Li) && check_routeduration(route, B, TH) 
        return true, A, B, W, D, Load, RT
    end
    for i in 2:N-1
        if  1 < route[i] <= n_cus+1 
            Fi = calc_Fi(i, route, W, B, RT, Li, li)
            Wp = sum(W[i+1:N-1])              
            W[i] += min(Fi, Wp)
            B[i] = A[i] + W[i]
            D[i] = B[i] + s_t[route[i]]
            for j in i+1:N
                A[j] = D[j-1] + dist[route[j-1], route[j]]
                B[j] = max(A[j], ei[route[j]])
                W[j] = B[j] - A[j]
                D[j] = B[j] + s_t[route[j]]
                if  n_cus+1 < route[j] < N_NODE
                    idx2 = findfirst(x->x==route[j]-n_cus, route)
                    # RT[j] = A[j] - D[idx2]
                    RT[j] = B[j] - D[idx2]
                end
            end
        end
    end
    #check whether the drive time is respected for all dropoff vertices after i, if so the route is feasible 
    if check_routeduration(route, B, TH) && check_ridetime(route, RT, Li)
        return true, A, B, W, D, Load, RT
    else           
        @show(check_routeduration(route, B, TH), check_ridetime(route, RT, Li))
        return false, -1,-1,-1,-1,-1,-1
    end
end




# insert charging events on a route. First, identify possible charging positions of the route, then insert randomly the charging event 
# on a random position and a random (random policy)/best(greedy policy) charger 
# function insert_charging_route(solution::Solution, darp, r, route::Vector{Int32}, dist_all, ec_route, fast_chg::Fast_chargers, parameter_energy, dt_r, dt_r_back, flag_init_sol)
        
#     degree_rand_chg_policy = parameter_energy.degree_rand_chg_policy
#     dist_orig, n_cus  = darp.dist_orig, darp.n_cus 
#     ei, li, s_t, dist, qi = darp.ei, darp.li, darp.s_t, darp.dist, darp.qi
#     set_physical = darp.set_physical 

#     v=Int32[]
#     max_num_chg = parameter_energy.max_num_chg # max num of charging events limits for a route 
#     success_chg= false; remain_e_chg =0; idx_loc =0
#     set_fast_chg = fast_chg.set_id_chg_installed
#     pos_insert_chg = zeros(Int32, 60, 3) # store the pair nodes (precedent and successive nodes) of the location of the charger
#     n_pos_insert_chg = 1; t_v_v1 = 0; t_arr_node_after_chg = 0; idx_pre_chg_pos =0
#     info_chg = zeros(Float32, 2+max_num_chg*4) # store at most 4 charging operations  
#     # find the positions (nodes) where the SOC (state of charge) of the vehicle is insufficient  
#     pos_insert_chg[n_pos_insert_chg,:] = [route[1], route[2], 1]
#     veh_id = solution.RI[r,5]
#     max_ec = parameter_energy.max_ec[veh_id]
#     i_ec_violate= 0 # will alwayse find a positive position 
#     dt_r_eng = dt_r[route] .* parameter_energy.beta[veh_id]
#     i_ec_violate = findfirst(x->x-max_ec>0, dt_r_eng)
#     # @show(r, route, dt_r[route], dt_r_eng, max_ec,i_ec_violate )
#     # vehicle go to recharge when there are no passengers onboard
#     cap_r_tmp = zeros(Int32, i_ec_violate-1)  
#     for (idx, vi) in enumerate(route[1:i_ec_violate-2])
#         cap_r_tmp[idx+1] =  cap_r_tmp[idx] + darp.qi[route[idx+1]]
#     end
#     idx_ri   = findall(x->x==0, cap_r_tmp) # find location of zero load    
#     # @show(route, i_ec_violate, cap_r_tmp,darp.qi[route],idx_ri )
#     for idx in idx_ri
#         if idx > 1
#             n_pos_insert_chg += 1
#             pos_insert_chg[n_pos_insert_chg, :] = [route[idx], route[idx+1], idx]
#         end
#     end
#     # @show(route, pos_insert_chg)
#     #randomly select a fast charger
#     idx_first_loc =  rand(collect(1:n_pos_insert_chg))
#     # chg is the fast charger index in the list of all fast chargers (vec_idx_chg)
#     v= pos_insert_chg[idx_first_loc,1:2]
#     # rand() > 0.5 ? rand_policy =true : rand_policy = false
#     # flag_init_sol && (rand_policy = true)
#     rand() > degree_rand_chg_policy ? rand_policy =true : rand_policy = false 

#     if rand_policy
#         idx_chg = rand(set_fast_chg)
#     else
#         idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg::Fast_chargers, veh_id, parameter_energy, ec_route)
#     end
#     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)

#     if (dt_r[v[1]] + t_access_1) * parameter_energy.beta[veh_id] > parameter_energy.max_ec[veh_id] && idx_first_loc > 1
#         idx_first_loc -= 1
#         v= pos_insert_chg[idx_first_loc,1:2]
#         t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
#     end

#     vec_delta_e  =  parameter_energy.E_min[veh_id] - (parameter_energy.E_init[veh_id] - (ec_route + parameter_energy.beta[veh_id] * (t_access_1+ t_access_2 - t_v_v1)))
#     # vec_delta_e <= 0 && error("vec_delta_e <=0 error!!")
#     t_chg = vec_delta_e / fast_chg.pw[idx_chg]
#     idxs_chg_loc = pos_insert_chg[1:n_pos_insert_chg, 3]
#     # @show(idx_loc, v, t_chg, idxs_chg_loc)
#     forward_time_slack, vec_t_start = get_forward_time_slack(route, ei, li, s_t, dist, n_cus, qi, idxs_chg_loc, dist_orig, darp)
#     # @show(forward_time_slack)
#     count_chg = 1; next_recharge = true; 
#     n_candicate_loc =  n_pos_insert_chg - idx_first_loc +1
#     # @show(pos_insert_chg, idxs_chg_loc)
#     remain_e_chg = vec_delta_e
#     if t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[idx_first_loc]
#         # @show("t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[1] !")
#         additional_time = t_access_1 + t_access_2 - t_v_v1 +  t_chg
#         # t_start_chg = vec_t_start[1] + t_access_1 + rand(1)[1]* (forward_time_slack[idx_first_loc]-t_chg) # random start time
#         delta_t = rand() *(forward_time_slack[idx_first_loc] - ( t_chg + t_access_1 + t_access_2 - t_v_v1))
#         t_start_chg = vec_t_start[idx_first_loc] + t_access_1 +  delta_t # random start time
#         info_chg[1] += additional_time;  info_chg[2] +=1
#         info_chg[3:6] = [v[1], idx_chg, t_start_chg, t_start_chg+t_chg]
#         success_chg = true
#         # @show( t_chg + t_access_1+ t_access_2 - t_v_v1 , forward_time_slack[1], info_chg[1:6])
#     else  
#         while next_recharge  
#             if count_chg >= n_candicate_loc  || count_chg >= max_num_chg
#                 success_chg = false
#                 break
#             else
#                 # @show(count_chg,v[1])
#                 next_recharge, remain_e_chg, t_arr_node_after_chg =  schedule_chg(info_chg, darp, parameter_energy, remain_e_chg, r, route, veh_id, forward_time_slack, vec_t_start, idx_pre_chg_pos, t_arr_node_after_chg, fast_chg, pos_insert_chg, idx_first_loc, count_chg, v[1], idx_chg, t_access_1, t_access_2, t_v_v1)
#                 idx_loc_pre = idx_first_loc - 1 + count_chg
#                 idx_pre_chg_pos = pos_insert_chg[idx_loc_pre, 3]
#                 # next location
#                 count_chg += 1
#                 if count_chg < n_candicate_loc +1
#                     v = pos_insert_chg[idx_first_loc-1+count_chg,1:2]
#                     if rand_policy    
#                         idx_chg = rand(set_fast_chg)
#                     else
#                         idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg::Fast_chargers, veh_id, parameter_energy, ec_route)
#                     end
#                     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
#                 end
#             end
#         end
#         # @show( "! next_recharge && success_chg = true : ")
#         ! next_recharge && (success_chg = true)
#     end
#     if success_chg              
#         solution.vec_chg_event[r,:] = info_chg
#         @show(success_chg, r, route, info_chg ) 
#         check_energy(solution, darp, route, r, parameter_energy, fast_chg)
#         return true,  pos_insert_chg, forward_time_slack, vec_t_start
#     else     
#         return false, pos_insert_chg, forward_time_slack, vec_t_start
#     end
# end


# # repair the charging operations by inserting charging operations starting from the first location after the depot
# function repair_charging(solution::Solution, darp, r, route::Vector{Int32}, dist_all, fast_chg::Fast_chargers, parameter_energy, set_physical, dt_r,  pos_insert_chg, forward_time_slack, vec_t_start)
    
#     degree_rand_chg_policy = parameter_energy.degree_rand_chg_policy
#     dist_orig, end_depot = darp.dist_orig, darp.end_depot 
#     v=Int32[]
#     success_chg= false
#     max_num_chg = parameter_energy.max_num_chg # max num of charging events on a route
#     set_fast_chg = fast_chg.set_id_chg_installed
#     n_pos_insert_chg = length(forward_time_slack); t_v_v1 = 0; t_arr_node_after_chg =0; idx_pre_chg_pos = 0
#     info_chg = zeros(Float32, 2+max_num_chg*4) # store at most 4 charging operations  
#     veh_id = solution.RI[r,5]
#     ec_route =  dt_r[end_depot] * parameter_energy.beta[veh_id]
#     compute_dt_r(route, dt_r, dist_orig)       
#     rand() > degree_rand_chg_policy ? rand_policy =true : rand_policy = false
#     v= pos_insert_chg[1,1:2]; idx_loc =1
#     if rand_policy    
#         idx_chg = rand(set_fast_chg)
#     else
#         idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg::Fast_chargers, veh_id, parameter_energy, ec_route)
#     end
#     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
    
#     vec_delta_e  =  parameter_energy.E_min[veh_id] - (parameter_energy.E_init[veh_id] - (ec_route + parameter_energy.beta[veh_id] * (t_access_1+ t_access_2 - t_v_v1)))
#     # vec_delta_e <= 0 && error("vec_delta_e <=0 error!!")
#     t_chg = vec_delta_e / fast_chg.pw[idx_chg]

#     count_chg = 1; next_recharge = true; idx_first_loc =2
#     # @show(pos_insert_chg, idxs_chg_loc)
#     remain_e_chg = vec_delta_e
#     if t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[1]
#         additional_time = t_access_1 + t_access_2 - t_v_v1 +  t_chg
#         # t_start_chg = vec_t_start[1] + t_access_1 + rand(1)[1]* (forward_time_slack[1]-t_chg) # random start time
#         delta_t = rand() *(forward_time_slack[1] - ( t_chg + t_access_1 + t_access_2 - t_v_v1))
#         t_start_chg = vec_t_start[1] + t_access_1 + delta_t # random start time
#         info_chg[1] += additional_time;  info_chg[2] +=1
#         info_chg[3:6] = [v[1], idx_chg, t_start_chg, t_start_chg+t_chg]
#         success_chg = true
#     else        
#         while next_recharge  
#             if count_chg >= n_pos_insert_chg || count_chg >= max_num_chg
#                 success_chg = false
#                 break
#             else
#                 next_recharge, remain_e_chg, t_arr_node_after_chg =  schedule_chg(info_chg, darp, parameter_energy, remain_e_chg, r, route, veh_id, forward_time_slack, vec_t_start, idx_pre_chg_pos, t_arr_node_after_chg, fast_chg, pos_insert_chg, idx_first_loc, count_chg, v[1], idx_chg, t_access_1, t_access_2, t_v_v1)
#                 idx_loc_pre = idx_first_loc - 1 + count_chg
#                 idx_pre_chg_pos = pos_insert_chg[idx_loc_pre, 3]
#                 count_chg += 1
#                 if count_chg < n_pos_insert_chg +1
#                     v = pos_insert_chg[count_chg,1:2]
#                     if rand_policy    
#                         idx_chg = rand(set_fast_chg)
#                     else
#                         idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg::Fast_chargers, veh_id, parameter_energy, ec_route)
#                     end
#                     # @show(count_chg, v, pos_insert_chg)
#                     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
#                 end
#             end
#         end
#         ! next_recharge && (success_chg = true)
#     end
#     if success_chg              
#         solution.vec_chg_event[r,:] = info_chg
#         # @show(" test 24: ", success_chg, r, route, info_chg ) 
#         check_energy(solution, darp, route, r, parameter_energy, fast_chg)   
        
#     end
#     return success_chg
# end


# function  check_energy_simple(solution::Solution, darp, route, r, parameter_energy, fast_chg )
    
#     veh_id = solution.RI[r,5]
#     if parameter_energy.is_electric[veh_id]
#         e_consum = parameter_energy.beta[veh_id] * (length_route_energy_check(route, darp)) 
#         _ , eng_chged = get_route_chg_time_energy(solution, r, fast_chg)
#         # @show(r,veh_id, parameter_energy.is_electric[veh_id], e_consum, eng_chged, parameter_energy.max_ec[veh_id])
#         if round(e_consum - eng_chged - parameter_energy.max_ec[veh_id], digits = DIGITS) > 0
#             @show(r, route, veh_id, e_consum, eng_chged, parameter_energy.max_ec[veh_id], e_consum - eng_chged - parameter_energy.max_ec[veh_id], solution.vec_chg_event )
#             error("check energy failed !")
#             return false
#         else
#             # @show("check_energy successful !" )
#             return true 
#         end
#     else
#         return true
#     end
# end
  
# function verify_solution(solution::Solution, Q)

     
#     K, T_max, n_charger = instance.K, instance.T_max, instance.n_charger
#     n_cus, start_depot, end_depot = darp.n_cus, darp.start_depot, darp.end_depot
#     parameter_energy = instance.parameter_energy
#     # set_physical, dist_all = darp.set_physical, darp.dist_all
#     ei, li, qi, Li, s_t, dist_orig, TH = darp.ei, darp.li, darp.qi, darp.Li, darp.s_t, darp.dist_orig, darp.TH

#     count_sol= zeros(Int32,2*n_cus+1)
#     infeasible_precedence = 0; infeasible_tw_cap_ride_time = 0; infeasible_energy=0; infeasible_double_tour = 0
#     n_route= solution.n_route
 
#     for _r in 1:n_route
#         route = get_route(solution, _r, start_depot, end_depot) 
#         infeasible_energy += !check_energy(solution, route, _r, parameter_energy, fast_chg)
#         infeasible_precedence += !check_precedence(route)
#         infeasible_tw_cap_ride_time += !eight_step(route, _r, ei, li, s_t, dist_orig, Q, Li, n_cus, TH)
#         for i in 2:length(route)-1
#             count_sol[route[i]] +=1 
#         end
#         nodes_route = nodes[route[2:end-1]]
#         infeasible_double_tour += !no_double_tour(nodes_route, D′) 
#     end
#     unserved_user = collect(solution.unserved_users)
#     for user in unserved_user
#         count_sol[user] +=1 ;  count_sol[user+n_cus] +=1
#     end
#     # check used fleet
#     fleet_error = !check_used_vehs(solution)

#     if sum(count_sol) != 2*n_cus || infeasible_precedence > 0 || infeasible_tw_cap_ride_time > 0 || fleet_error || infeasible_energy > 0 || infeasible_double_tour > 0
#         @show("verify_solution failed: ",count_sol, infeasible_precedence, infeasible_tw_cap_ride_time, fleet_error, infeasible_energy, infeasible_double_tour)
#         if infeasible_precedence > 0
#             @show("check precedence failed !!", infeasible_precedence)
#             for r in 1:n_route
#                 route = get_route(solution, r, start_depot, end_depot) 
#                 @show(r, route)
#             end
#         end
#         if infeasible_tw_cap_ride_time > 0
#             for _r in 1:solution.n_route 
#                 if ! verify_eight_step(solution, _r, Q, Li)
#                     route = get_route(solution, _r, start_depot, end_depot)
#                     A, B, W, D, Load, RT = calculate_ads(route, 0, ei, li, s_t, dist_orig, n_cus, qi)
#                     @show("tw check failed : ",_r,  route ) 
#                     @show(ei[route], li[route], A, li[route]-B, W, D, Load, RT, Li[route])
#                 end 
#             end
#         end
#         fleet_error && @show("used vehicle error! ", solution.RI[1:n_route,5])
#         infeasible_energy > 0 && @show("infeasible_energy error : ")
#         infeasible_double_tour > 0 && @show("infeasible_double_tour > 0: ", infeasible_double_tour)
#         error("verify_solution failed")
#     else
#         @show("valid solution")
#         # @show(count_sol, infeasible_precedence, infeasible_tw_cap_ride_time, fleet error, infeasible_energy)
#     end
# end

# #schedule a charging event (start time, end time, duration and the energy to be recharged at next charging location of the route)
# function schedule_chg(n_pos_insert_chg, info_chg, darp, parameter_energy, remain_e_chg, r, route, veh_id, forward_time_slack, vec_t_start, idx_pre_chg_pos, t_arr_node_after_chg, fast_chg::Fast_chargers, pos_insert_chg, idx_first_loc, count_chg, v1, idx_chg, t_access_1, t_access_2, t_v_v1)
    
#     li = darp.li
#     add_access_chg_time = t_access_1+ t_access_2 - t_v_v1
#     if count_chg > 1
#         remain_e_chg += (add_access_chg_time * parameter_energy.beta[veh_id])
#     end
#     idx_loc = idx_first_loc - 1 + count_chg
#     idx_pos = pos_insert_chg[idx_loc,3]
#     if count_chg > 1 # update forward_time_slack calculation
#         W, B = calculate_ads_chg_schedule(route[idx_pre_chg_pos+1:end], t_arr_node_after_chg, darp)  # re_calculate A,B,W,D after the precedent chg event
#         idx_Fi = idx_pos - idx_pre_chg_pos
#         if idx_Fi < 0
#             @show(idx_first_loc, count_chg, idx_loc, n_pos_insert_chg, n_candicate_loc =  n_pos_insert_chg - idx_first_loc +1)
#             @show(idx_loc, pos_insert_chg[idx_loc, :], pos_insert_chg,  pos_insert_chg[idx_loc,3], idx_pos, idx_pre_chg_pos)
#         end
#         forward_time_slack_new, t_start_fw_slack = get_forward_time_slack_revise(route[idx_pre_chg_pos+1:end], idx_Fi, W, B, li)
#         # t_chg_constr = forward_time_slack_new - add_access_chg_time
#         t_chg_constr = max(0, forward_time_slack_new - add_access_chg_time)
#     else
#         # t_chg_constr = forward_time_slack[idx_loc] - add_access_chg_time
#         t_chg_constr =  max(0,forward_time_slack[idx_loc] - add_access_chg_time,0)
#         t_start_fw_slack = vec_t_start[idx_loc] 
#     end
    
#     t_chg_desired = remain_e_chg / fast_chg.pw[idx_chg]
#     t_chg= min(t_chg_constr, t_chg_desired)
#     additional_time = add_access_chg_time +t_chg 
#     info_chg[1] += additional_time;  info_chg[2] +=1
#     idxs = collect(2+(count_chg-1)*4+1 : 2+count_chg*4)
#     # @show(forward_time_slack[count_chg] , add_access_chg_time,t_chg_constr,t_chg_desired)
#     # @show(count_chg, t_chg_constr, t_chg_desired, t_chg )
#     if  t_chg_desired > t_chg_constr
#         t_chg_start = t_start_fw_slack + t_access_1
#         t_chg_end = t_chg_start + t_chg
#         info_chg[idxs] = [v1, idx_chg, t_chg_start, t_chg_end]
#         remain_e_chg -= (t_chg * fast_chg.pw[idx_chg])
#         t_arr_node_after_chg = t_chg_end + t_access_2
#         return true, remain_e_chg, t_arr_node_after_chg
#     else
#         # @show(remain_e_chg=0)
#         t_chg_start = t_start_fw_slack + t_access_1 + rand(1)[1]*(t_chg_constr - t_chg) # random start time
#         t_chg_end = t_chg_start + t_chg
#         info_chg[idxs] = [v1, idx_chg, t_chg_start, t_chg_end]
#         return false, 0, 0
#     end 
# end

# function gain_change_route_type(solution, route, r, dist_route, new_type_veh, darp, parameter_energy)

#     current_cost = cost_route(solution, r, start_depot, darp, parameter_energy)
#     cap_new_type = parameter_energy.veh_info.cap_veh
#     check_capacty(route,cap_new_type)

# end
 
# # swap electric to gasoline if have cost saving (based on the input setting)
# function veh_exchange(solution::Solution, darp, fast_chg, parameter_energy, instance, global_parameter, dt_r, dt_r_back)
    
#     start_depot = darp.start_depot; k_max = instance.K_MAX
#     # @show(k_max,  instance.K_MAX)
#     vec_n_veh = instance.vec_n_veh # fleet_size of each veh type
#     n_route = solution.n_route
#     co2_threshold = global_parameter.co2_threshold
#     penalty = global_parameter.penalty 
#     for r in 1:n_route 
#         veh_id= solution.RI[r,5]  
#         if parameter_energy.is_electric[veh_id]       
#             # @show(length_route(solution, r, start_depot, darp), global_parameter.threshold_advantage_ev)
#             veh_type = parameter_energy.veh_type[veh_id]  
#             if length_route(solution, r, start_depot, darp) > global_parameter.threshold_advantage_ev[veh_type]
#                 # find veh_id that is a gasoline veh and not on the current used veh list
#                 co2= solution.co2_emission + co2_route_gv(solution, r, start_depot, darp, parameter_energy)
#                 # @show(co2, co2_threshold, solution.co2_emission , co2_route_gv(solution, r, start_depot, darp, parameter_energy))
#                 if co2 <= co2_threshold
#                     # @show("co2 <= co2_threshold !!")
#                     veh_ids_sol= solution.RI[1:n_route,5]
#                     veh_ids_candidate = collect(1:k_max) 
#                     setdiff!(veh_ids_candidate, veh_ids_sol) 
#                     idx_tmp = findfirst(x->x==false, parameter_energy.is_electric[veh_ids_candidate])
#                     if isnothing(idx_tmp) 
#                         error("not sufficient size of gasoline fleet, please increase K_MAX !!") 
#                     end
#                     veh_id_new = veh_ids_candidate[idx_tmp]  # list gv
#                     # @show(solution.RI[r,5],veh_id_new, parameter_energy.is_electric[veh_id_new])
#                     solution.RI[r,5] =  veh_id_new
#                     solution.total_cost, solution.total_cost_with_penalty = cost_solution(solution, global_parameter, darp, parameter_energy, fast_chg, instance, dt_r, dt_r_back)
#                     solution.vec_chg_event[r,:] .*= 0
#                 end
#             end
#         end 
#     end   
# end


# function show_eight_step_test(solution::Solution, r)

#         route = get_route(solution, r, 1, 2*n_cus+2)
#         success, A, B, W, D, Load, RT = eight_step_test(route, r, ei, li, s_t, dist, Q, Li, n_cus, TH)
#         @show("show_eight_step_test : ", r, route,  success)
#         tw_detail = DataFrame(node = route, A=A, B=B, W=W, D=D, RT=RT, Load=Load, ei = ei[route], li = li[route], st = s_t[route], Li=Li[route])
#         XLSX.writetable("tw_ride_time_detail_orig_r" * string(r) * ".xlsx", overwrite=true, tw_detail)
# end

# function SP(soluT::Solution, sol_best::Solution, sol_ref::Solution, pool_routes, lookup_route, costs_route, Ri, Rt, Nt, e_r, cap_r, l_r, maxcap_r)

#     # @show(length(pool_routes))
#     n_route = length(pool_routes)
#     nt      = length(Nt)
#     # @show(n_route, nt,costs, length(pool_routes))
#     costs= copy(costs_route)
#     model = direct_model(Gurobi.Optimizer(GRB_ENV))
#     # set_optimizer_attribute(model, "Presolve", 0)
#     set_optimizer_attribute(model, "OutputFlag", 0)
#     set_silent(model);set_route=Set()
#     # setup warmstart
#     y_start=zeros(Bool, n_route)
#     has_key = true
#     for r in 1:sol_best.n_route
#         route = get_route(sol_best,r)
#         cap_veh = sol_best.RI[r,4] + 10000
#         set_route = union(Set(route), cap_veh)
#         if haskey(pool_routes, set_route)
#             y_start[pool_routes[set_route][1]]=1
#         else
#             has_key = false
#             break
#         end
#     end
#     if ! has_key
#         y_start = y_start .* 0
#         for r in 1:sol_ref.n_route
#             route = get_route(sol_ref,r)
#             cap_veh = sol_ref.RI[r,4] + 10000
#             set_route = union(Set(route), cap_veh)
#             y_start[pool_routes[set_route][1]]=1
#         end
#     end
    
#     @variable(model, y[1:n_route], Bin)
#     set_start_value.(y, y_start) # warm start
#     @objective(model, Min, (costs' * y))
#     @constraint(model, [i = 1:n_cus], sum(y[j] for j in Ri[i]) == 1)
#     @constraint(model, [t = 1:nt], sum(y[j] for j in Rt[t]) <= Nt[t])
#     optimize!(model)
#     # @show(lookup_route)
#     if termination_status(model) == MOI.OPTIMAL
#         y= round.(Int, value.(y))
#         idxs= findall(x ->x>0, y)
#         # @show(y,idxs)
#         soluT.n_route=length(idxs)
#         soluT.total_cost =  objective_value(model)
#         for r in 1:length(idxs)
#             # @show(idxs[r])
#             route = lookup_route[idxs[r]]
#             # @show(route)
#             route_auxiliary_update(soluT, route, r, e_r, cap_r, l_r, maxcap_r)
#         end

#         # if abs(objective_value(model) - length_solution(soluT)) > 0.001 
#         #     @show(objective_value(model), length_solution(soluT))
#         #     error(" objective_value(model) ! = length_solution(soluT) error !!")
#         # end
#         # verify_solution(soluT, Q)
         
#         return  objective_value(model)
#     else
#         return  -1
#     end
# end

# function record_chg_event(solution::Solution, r, v1, fast_chg_id, t_start, t_end, additional_time)
#     solution.vec_chg_event[r, 1] += additional_time
#     solution.vec_chg_event[r, 2] += 1
#     idx_chg_event =  solution.vec_chg_event[r, 2] 
#     idxs = 2 + (idx_chg_event-1) * 4 + 1: 2 + idx_chg_event * 4
#     solution.vec_chg_event[r,idxs] = [v1, fast_chg_id, t_start, t_end]
# end


# #########################################"
# # Breakers et al. (2014)'s SA algorithm, not used for flexbus problem
# #######################################"
# function SA(solution::Solution, soluT::Solution, solu_init::Solution, homogenous_fleet, Nt, type_veh, e_r, cap_r, l_r, 
#     maxcap_r, N_ITER, avg_dist_between_nodes, Q ,vec_total_cost, start_depot, end_depot, qi, Li, layer_nodes, lyrs_compatible)
    
#     # solution = solu
#     sol_best=deepcopy(solution)  # x_best
#     sol    = deepcopy(solution)  # x, x' => solution 
    
#     copy_e_r   = copy(e_r);     copy_cap_r = copy(cap_r)
#     copy_l_r   = copy(l_r);     copy_maxcap_r = copy(maxcap_r)
#     copy_dt_r  = copy(dt_r);    copy_dt_r_back = copy(dt_r_back)
#     best_e_r   = copy(e_r);     best_cap_r = copy(cap_r)
#     best_l_r   = copy(l_r);     best_maxcap_r = copy(maxcap_r)
#     best_dt_r  = copy(dt_r);    best_dt_r_back = copy(dt_r_back)
    
#     # SA parameter
#     T_red = 300
#     t_max = 1.2 
    
#     T_max = t_max *avg_dist_between_nodes
#     n_imp = 400
#     T = T_max
    
#     i_imp=0 ; 
#     if homogenous_fleet
#         ls_name= [relocate, two_opt_star, two_opt, swap_two_users, four_opt_intra, swap_seg, relocate_worst, relocate_rand, create_route]
#         # ls_name= [relocate, two_opt_star, two_opt, swap_two_users, four_opt_intra, remove_route]
#     else
#         ls_name= [swap_two_vehs, relocate, two_opt_star, two_opt, swap_two_users, four_opt_intra, swap_seg, relocate_worst, relocate_rand, create_route]
#     end
    
#     rho = 0.5
#     n_op = length(ls_name)
#     # delta=[3,2,1] # new best , improve,  not improve 
#     delta=[9,5,1] #  new best , improve,  not improve
#     score = ones(Float32,n_op)*delta[2]
#     score_pre = copy(score)
#     count_score = zeros(Int32,n_op)
#     weights = Weights(score ./ sum(score)) #init score
#     idx_op =collect(1:n_op)
#     step_size = 100
#     #  step_size_Gurobi=100000 # SP procedure not used
    
#     # set partitioning problem setup
#     # pool_routes  = Dict()
#     # lookup_route = Dict()
#     # dict_pool_routes = Dict()
#     # costs_route  = Float32[]
#     # lookup_typeveh = Int32[]
#     # key_elite_routes = Set()
#     # Ri = Vector{Set{Int32}}(undef, n_cus)
#     # Rt = Vector{Set{Int32}}(undef, length(Nt))
#     # for i in 1:n_cus
#     #     Ri[i] = Set()
#     # end
#     # for i in 1:length(Nt)
#     #     Rt[i] = Set()
#     # end
    
#     # total_cost=-1; key=[0]
#     # idx_CP =1; elite_sol = false ; idx_stagnant =0
#     # pool_pre = Any[]; pool_succ=Any[]; pool_RI = Any[]; pool_TT_walk =Any[] ;
#     # pool_unserved = Any[]; pool_penalty= Float32[]; pool_n_route = Int32[]; pool_total_cost = Float32[]
#     pre_best_cost=100000; top_elite = 0.5 ; max_stagnant =3

#     # add_sol_routes_pool(solution, pool_routes, dict_pool_routes, lookup_route, lookup_typeveh, 
#     #        type_veh, costs_route, key, start_depot, end_depot)    
#     # add_sol_pool(solution, pool_pre, pool_succ, pool_RI, pool_TT_walk, pool_n_route, pool_total_cost, pool_unserved, pool_penalty)
                       
#     count_stagnant = 0
#     total_cost_ref=solution.total_cost
#     #rand weight
#     weights_rand = ones(Float32, length(ls_name)) .* 1/length(ls_name)
#     weights = Weights(weights_rand)

#     for idx in 1:N_ITER
#         # if idx == N_ITER
#         #     @show( length(pool_routes) )
#         # end
        
#         if idx % step_size == 0 
#             # @show(count_score)
#             # weights, score, count_score, score_pre = update_weight(weights, score, count_score, score_pre, rho) 
#             # @show(weights, score,  score_pre)
#             if  abs(pre_best_cost - sol_best.total_cost) < 0.0005 * pre_best_cost 
#                 count_stagnant += 1
#             else
#                 count_stagnant = 0
#             end
#             pre_best_cost = sol_best.total_cost
#             push!(vec_total_cost, sol_best.total_cost) 
#             # @show(length(pool_routes))
#         end
     
#         if count_stagnant == 300
#             @show(idx, count_stagnant)
#             return sol_best              
#         end
 
#         i_imp += 1         
#         for i in 1:n_op
#             # Random selecting LS operator
#             ls_idx = sample(idx_op, weights)
#             ls = ls_name[ls_idx]
#             done = LS(ls, solution, soluT, e_r, cap_r, l_r, maxcap_r, Q, T, Li, layer_nodes, lyrs_compatible,, type_veh, cap_type_veh)
#             # done && (count_score[ls_idx] += 1)
#             # solution.total_cost = length_solution(solution)
#             solution.total_cost = length_solution(solution, λ, penalty, start_depot)
#             # elite_sol = false
#             # if (sol.n_route > K && solution.n_route < sol.n_route) || solution.total_cost < sol.total_cost + T
#             if solution.total_cost < sol.total_cost + T                
#                 # done && (score[ls_idx] += delta[2])
#                 update_sol(sol, solution)
#                 copy_route_info(copy_e_r, copy_cap_r, copy_l_r, copy_maxcap_r, copy_dt_r, copy_dt_r_back, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
#                 # sol.total_cost < sol_best.total_cost * 1.05 && (elite_sol = true)
#                 # add_sol_routes_pool(solution, pool_routes, dict_pool_routes, lookup_route, lookup_typeveh, 
#                 #      type_veh, costs_route, key, start_depot, end_depot)
#                 # add_sol_pool(solution, pool_pre, pool_succ, pool_RI, pool_TT_walk, pool_n_route, pool_total_cost, pool_unserved, pool_penalty)
               
#                 if (sol_best.n_route > K && sol.n_route <sol_best.n_route) || sol.total_cost < sol_best.total_cost 
#                     # done && (score[ls_idx] += (delta[1]-delta[2]))     
#                     update_sol(sol_best,sol) 
#                     copy_route_info(best_e_r, best_cap_r, best_l_r, best_maxcap_r, best_dt_r, best_dt_r_back, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back) 
#                     i_imp = 0
#                 end
#             else
#                 update_sol(solution,sol) 
#                 copy_route_info(e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, copy_e_r, copy_cap_r, copy_l_r, copy_maxcap_r, copy_dt_r, copy_dt_r_back)
#                 # done && (score[ls_idx] += delta[3])
#             end
#         end
#         if i_imp > 0            
#             T -= T_max / T_red
#             if T<0
#                 T = rand()*T_max
#                 if i_imp > n_imp * sol_best.n_route
#                     update_sol(sol,sol_best)
#                     update_sol(solution,sol_best)
#                     copy_route_info(e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, best_e_r, best_cap_r, best_l_r, best_maxcap_r, best_dt_r, best_dt_r_back)
#                     i_imp = 0
#                 end
#             end            
#         end
#     end
#     return sol_best
# end


# function SA_braekers(solution::Solution, soluT::Solution, homogenous_fleet, Nt, type_veh, e_r, cap_r, l_r, 
#            maxcap_r, N_ITER, avg_dist_between_nodes, Q , vec_total_cost, parameter_eng)
    
#     # solution.total_cost = length_solution(solution)
#     solution.total_cost = length_solution(solution, λ, penalty, start_depot)
#     sol_best=deepcopy(solution)  # x_best
#     sol    = deepcopy(solution)  # x, x' => solution 

#     copy_e_r   = copy(e_r);     copy_cap_r = copy(cap_r)
#     copy_l_r   = copy(l_r);     copy_maxcap_r = copy(maxcap_r)
#     best_e_r   = copy(e_r);     best_cap_r = copy(cap_r)
#     best_l_r   = copy(l_r);     best_maxcap_r = copy(maxcap_r)
#     # parameter
#     T_red = 300
#     t_max = 1.2 

#     T_max = t_max *avg_dist_between_nodes
#     n_imp = 400
#     T = T_max

#     i_imp=0 

#     if homogenous_fleet
#         ls_name= [relocate, two_opt_star, two_opt, four_opt_intra, remove_route]
#     else
#         ls_name= [ relocate, two_opt_star, two_opt, four_opt_intra, remove_route, swap_two_vehs]
#     end

#     order=collect(1:length(ls_name))
#     step_size = 100 
#     for idx in 1:N_ITER   
        
#         if idx % step_size == 0 
#             push!(vec_total_cost, sol_best.total_cost) 
#         end

#         i_imp += 1        
#         shuffle!(order)
#         for i in order
#             LS(ls_name[i], solution, soluT, e_r, cap_r, l_r, maxcap_r, Q,T, Li, layer_nodes, lyrs_compatible)
#             # solution.total_cost = length_solution(solution)
#             solution.total_cost = length_solution(solution, λ, penalty, start_depot)
#             if (sol.n_route > K && solution.n_route < sol.n_route) || solution.total_cost < sol.total_cost + T
#                 update_sol(sol, solution)
#                 copy_route_info(copy_e_r, copy_cap_r, copy_l_r, copy_maxcap_r, e_r, cap_r, l_r, maxcap_r)
#                 if (sol_best.n_route > K && sol.n_route <sol_best.n_route) || sol.total_cost < sol_best.total_cost                  
#                     update_sol(sol_best,sol)
#                     copy_route_info(best_e_r, best_cap_r, best_l_r, best_maxcap_r, e_r, cap_r, l_r, maxcap_r)
#                     i_imp = 0
#                 end
#             else
#                 update_sol(solution,sol) 
#                 copy_route_info(e_r, cap_r, l_r, maxcap_r, copy_e_r, copy_cap_r, copy_l_r, copy_maxcap_r)
#             end
#         end
#         if i_imp > 0            
#             T -= T_max / T_red
#             if T<0
#                 T = rand()*T_max
#                 if i_imp > n_imp * sol_best.n_route
#                     update_sol(sol,sol_best)
#                     update_sol(solution,sol_best)
#                     copy_route_info(e_r, cap_r, l_r, maxcap_r, best_e_r, best_cap_r, best_l_r, best_maxcap_r)
#                     i_imp = 0
#                 end
#             end            
#         end
#     end
#     return sol_best
# end

# SA hybrid with the integrated set partitioning procedure
# function SA_hybrid(solution::Solution, soluT::Solution, solu_init::Solution, bigM, homogenous_fleet, Nt, type_veh, e_r, cap_r, l_r, 
#     maxcap_r, dt_r, dt_r_back, N_ITER, max_stagnant_multiplier, avg_dist_between_nodes, Q ,vec_total_cost, start_depot, end_depot, qi, Li, layer_nodes, lyrs_compatible, flag_dev, nodes, D′, cap_type_veh, vec_n_veh, occ_state_chg_init, discrte_t_interval)
    
#     sol_best = deepcopy(solution)  # x_best
#     sol      = deepcopy(solution)  # x, x' => solution 
#     copy_e_r   = copy(e_r);     copy_cap_r = copy(cap_r)
#     copy_l_r   = copy(l_r);     copy_maxcap_r = copy(maxcap_r)
#     copy_dt_r  = copy(dt_r);    copy_dt_r_back = copy(dt_r_back)
#     best_e_r   = copy(e_r);     best_cap_r = copy(cap_r)
#     best_l_r   = copy(l_r);     best_maxcap_r = copy(maxcap_r)
#     best_dt_r = copy(dt_r);     best_dt_r_back = copy(dt_r_back)
    
#    =Int32[] # store capacity of veh info for the removed veh and use it to create new veh in local search
#     # SA parameter
#     T_red = 300
#     t_max = 1.2
    
#     T_max = t_max *avg_dist_between_nodes
#     n_imp = 400
#     T = T_max
    
#     i_imp=0 ; 
#     ls_name= [relocate, two_opt_star, two_opt, swap_two_users, swap_seg, four_opt_intra, relocate_worst, relocate_rand, create_route]
    
#     n_op = length(ls_name)
#     idx_op =collect(1:n_op)
#     step_size = 100; total_cost=-1; pre_best_cost=100000

#     ## Set Partitioning (SP), not used
#     # step_size_Gurobi=1000; max_num_route_with_chg = 10*1000 
#     # pool_routes  = Dict()
#     # pool_chg_event = Dict() # store the charging events info of routes
#     # lookup_route = Dict()
#     # dict_pool_routes = Dict()
#     # pool_veh_id = Dict() # veh_id of routes
#     # checked_table = Dict() # lookup table to see whether a pair of routes have been checked for chg occ conflict between them
#     # conflict_set = Set()    # set of routes having chg conflicts with other routes
#     # costs_route  = Float32[]
#     # lookup_typeveh = Int32[]
#     # # set partitioning problem setup
#     # Ri = Vector{Set{Int32}}(undef, n_cus)
#     # Rt = Vector{Set{Int32}}(undef, length(Nt))
#     # Rk = Vector{Set{Int32}}(undef, K) 
#     # R_conflict = Vector{Set{Int32}}(undef, max_num_route_with_chg)
#     # for i in 1:n_cus
#     #     Ri[i] = Set()
#     # end
#     # for i in 1:length(Nt)
#     #     Rt[i] = Set()
#     # end
#     # for i in 1:max_num_route_with_chg
#     #     R_conflict[i] = Set()
#     # end
#     # for i in 1:K
#     #     Rk[i] = Set()
#     # end
#     # key=[0]; idx_CP =1; idx_stagnant =0
#     # pool_pre = Any[]; pool_succ=Any[]; pool_RI = Any[]; pool_TT_walk =Any[] ;pool_sol_chg_event =Any[]
#     # pool_unserved = Any[]; pool_penalty= Float32[]; pool_n_route = Int32[]; pool_total_cost = Float32[]
#     # top_elite = 0.5 ; max_stagnant =3
#     # total_cost_ref=solution.total_cost
    
#     # add_sol_routes_pool(solution, pool_routes, pool_chg_event, pool_veh_id, dict_pool_routes, lookup_route, lookup_typeveh, 
#     # type_veh, costs_route, key, start_depot, end_depot)    
#     # add_sol_pool(solution, pool_pre, pool_succ, pool_RI, pool_TT_walk, pool_n_route, pool_total_cost, pool_unserved, pool_penalty, pool_sol_chg_event)
    
#     # @show(pool_routes, lookup_route, costs_route, lookup_typeveh, type_veh,key)
#     count_stagnant = 0
#     #rand weight
#     weights_rand = ones(Float32, length(ls_name)) .* 1/length(ls_name)
#     weights = Weights(weights_rand)
#     for idx in 1:N_ITER
#         # if idx == N_ITER && flag_dev == 1
#         #     @show( length(pool_routes) )
#         # end       
#         if idx % step_size == 0 
#             if  abs(pre_best_cost - sol_best.total_cost) < 0.0005 * pre_best_cost 
#                 count_stagnant += 1
#             else
#                 count_stagnant = 0
#             end
#             pre_best_cost = sol_best.total_cost
#         end
        
#         if count_stagnant ==  max_stagnant_multiplier
#             flag_dev ==1 && @show(idx, count_stagnant)
#             return sol_best
#         end
#         # idx % step_show == 0 && @show(sol_best.total_cost)
#         # if length(pool_routes) > step_size_Gurobi * idx_CP
#         #     # @show("iter = ",idx, " N_ITER = ", N_ITER)
#         #     idx_CP += 1
#         #     set_Ri_Rt(lookup_route, lookup_typeveh, Ri, Rt, Rk, pool_veh_id, R_conflict, key, pool_chg_event, occ_state_chg_init, discrte_t_interval, checked_table, conflict_set)
#         #     sol_ref = set_sol_ref(soluT, pool_pre, pool_succ, pool_RI, pool_n_route, pool_total_cost, pool_TT_walk, pool_unserved, pool_penalty, pool_sol_chg_event)
#         #     # val_tmp = sol_best.total_cost
#         #     soluT = deepcopy(solu_init)
#         #     total_cost = SP(soluT, sol_best, sol_ref, pool_routes, pool_chg_event, pool_veh_id, lookup_route, costs_route, Ri, Rt, Rk, R_conflict, conflict_set, Nt, e_r,cap_r,l_r,maxcap_r, dt_r, dt_r_back, lookup_typeveh, cap_type_veh)
#         #     if total_cost != -1 && total_cost <= total_cost_ref
#         #         update_sol(sol_best,soluT) 
#         #         update_sol(solution,soluT)
#         #         copy_route_info(best_e_r, best_cap_r, best_l_r, best_maxcap_r, best_dt_r, best_dt_r_back, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
#         #         if  abs(pre_best_cost - total_cost) < 0.005 * pre_best_cost 
#         #             idx_stagnant += 1
#         #         else 
#         #             idx_stagnant = 0
#         #         end
#         #         # @show(idx_stagnant)
#         #         pre_best_cost =   total_cost 
#         #     end                 
#         # end
        
#         # if idx_stagnant == max_stagnant
#         #     pool_pre, pool_succ, pool_RI, pool_unserved, pool_penalty, pool_TT_walk, pool_n_route, pool_total_cost, pool_sol_chg_event = restart(pool_pre, pool_succ, pool_sol_chg_event, pool_RI, pool_unserved, pool_penalty, pool_n_route, pool_total_cost, pool_TT_walk, top_elite)
#         #     reset_pool_routes(soluT, pool_routes, pool_chg_event, pool_sol_chg_event, dict_pool_routes, lookup_route, lookup_typeveh, 
#         #                type_veh, costs_route, key, Ri, Rt, Rk, R_conflict, Nt, pool_pre, pool_succ, pool_RI, pool_veh_id, pool_unserved, pool_penalty, pool_TT_walk, pool_n_route, pool_total_cost, start_depot, end_depot, checked_table, conflict_set)
#         #     pre_best_cost = solution.total_cost
#         #     idx_stagnant = 0; idx_CP =1
#         #     max_stagnant += 1
#         # end
        
#         i_imp += 1         
#         for i in 1:n_op
#             # Random selecting LS operator
#             ls_idx = sample(idx_op, weights)
#             ls = ls_name[ls_idx]        
#             LS(ls, solution, soluT, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, Q, T, Li, layer_nodes, lyrs_compatible,, type_veh, cap_type_veh)
#             verify_solution(solution, Q)
#             solution.total_cost = length_solution(solution, λ, penalty, start_depot)
#             if solution.total_cost < sol.total_cost + T && check_chg_occ_constr(solution, occ_state_chg_init, discrte_t_interval)              
#                 update_sol(sol, solution)
#                 copy_route_info(copy_e_r, copy_cap_r, copy_l_r, copy_maxcap_r, copy_dt_r, copy_dt_r_back, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back)
#                 ## SP
#                 # add_sol_routes_pool(solution, pool_routes, pool_chg_event, pool_veh_id, dict_pool_routes, lookup_route, lookup_typeveh, 
#                 #      type_veh, costs_route, key, start_depot, end_depot)
#                 # add_sol_pool(solution, pool_pre, pool_succ, pool_RI, pool_TT_walk, pool_n_route, pool_total_cost, pool_unserved, pool_penalty, pool_sol_chg_event)
                
#                 if sol.total_cost < sol_best.total_cost && (sol.n_route < K+1) 
#                     update_sol(sol_best, sol) 
#                     copy_route_info(best_e_r, best_cap_r, best_l_r, best_maxcap_r, best_dt_r, best_dt_r_back, e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back) 
#                     i_imp = 0
#                 end
#             else
#                 update_sol(solution, sol) 
#                 copy_route_info(e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, copy_e_r, copy_cap_r, copy_l_r, copy_maxcap_r, copy_dt_r, copy_dt_r_back)
#             end
#         end
#         if i_imp > 0            
#             T -= T_max / T_red
#             if T<0
#                 T = rand()*T_max
#                 if i_imp > n_imp * sol_best.n_route
#                     update_sol(sol, sol_best)
#                     update_sol(solution, sol_best)
#                     copy_route_info(e_r, cap_r, l_r, maxcap_r, dt_r, dt_r_back, best_e_r, best_cap_r, best_l_r, best_maxcap_r, best_dt_r, best_dt_r_back)
#                     i_imp = 0
#                 end
#             end            
#         end
#     end
#     return sol_best
# end


# function set_layered_graph(instance, global_parameter, time_table)
    
#     set_l_r = instance.set_l_r
#     bigM, detour_factor =global_parameter.bigM, global_parameter.detour_factor
#     duplcate_charger = global_parameter.duplcate_charger
#     n_c, n_bs, n_ts, n_charger, n_depot = instance.n_c, instance.n_bs, instance.n_ts, instance.n_charger, instance.n_depot;
#     parameter_energy, coord_cus, coord_chg, coord_bus, coord_transit, coord_depot = instance.parameter_energy, instance.coord_cus, instance.coord_chg, instance.coord_bus, instance.coord_transit, instance.coord_depot
#     T_max, t_s, speed_bus = instance.T_max, instance.t_s, instance.v_k
#     avg_wlk_speed  = global_parameter.avg_wlk_speed
#     vec_used_lyr=sort!(union(instance.set_l_r))  # used layers, e.g. 1,5,6
#     n_c_lyr=size(vec_used_lyr)[1] # total number of used layers
#     ts_time_table = time_table[:,1]  
#     ######################
#     # created layered graph
#     ######################
    
#     lyr_bus = zeros(Int32, n_bs * n_c_lyr) # layers of bus nodes  
#     G′=zeros(Int32, n_c_lyr, n_bs)  
#     D′=zeros(Int32, n_c_lyr, n_ts)  
#     G=zeros(Int32,n_bs)             
#     D=zeros(Int32,n_ts)            
#     S=zeros(Int32,n_charger)  
#     S′=zeros(Int32, n_charger * duplcate_charger)       
#     n_V = n_c_lyr*(n_bs+n_ts) + n_charger + n_depot + n_c + n_charger*duplcate_charger
#     lyrs_nodes = ones(Int32, n_V)
#     set_physical=zeros(Int32,n_V)
#     # generate all vertices and export to excel file for verification
#     V_DF = DataFrame(id=collect(Int32, 1:n_V), v_type=Vector{String}(undef, n_V), x=zeros(Float32,n_V), y=zeros(Float32,n_V),
#     lyr=zeros(Int32,n_V), ei = zeros(Float32, n_V),li=zeros(Float32,n_V),chg_power=zeros(Float32,n_V), set_physical = zeros(Float32,n_V), s_t = zeros(Float32,n_V), Li = zeros(Float32,n_V))
    
#     #each layer has n_bs + n_ts vertices (cs: customer, bs: bus, ts: transit stop)
#     _set_Gl=collect(1:n_bs)
#     _set_Dl=collect(n_bs+1:n_bs+n_ts)
#     nodes_bus = zeros(Int32, n_bs * n_c_lyr)
#     for _lyr in 1:n_c_lyr  
#         # dummy bus vertices    
#         temp = (_lyr-1)*(n_bs+n_ts)
#         temp_2 = collect((_lyr-1)*n_bs+1 : _lyr*n_bs)
#         _vec = collect(temp+1 : temp+ n_bs)
#         V_DF[_vec, 2]   = ["bus stop" for ss in 1:n_bs]
#         V_DF[_vec, 3:4] = coord_bus
#         V_DF[_vec, 5]   = [_lyr for ss in 1:n_bs]
#         G′[_lyr, :]  = _vec
#         nodes_bus[temp_2] = _vec
#         lyr_bus[temp_2] = [_lyr for i in 1:n_bs]
#         set_physical[_vec] = _set_Gl
#         # dummy transit vertices
#         _vec = collect( temp + n_bs+1 :  temp + n_bs + n_ts)
#         D′[_lyr,:]   = _vec
#         V_DF[_vec, 2]   = ["transit stop" for ss in 1:n_ts]
#         V_DF[_vec, 3:4] = coord_transit
#         V_DF[_vec, 5]   = [_lyr for ss in 1:n_ts]
#         _ei,_li     = time_table[vec_used_lyr[_lyr],2:3]
#         V_DF[_vec, 6]   = [ _ei for ss in 1:n_ts]
#         V_DF[_vec, 7]   = [ _li for ss in 1:n_ts]
#         set_physical[_vec] = _set_Dl
#         _vec = collect(temp+1 : temp+ n_bs+n_ts)
#         lyrs_nodes[_vec] = [_lyr for ss in 1:n_bs+n_ts]
#     end
    
#     G=G′[1,:] ; D=D′[1,:]
#     # @show(G′,D′)
#     # chargers
#     _count = n_c_lyr*(n_bs+n_ts); count_bus_ts = _count 
#     S = collect(_count+1: _count+n_charger)
#     set_physical[S] =collect(n_bs+n_ts+1:n_bs+n_ts+n_charger)
#     V_DF[S, 2] = ["charger" for ss in 1:n_charger]
#     V_DF[S, 3:4] = coord_chg ; V_DF[S, 8] = copy(parameter_energy.pw_charger)
#     #depot
#     v_depot = _count + n_charger+1 
#     set_physical[v_depot] = n_bs+n_ts+n_charger+1
#     V_DF[v_depot,2:7]=["depot",coord_depot[1],coord_depot[2],0, 0,T_max]
#     #customers
#     # reset the layer index for customers on used layers only
#     set_l_r_new = zeros(Int32, n_c)
#     for i in 1:n_c_lyr # n_c_lyr is the total number of used layers
#         temp = findall(x -> x==vec_used_lyr[i], set_l_r)
#         for j in eachindex(temp)
#             set_l_r_new[temp[j]] = i
#         end
#     end
#     # @show(set_l_r_new) # new layer idx of each customer, eg 1,2,3
#     set_d_r = [D′[set_l_r_new[i],ts_time_table[set_l_r[i]]] for i in 1:n_c] # set of d_r, d_r is the dummy transit node of cus r
#     set_d_r_physical = [n_bs + ts_time_table[set_l_r[i]] for i in 1:n_c ]
#     set_d_l =zeros(Int32,n_c_lyr)# node id for active transit node from lyr 1 to n_lyr 
#     set_s_l =zeros(Int32,n_c_lyr)# station id for active transit node from lyr 1 to n_lyr 
#     for i=1:n_c
#         _lr=set_l_r_new[i]
#         set_d_l[_lr]=set_d_r[i]
#         set_s_l[_lr]=ts_time_table[_lr] # transit station id
#     end
#     # @show( set_l_r_new, set_d_r, set_d_l,set_s_l)
#     R= collect(v_depot+1: v_depot+n_c)
#     V_DF[R[:,1],2]    = ["customer" for ss in 1:n_c]
#     V_DF[R[:,1], 3:4] = coord_cus; V_DF[R[:,1], 5] = set_l_r_new
    
#     # dummy charger nodes
#     tmp = v_depot+n_c
#     S′ = collect(tmp+1: tmp+duplcate_charger*n_charger)
#     for chg= 1:n_charger   
#         for j in 1:duplcate_charger
#             idx = tmp + (chg-1)*duplcate_charger+j
#             V_DF[idx, 2:4] = ["dummy charger", coord_chg[chg,1], coord_chg[chg,2]]
#             V_DF[idx, 8] = parameter_energy.pw_charger[chg]
#             set_physical[idx] = set_physical[S[chg]]
#         end
#     end
    
#     vec_st= zeros(Float32,size(V_DF)[1])
#     vec_st[nodes_bus] .= t_s
#     V_DF[:,9]  = set_physical; V_DF[:,10] = vec_st
    
#     # @show(set_physical)
#     # @show(G′);@show(D′);@show(G);@show(D);@show(S);@show(R);@show(set_l_r);@show(set_d_r);@show(set_s_l); @show(v_depot);
#     # @show(set_d_r_physical);show(set_physical);@show(lyrs_nodes)
    
#     ###############################################################################
#     # identify compatible layers, update on 10-18 
#     #############################################################################"
#     lyrs_compatible= falses(n_c_lyr,n_c_lyr) 
#     for i in 1:n_c_lyr
#         lyrs_compatible[i,i]= true
#     end
#     for i in 1:n_c_lyr-1 # lower lyer
#         ts_i= time_table[vec_used_lyr[i],1]
#         for j in i+1:n_c_lyr # upper lyer
#             ts_j=time_table[vec_used_lyr[j],1]
#             _temp = coord_transit[ts_i,:]-coord_transit[ts_j,:]
#             _dist_tt= sqrt(sum(_temp.^2)) / speed_bus
#             # @show(vec_used_lyr[i],vec_used_lyr[j], time_table[vec_used_lyr[i],3],time_table[vec_used_lyr[j],2],_dist_tt)
#             # ei+tij =<lj && li+tij >= ej (the criteria if two layers are compatible)
#             if  time_table[vec_used_lyr[i],2] + _dist_tt <= time_table[vec_used_lyr[j],3] && (time_table[vec_used_lyr[i],3] + _dist_tt >= time_table[vec_used_lyr[j],2])
#                 lyrs_compatible[i,j]=true; lyrs_compatible[j,i]=true
#             end
#         end
#     end
#     lyrs_compatible_extend = falses(n_c_lyr,n_c_lyr)
#     lyrs_compatible_extend[1,:]= lyrs_compatible[1,:]
#     tmp1= findall(x->x==true, lyrs_compatible[1,:])
#     tmp2= findall(x->x==true, lyrs_compatible[2,:])
#     tmp= union(tmp1, tmp2)
#     lyrs_compatible_extend[2, tmp] .= true
#     if n_c_lyr >2
#         for i in 1:n_c_lyr-2
#             tmp3= findall(x->x==true, lyrs_compatible[i+2,:])
#             tmp= union(tmp1, tmp2, tmp3)
#             lyrs_compatible_extend[i+2, tmp] .= true
#             tmp1= tmp2
#             tmp2= tmp3
#         end
#     end
#     coord_V = hcat(V_DF.x,V_DF.y); ei_li_V = hcat(V_DF.ei,V_DF.li)
#     G_D=vcat(G,D) ; V_l=vcat(G,D,S,v_depot); n_node_G_D = length(G_D)
#     # please be aware that dist_GD, dist, dist_all are measured in minutes, not in kilometer 
    
    
#     dist_GD, dist_all= comp_dist(coord_V, G_D, V_l, speed_bus)  
    
    
#     # walking time from customer's origin to bus stop
#     dist_r_bus = comp_dist_cus_bus(coord_V, nodes_bus, R, set_l_r_new , lyr_bus, avg_wlk_speed, bigM) # in minutes
#     Li_tmp = zeros(Float32,size(V_DF)[1])
#     setup_ei_li_V_DF(V_DF, n_c_lyr, set_d_l, set_physical, detour_factor, dist_GD, t_s, T_max, n_bs, n_ts)
#     set_Li_V(Li_tmp, n_bs, n_c_lyr, set_d_l, set_physical, detour_factor, dist_GD, n_ts)
#     V_DF[:,11] = Li_tmp
    
#     lgraph=LGraph(n_c_lyr, vec_used_lyr, nodes_bus, G′, D′, G, D, S, v_depot, R, S′, set_d_l, set_s_l, set_l_r, set_l_r_new, set_d_r, set_d_r_physical,
#     set_physical, lyr_bus, lyrs_nodes, lyrs_compatible, lyrs_compatible_extend, dist_GD, dist_all, dist_r_bus, coord_V)
    
#     return lgraph, V_DF, dist_r_bus
    
# end


# function setup_lyr_id(set_l_r, is_last_mile)
#     set_l_r_new      = zeros(Int32, 1000)
#     is_last_mile_lyr = falses(1000)
#     count = 0; dict_lyr=Dict();dict_tmp=Dict()  
#     for i in 1:n_c 
#         key=(set_l_r[i],is_last_mile[i])
#         if !haskey(dict_tmp, key)
#             count += 1
#             dict_lyr[count]=key
#             dict_tmp[key]=count
#             set_l_r_new[count] = count
#             is_last_mile_lyr[count] = is_last_mile[i]
#         end
#     end
#     return set_l_r_new[1:count], is_last_mile_lyr[1:count], dict_lyr
# end

# # insert charging events on a route. First, identify possible charging positions of the route, then insert randomly the charging event 
# # on a random position and a random (random policy)/best(greedy policy) charger 
# function insert_charging_route_test(solution::Solution, r, route::Vector{Int32}, dist_all, ec_route, fast_chg::Fast_chargers, parameter_energy, set_physical, dt_r, dt_r_back, flag_init_sol)
    
#     end_depot = 2*n_cus+2; v=Int32[]
#     max_num_chg = 4 # max num of charging events limits for a route 
#     success_chg= false; remain_e_chg =0; idx_loc =0
#     n_fast_chg = fast_chg.n_fast_chg
#     pos_insert_chg = zeros(Int32, 30, 3) # store the pair nodes (precedent and successive nodes) of the location of the charger
#     n_pos_insert_chg = 1; t_v_v1 = 0; t_arr_node_after_chg = 0; idx_pre_chg_pos =0
#     info_chg = zeros(Float32, 2+4*4) # store at most 4 charging operations  
#     # find the positions (nodes) where the SOC (state of charge) of the vehicle is insufficient  
#     pos_insert_chg[n_pos_insert_chg,:] = [route[1], route[2], 1]
#     veh_id = solution.RI[r,5]
#     max_ec = parameter_energy.max_ec[veh_id]
#     i_ec_violate= 0 # will alwayse find a positive position 
#     dt_r_eng = dt_r[route] .* parameter_energy.beta[veh_id]
#     i_ec_violate = findfirst(x->x-max_ec>0, dt_r_eng)
#     # @show(r, route, dt_r[route], dt_r_eng, max_ec,i_ec_violate )
#     # need to check whether all customers on board have beed drop-offed to be a legitimate candidate charging position
#     for (idx, v) in enumerate(route[1:i_ec_violate-1])
#         if  v > n_cus+1 && (route[idx+1] < n_cus+2 || route[idx+1] == end_depot)
#             if (idx-1)%2 ==0 # check if all customers on board have beed drop-offed
#                 n_pos_insert_chg += 1
#                 pos_insert_chg[n_pos_insert_chg, :] = [v, route[idx+1], idx]
#             end
#         end
#     end 
  
#     # @show(route, pos_insert_chg)
#     #randomly select a fast charger
#     idx_first_loc =  1 
#     # chg is the fast charger index in the list of all fast chargers (vec_idx_chg)
#     v= pos_insert_chg[idx_first_loc,1:2]
#     rand_policy = false
    
#     if rand_policy
#         idx_chg = rand(collect(1:n_fast_chg))
#     else
#         idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg)
#     end
#     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)

#     if (dt_r[v[1]] + t_access_1) * parameter_energy.beta[veh_id] > parameter_energy.max_ec[veh_id] && idx_first_loc > 1
#         idx_first_loc -= 1
#         v= pos_insert_chg[idx_first_loc,1:2]
#         t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
#     end

#     vec_delta_e  =  parameter_energy.E_min[veh_id] - (parameter_energy.E_init[veh_id] - (ec_route + parameter_energy.beta[veh_id] * (t_access_1+ t_access_2 - t_v_v1)))
#     # vec_delta_e <= 0 && error("vec_delta_e <=0 error!!")
#     t_chg = vec_delta_e / fast_chg.pw[idx_chg]
#     idxs_chg_loc = pos_insert_chg[1:n_pos_insert_chg, 3]
#         # @show(idx_loc, v, t_chg, idxs_chg_loc)
#     forward_time_slack, vec_t_start = get_forward_time_slack(route, ei, li, s_t, dist, n_cus, qi, idxs_chg_loc,dist_orig, darp)
#     # @show(forward_time_slack)
#     count_chg = 1; next_recharge = true; 
#     n_candicate_loc =  n_pos_insert_chg - idx_first_loc +1
#     # @show(pos_insert_chg, idxs_chg_loc)
#     remain_e_chg = vec_delta_e
#     if t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[idx_first_loc]
#         # @show("t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[1] !")
#         additional_time = t_access_1 + t_access_2 - t_v_v1 +  t_chg
#         t_start_chg = vec_t_start[1] + t_access_1 + rand(1)[1]* (forward_time_slack[idx_first_loc]-t_chg) # random start time
#         info_chg[1] += additional_time;  info_chg[2] +=1
#         info_chg[3:6] = [v[1], idx_chg, t_start_chg, t_start_chg+t_chg]
#         success_chg = true
#         # @show( t_chg + t_access_1+ t_access_2 - t_v_v1 , forward_time_slack[1], info_chg[1:6])
#     else        
#         while next_recharge  
#             if count_chg >= n_candicate_loc  || count_chg >= max_num_chg
#                 success_chg = false
#                 break
#             else
#                 # @show(count_chg,v[1])
#                 next_recharge, remain_e_chg, t_arr_node_after_chg =  schedule_chg(info_chg, darp, remain_e_chg, r, route, veh_id, forward_time_slack, vec_t_start, idx_pre_chg_pos, t_arr_node_after_chg, fast_chg, pos_insert_chg, idx_first_loc, count_chg, v[1], idx_chg, t_access_1, t_access_2, t_v_v1)
#                 idx_loc_pre = idx_first_loc - 1 + count_chg
#                 idx_pre_chg_pos = pos_insert_chg[idx_loc_pre, 3]
#                 # next location
#                 count_chg += 1
#                 if count_chg < n_candicate_loc +1
#                     v = pos_insert_chg[idx_first_loc-1+count_chg,1:2]
#                     if rand_policy    
#                         idx_chg = rand(collect(1:n_fast_chg))
#                     else
#                         idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg)
#                     end
#                     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
#                 end
#             end
#         end
#         # @show( "! next_recharge && success_chg = true : ")
#         ! next_recharge && (success_chg = true)
#     end
#     if success_chg              
#         solution.vec_chg_event[r,:] = info_chg
#         # @show( success_chg, r, route, info_chg ) 
#         # check_energy(solution, route, r, parameter_energy, fast_chg)
#         return true,  pos_insert_chg, forward_time_slack, vec_t_start
#     else     
#         return false, pos_insert_chg, forward_time_slack, vec_t_start
#     end
# end

# # insert charging events on a route. First, identify possible charging positions of the route, then insert randomly the charging event 
# # on a random position and a random (random policy)/best(greedy policy) charger 
# function insert_charging_route_exchange(solution::Solution, darp, r, route::Vector{Int32}, dist_all, ec_route, fast_chg::Fast_chargers, parameter_energy, set_physical, dt_r, dt_r_back, flag_init_sol)
    
        
#     dist_orig, n_cus, end_depot = darp.dist_orig, darp.n_cus, darp.end_depot
#     ei, li, s_t, dist, qi = darp.ei, darp.li, darp.s_t, darp.dist, darp.qi

#     v=Int32[]
#     max_num_chg = 4 # max num of charging events limits for a route 
#     success_chg= false; remain_e_chg =0; idx_loc =0
#     n_fast_chg = fast_chg.n_fast_chg
#     pos_insert_chg = zeros(Int32, 30, 3) # store the pair nodes (precedent and successive nodes) of the location of the charger
#     n_pos_insert_chg = 1; t_v_v1 = 0; t_arr_node_after_chg = 0; idx_pre_chg_pos =0
#     info_chg = zeros(Float32, 2+4*4) # store at most 4 charging operations  
#     # find the positions (nodes) where the SOC (state of charge) of the vehicle is insufficient  
#     pos_insert_chg[n_pos_insert_chg,:] = [route[1], route[2], 1]
#     veh_id = solution.RI[r,5]
#     max_ec = parameter_energy.max_ec[veh_id]
#     i_ec_violate= 0 # will alwayse find a positive position 
#     dt_r_eng = dt_r[route] .* parameter_energy.beta[veh_id]
#     i_ec_violate = findfirst(x->x-max_ec>0, dt_r_eng)
#     # @show(r, route, dt_r[route], dt_r_eng, max_ec,i_ec_violate )
#     # need to check whether all customers on board have beed drop-offed to be a legitimate candidate charging position
#     for (idx, v) in enumerate(route[1:i_ec_violate-1])
#         if  v > n_cus+1 && (route[idx+1] < n_cus+2 || route[idx+1] == end_depot)
#             if (idx-1)%2 ==0 # check if all customers on board have beed drop-offed
#                 n_pos_insert_chg += 1
#                 pos_insert_chg[n_pos_insert_chg, :] = [v, route[idx+1], idx]
#             end
#         end
#     end 
  
#     # @show(route, pos_insert_chg)
#     #randomly select a fast charger
#     idx_first_loc =  1 # from the first location
#     # chg is the fast charger index in the list of all fast chargers (vec_idx_chg)
#     v= pos_insert_chg[idx_first_loc,1:2]
#     rand() > 0.5 ? rand_policy =true : rand_policy = false
#     flag_init_sol && (rand_policy = true)
#     if rand_policy
#         idx_chg = rand(collect(1:n_fast_chg))
#     else
#         idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg)
#     end
#     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)

#     if (dt_r[v[1]] + t_access_1) * parameter_energy.beta[veh_id] > parameter_energy.max_ec[veh_id] && idx_first_loc > 1
#         idx_first_loc -= 1
#         v= pos_insert_chg[idx_first_loc,1:2]
#         t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
#     end

#     vec_delta_e  =  parameter_energy.E_min[veh_id] - (parameter_energy.E_init[veh_id] - (ec_route + parameter_energy.beta[veh_id] * (t_access_1+ t_access_2 - t_v_v1)))
#     # vec_delta_e <= 0 && error("vec_delta_e <=0 error!!")
#     t_chg = vec_delta_e / fast_chg.pw[idx_chg]
#     idxs_chg_loc = pos_insert_chg[1:n_pos_insert_chg, 3]
#         # @show(idx_loc, v, t_chg, idxs_chg_loc)
#     forward_time_slack, vec_t_start = get_forward_time_slack(route, ei, li, s_t, dist, n_cus, qi, idxs_chg_loc,dist_orig, darp)
#     # @show(forward_time_slack)
#     count_chg = 1; next_recharge = true; 
#     n_candicate_loc =  n_pos_insert_chg - idx_first_loc +1
#     # @show(pos_insert_chg, idxs_chg_loc)
#     remain_e_chg = vec_delta_e
#     if t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[idx_first_loc]
#         # @show("t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[1] !")
#         additional_time = t_access_1 + t_access_2 - t_v_v1 +  t_chg
#         t_start_chg = vec_t_start[1] + t_access_1 + rand(1)[1]* (forward_time_slack[idx_first_loc]-t_chg) # random start time
#         info_chg[1] += additional_time;  info_chg[2] +=1
#         info_chg[3:6] = [v[1], idx_chg, t_start_chg, t_start_chg+t_chg]
#         success_chg = true
#         # @show( t_chg + t_access_1+ t_access_2 - t_v_v1 , forward_time_slack[1], info_chg[1:6])
#     else        
#         while next_recharge  
#             if count_chg >= n_candicate_loc  || count_chg >= max_num_chg
#                 success_chg = false
#                 break
#             else
#                 # @show(count_chg,v[1])
#                 next_recharge, remain_e_chg, t_arr_node_after_chg =  schedule_chg(info_chg, darp, remain_e_chg, r, route, veh_id, forward_time_slack, vec_t_start, idx_pre_chg_pos, t_arr_node_after_chg, fast_chg, pos_insert_chg, idx_first_loc, count_chg, v[1], idx_chg, t_access_1, t_access_2, t_v_v1)
#                 idx_loc_pre = idx_first_loc - 1 + count_chg
#                 idx_pre_chg_pos = pos_insert_chg[idx_loc_pre, 3]
#                 # next location
#                 count_chg += 1
#                 if count_chg < n_candicate_loc +1
#                     v = pos_insert_chg[idx_first_loc-1+count_chg,1:2]
#                     if rand_policy    
#                         idx_chg = rand(collect(1:n_fast_chg))
#                     else
#                         idx_chg = get_chg_greedy(v, dist_all, set_physical, fast_chg)
#                     end
#                     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
#                 end
#             end
#         end
#         # @show( "! next_recharge && success_chg = true : ")
#         ! next_recharge && (success_chg = true)
#     end
#     if success_chg              
#         solution.vec_chg_event[r,:] = info_chg
#         # @show(" success_chg, r, route, info_chg ) 
#         # check_energy(solution, route, r, parameter_energy, fast_chg)
#         return true,  pos_insert_chg, forward_time_slack, vec_t_start
#     else     
#         return false, pos_insert_chg, forward_time_slack, vec_t_start
#     end
# end


# function verify_solution_detail(solution::Solution, Q, instance, darp, lgraph, fast_chg)

#     start_depot = 1
#     end_depot = 2*n_cus+2
#     # disp arr_t and dep_t of of routes
#     verify_solution(solution, Q)
#     # for r in 1:solution.n_route
#     #     route = get_route(solution, r, start_depot, end_depot)
#     #     disp_route(route, r, ei, li, s_t, dist, Q, Li, n_cus, TH)
#     # end

# end


# insert charging operations on a route given the assigned charger 
# function insert_charging_assigned(solution::Solution, r, route::Vector{Int32}, idx_chg, dist_all, ec_route, fast_chg::Fast_chargers, parameter_energy, set_physical, nodes, dt_r, dt_r_back)
    
#     end_depot = 2*n_cus+2; v=Int32[]
#     max_num_chg = 4 # max num of charging events limits for a route 
#     success_chg= false; remain_e_chg =0
#     pos_insert_chg = zeros(Int32, 30, 3) # store the pair nodes (precedent and successive nodes) of the location of the charger
#     n_pos_insert_chg = 1; t_v_v1 = 0; t_arr_node_after_chg = 0; idx_pre_chg_pos =0
#     info_chg = zeros(Float32, 2+4*4) # store at most 4 charging operations  
#     # find the positions (nodes) where the SOC (state of charge) of the vehicle is insufficient  
#     pos_insert_chg[n_pos_insert_chg,:] = [route[1], route[2], 1]
#     veh_id = solution.RI[r,5]
#     max_ec = parameter_energy.max_ec[veh_id]
#     dt_r_eng = dt_r[route] .* parameter_energy.beta[veh_id]
#     i_ec_violate = findfirst(x->x-max_ec>0, dt_r_eng)
#     for (idx, v) in enumerate(route[1:i_ec_violate-1])
#         if  v > n_cus+1 && (route[idx+1] < n_cus+2 || route[idx+1] == end_depot)
#             if (idx-1)%2 ==0 # check if all customers on board have beed drop-offed
#                 n_pos_insert_chg += 1
#                 pos_insert_chg[n_pos_insert_chg, :] = [v, route[idx+1], idx]
#             end
#         end
#     end 
#     idx_first_loc =  rand(collect(1:n_pos_insert_chg))
#     v= pos_insert_chg[idx_first_loc,1:2]

#     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)

#     if (dt_r[v[1]] + t_access_1) * parameter_energy.beta[veh_id] > parameter_energy.max_ec[veh_id] && idx_first_loc > 1
#         idx_first_loc -= 1
#         v= pos_insert_chg[idx_first_loc,1:2]
#         t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
#     end

#     vec_delta_e  =  parameter_energy.E_min[veh_id] - (parameter_energy.E_init[veh_id] - (ec_route + parameter_energy.beta[veh_id] * (t_access_1+ t_access_2 - t_v_v1)))
#     t_chg = vec_delta_e / fast_chg.pw[idx_chg]
#     idxs_chg_loc = pos_insert_chg[1:n_pos_insert_chg, 3]
#     forward_time_slack, vec_t_start = get_forward_time_slack(route, ei, li, s_t, dist, n_cus, qi, idxs_chg_loc, dist_orig, darp)
#     count_chg = 1; next_recharge = true; 
#     n_candicate_loc =  n_pos_insert_chg - idx_first_loc +1 
#     remain_e_chg = vec_delta_e
#     if t_chg + t_access_1+ t_access_2 - t_v_v1 < forward_time_slack[idx_first_loc] 
#         additional_time = t_access_1 + t_access_2 - t_v_v1 +  t_chg
#         t_start_chg = vec_t_start[1] + t_access_1 + rand(1)[1]* (forward_time_slack[idx_first_loc]-t_chg) # random start time
#         info_chg[1] += additional_time;  info_chg[2] +=1
#         info_chg[3:6] = [v[1], idx_chg, t_start_chg, t_start_chg+t_chg]
#         success_chg = true 
#     else        
#         while next_recharge  
#             if count_chg > n_candicate_loc  || count_chg > max_num_chg
#                 success_chg = false
#                 break
#             else
#                 next_recharge, remain_e_chg, t_arr_node_after_chg =  schedule_chg(info_chg, darp, remain_e_chg, r, route, veh_id, forward_time_slack, vec_t_start, idx_pre_chg_pos, t_arr_node_after_chg, fast_chg, pos_insert_chg, idx_first_loc, count_chg, v[1], idx_chg, t_access_1, t_access_2, t_v_v1)
#                 idx_loc_pre = idx_first_loc - 1 + count_chg
#                 idx_pre_chg_pos = pos_insert_chg[idx_loc_pre, 3]
#                 # next location
#                 count_chg += 1
#                 if count_chg < n_candicate_loc +1
#                     v = pos_insert_chg[idx_first_loc-1+count_chg,1:2]
#                     t_access_1, t_access_2, t_v_v1 =  compute_t_access(v, idx_chg, dist_all, set_physical, fast_chg)
#                 end
#             end
#         end 
#         ! next_recharge && (success_chg = true)
#     end
#     if success_chg              
#         solution.vec_chg_event[r,:] = info_chg 
#         return true
#     else     
#         return false
#     end
# end
