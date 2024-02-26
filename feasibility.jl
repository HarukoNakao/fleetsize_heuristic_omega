function calculate_ads(route::Vector{Int32}, t_start, ei, li, s_t, dist, n_cus, qi)
    # 1 and 2n+2 are the depots, 2...n+1 (pickup), n+2...2n+1 dropoff
    N = length(route)
    A = zeros(Float32, N)  #Arrival time at vertex i 
    B = zeros(Float32, N)  #Beginning of service at vertex i 
    W = zeros(Float32, N)  #Wait time before starting service at vertex i 
    D = zeros(Float32, N)  #Departure time at vertex i 
    Load = zeros(Int32, N) #cummulative laod after leaving vertex i 
    RT = zeros(Float32, N) #ride time, only specified for pickup vertices else 0 
    A[1] = t_start
    B[1] = max(t_start, ei[route[1]])
    D[1] = B[1] + s_t[route[1]]
    W[1] = B[1] - A[1]
    for idx in 2:N
        A[idx] = D[idx-1] + dist[route[idx-1], route[idx]]
        B[idx] = max(A[idx], ei[route[idx]])
        W[idx] = B[idx] - A[idx]
        D[idx] = B[idx] + s_t[route[idx]]
        Load[idx] = Load[idx-1] + qi[route[idx]]
    end
    for idx in 2:N-1  # pickup or dropoff
        if route[idx] > n_cus + 1
            idx2 = findfirst(x -> x == route[idx] - n_cus, route)
            # RT[idx] = A[idx] - D[idx2]
            RT[idx] = B[idx] - D[idx2]
        end
    end

    return A, B, W, D, Load, RT
end
function calculate_t_wait_no_chg(route::Vector{Int32}, t_start, ei, li, s_t, dist_orig, n_cus)
    
    t_wait = 0 
    N = length(route)
    A = zeros(Float32, N)  #Arrival time at vertex i 
    B = zeros(Float32, N)  #Beginning of service at vertex i 
    W = zeros(Float32, N)  #Wait time before starting service at vertex i 
    D = zeros(Float32, N)  #Departure time at vertex i  
    A[1] = t_start
    B[1] = max(t_start, ei[route[1]])
    D[1] = B[1] + s_t[route[1]]
    W[1] = B[1] - A[1]
    for idx in 2:N
        A[idx] = D[idx-1] + dist_orig[route[idx-1], route[idx]]
        B[idx] = max(A[idx], ei[route[idx]])
        W[idx] = B[idx] - A[idx]
        D[idx] = B[idx] + s_t[route[idx]] 
        if W[idx] > 0 && route[idx] > n_cus+1
            idx2 = findfirst(x -> x == route[idx] - n_cus, route)
            if B[idx2]+W[idx] > li[route[idx2]]
                t_wait += (W[idx] - (li[route[idx2]]-B[idx2])) # delay the beginning of service time at li 
            end
        end 
    end
    return t_wait
end

function calculate_t_wait_chg(route::Vector{Int32}, t_start, ei_route, li_route, s_t_route, dist_orig, n_cus, fast_chg)
  
    t_wait = 0
    n_node_darp = 2 * n_cus + 2
    N = length(route)
    A = zeros(Float32, N)  #Arrival time at vertex i 
    B = zeros(Float32, N)  #Beginning of service at vertex i 
    W = zeros(Float32, N)  #Wait time before starting service at vertex i 
    D = zeros(Float32, N)  #Departure time at vertex i   
    A[1] = t_start
    B[1] = max(t_start, ei_route[1])
    D[1] = B[1] + s_t[route[1]]
    W[1] = B[1] - A[1]
    for idx in 2:N
        if route[idx-1] > n_node_darp
            v1, v2 = fast_chg.v_physical_chg[route[idx-1]-10000], set_physical[route[idx]]
            A[idx] = D[idx-1] + dist_all[v1, v2] 
            B[idx] = max(A[idx], ei_route[idx])
            W[idx] = B[idx] - A[idx]
            D[idx] = B[idx] + s_t_route[idx] 
        elseif route[idx] > n_node_darp
            v1, v2 = set_physical[route[idx-1]], fast_chg.v_physical_chg[route[idx]-10000]
            A[idx] = D[idx-1] + dist_all[v1, v2] 
            B[idx] = max(A[idx], ei_route[idx])
            W[idx] = B[idx] - A[idx]
            D[idx] = B[idx] + s_t_route[idx] 
        else
            A[idx] = D[idx-1] + dist_orig[route[idx-1], route[idx]]
            B[idx] = max(A[idx], ei_route[idx])
            W[idx] = B[idx] - A[idx]
            D[idx] = B[idx] + s_t_route[idx] 
            if W[idx] > 0 && route[idx] > n_cus+1
                idx2 = findfirst(x -> x == route[idx] - n_cus, route)
                if B[idx2]+W[idx] > li_route[idx2]
                    t_wait += (W[idx] - (li_route[idx2]-B[idx2])) # delay the beginning of service time at li 
                end
            end
        end
    end 
    return t_wait
end

# for tw and ride time check when involving charging operations on a route
function calculate_ads_t_wait_chg(route::Vector{Int32}, t_start, ei_route, s_t_route, dist_orig, dist_all, n_cus, fast_chg)
    # 1 and 2n+2 are the depots, 2...n+1 (pickup), n+2...2n+1 dropoff
    n_node_darp = 2 * n_cus + 2
    N = length(route)
    A = zeros(Float32, N)  #Arrival time at vertex i 
    B = zeros(Float32, N)  #Beginning of service at vertex i 
    W = zeros(Float32, N)  #Wait time before starting service at vertex i 
    D = zeros(Float32, N)  #Departure time at vertex i  
    RT = zeros(Float32, N) #ride time, only specified for pickup vertices else 0 
    A[1] = t_start
    B[1] = max(t_start, ei_route[1])
    D[1] = B[1] + s_t[route[1]]
    W[1] = B[1] - A[1]
    for idx in 2:N
        if route[idx-1] > n_node_darp
            v1, v2 =  fast_chg.v_physical_chg[route[idx-1]-10000], set_physical[route[idx]]
            A[idx] = D[idx-1] + dist_all[v1, v2] 
        elseif route[idx] > n_node_darp
            v1, v2 = set_physical[route[idx-1]],  fast_chg.v_physical_chg[route[idx]-10000]
            A[idx] = D[idx-1] + dist_all[v1, v2] 
        else
            A[idx] = D[idx-1] + dist_orig[route[idx-1], route[idx]]
        end
        B[idx] = max(A[idx], ei_route[idx])
        W[idx] = B[idx] - A[idx]
        D[idx] = B[idx] + s_t_route[idx] 
    end
    for idx in 2:N-1  # pickup or dropoff
        if route[idx] > n_cus + 1 && (route[idx] < 2 * n_cus + 2)
            idx2 = findfirst(x -> x == route[idx] - n_cus, route)
            # RT[idx] = A[idx] - D[idx2]
            RT[idx] = B[idx] - D[idx2]
        end
    end

    return A, B, W, D, RT
end
                             
# for tw and ride time check when involving charging operations on a route
function calculate_ads_chg(fast_chg, darp, route::Vector{Int32}, t_start, ei_route, li_route, s_t_route, dist_orig, n_cus, qi_route, set_physical, dist_all)
 
    s_t = darp.s_t
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
        elseif route[idx] > n_node_darp
            v1, v2 = set_physical[route[idx-1]], fast_chg.v_physical_chg[route[idx]-10000]
            A[idx] = D[idx-1] + dist_all[v1, v2] 
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
            RT[idx] = B[idx] - D[idx2]
            # RT[idx] = A[idx] - D[idx2]
        end
    end

    return A, B, W, D, Load, RT
end

function check_tw_cap_chg(route::Vector{Int32}, B, li_route, Load, ri, Q)

    for i in 2:length(route)
        if B[i] > li_route[i] || Load[i] > Q[ri]
            return false
        end
    end
    return true
end

#@@@@@@@@@@@@@@@@@@ PARTIAL CHECKS #@@@@@@@@@@@@@@@@@@
function check_tw_cap(route::Vector{Int32}, B, li, Load, ri, Q)

    for i in 2:length(route)
        if round(B[i] - li[route[i]], digits=DIGITS) > 0 || Load[i] > Q[ri]
            return false
        end
    end
    return true
end

function check_ridetime_chg(route::Vector{Int32}, RT, Li_route)

    for i in 1:length(route)
        if round(RT[i] - Li_route[i], digits=DIGITS) > 0
            return false
        end
    end
    return true
end


function check_ridetime(route::Vector{Int32}, RT, Li)

    for i in 1:length(route)
        if round(RT[i] - Li[route[i]], digits=DIGITS) > 0
            return false
        end
    end
    return true
end

function check_routeduration(route::Vector{Int32}, B, TH)

    round(B[end] - B[1], digits=DIGITS) > TH && return false
    return true
end

function calculate_ads_chg_schedule(route::Vector{Int32}, t_start, darp)

    ei, s_t, dist_orig = darp.ei, darp.s_t, darp.dist_orig
    # 1 and 2n+2 are the depots, 2...n+1 (pickup), n+2...2n+1 dropoff
    N = length(route)
    A = zeros(Float32, N)  #Arrival time at vertex i 
    B = zeros(Float32, N)  #Beginning of service at vertex i 
    W = zeros(Float32, N)  #Wait time before starting service at vertex i 
    D = zeros(Float32, N)  #Departure time at vertex i 
    A[1] = t_start
    B[1] = max(t_start, ei[route[1]])
    D[1] = B[1] + s_t[route[1]]
    W[1] = B[1] - A[1]
    for idx in 2:N
        A[idx] = D[idx-1] + dist_orig[route[idx-1], route[idx]] # here is the difference (using non penalized dist matrix) with calculate_ads
        B[idx] = max(A[idx], ei[route[idx]])
        W[idx] = B[idx] - A[idx]
        D[idx] = B[idx] + s_t[route[idx]]
    end
    return W, B
end

# get the slack time as charging time constraint which satisfies the TW and ride time constraints of the (remaining) route
# t_start is the time arriving at the first node of route
function get_forward_time_slack_revise(route::Vector{Int32}, idx_node, W, B, li)

    N = length(route)  
    Fi = calc_Fi_chg_scheduling(idx_node, route, W, B, li)
    Wp = sum(W[idx_node+1:N-1])
    forward_time_slack = min(Fi, Wp)
    return forward_time_slack, B[idx_node]
end

function get_forward_time_slack(route::Vector{Int32}, ei, li, s_t, dist, n_cus, qi, idxs_chg_loc, dist_orig, darp)
                                              
    N = length(route)
    n_loc = length(idxs_chg_loc)
    forward_time_slack = zeros(Float32, n_loc)
    # t_start = ei[route[1]]
    # A, B, W, D, Load, RT = calculate_ads(route, t_start, ei, li, s_t, dist, n_cus, qi)
    W, B = calculate_ads_chg_schedule(route, ei[route[1]], darp)
    for idx in 1:n_loc
        i = idxs_chg_loc[idx]
        Fi = calc_Fi_chg_scheduling(i, route, W, B, li)
        Wp = sum(W[i+1:N-1])
        forward_time_slack[idx] = min(Fi, Wp)
        # @show(i,Fi,Wp,forward_time_slack[idx] )
    end
    # @show(forward_time_slack, B[idxs_chg_loc])
    return forward_time_slack, B[idxs_chg_loc]
end

 
# calculate Fi for charging scheduling
# see the paper of Cordeau and Laporte, 2003, Transp. Res. Part B
function calc_Fi_chg_scheduling(i_start, route::Vector{Int32}, W::Vector{Float32}, B::Vector{Float32}, li)
    #Li needs to be changed to a vector of Li for which each user's max ridetime is different
    N = length(route)
    _bigM = 100000
    slacks = ones(Float32, N - i_start + 1) * _bigM
    idx = 1 
    for j in i_start:N
        term1 = sum(W[i_start+1:j]) 
        term2 = li[route[j]] - B[j]
        slack = term1 + term2
        slacks[idx] = slack
        idx += 1
    end
    return minimum(slacks)
end

# for tw and ride time check when there are charging operations on a route
function calc_Fi_chg(i_start, route::Vector{Int32}, W::Vector{Float32}, B::Vector{Float32}, RT::Vector{Float32}, Li_route, li_route)
    N = length(route)
    _bigM = 100000
    slacks = ones(Float32, N - i_start + 1) * _bigM
    idx = 1
    for j in i_start:N 
        term1 = sum(W[i_start+1:j]) # modified on 10-3-2022 
        term2 = max(0, min(li_route[j] - B[j], Li_route[j] - RT[j])) #or RT[j] equal to zero if j-n not visited before i 
        slack = term1 + term2
        slacks[idx] = slack
        idx += 1
    end
    return minimum(slacks)
end


#@@@@@@@@@@@@@@@@@@ FORWARD TIME SLACK #@@@@@@@@@@@@@@@@@@
function calc_Fi(i_start, route::Vector{Int32}, W::Vector{Float32}, B::Vector{Float32}, RT::Vector{Float32}, Li::Vector{Float32}, li::Vector{Float32})
    #Li needs to be changed to a vector of Li for which each user's max ridetime is different
    N = length(route)
    _bigM = 100000
    slacks = ones(Float32, N - i_start + 1) * _bigM
    idx = 1
    for j in i_start:N 
        term1 = sum(W[i_start+1:j]) # modify on 10-3-2022 
        term2 = max(0, min(li[route[j]] - B[j], Li[route[j]] - RT[j])) #or RT[j] equal to zero if j-n not visited before i 
        slack = term1 + term2
        slacks[idx] = slack
        idx += 1
    end
    return minimum(slacks)

end


function eight_step(route::Vector{Int32}, ri, ei, li, s_t, dist, Q, Li, n_cus, TH, qi) 

    N_NODE = 2 * n_cus + 2
    N = length(route)
    t_start = ei[route[1]]
    A, B, W, D, Load, RT = calculate_ads(route, t_start, ei, li, s_t, dist, n_cus, qi)
    if !check_tw_cap(route, B, li, Load, ri, Q)
        return false
    end

    F0 = calc_Fi(1, route, W, B, RT, Li, li)
    Wp = sum(W[2:N-1])
    t_start_new = ei[route[1]] + min(F0, Wp)
    A, B, W, D, Load, RT = calculate_ads(route, t_start_new, ei, li, s_t, dist, n_cus, qi)

    if check_ridetime(route, RT, Li) && check_routeduration(route, B, TH)
        return true
    end

    for i in 2:N-1
        if 1 < route[i] <= n_cus + 1
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
                if n_cus + 1 < route[j] < N_NODE
                    idx2 = findfirst(x -> x == route[j] - n_cus, route)
                    # RT[j] = A[j] - D[idx2]
                    RT[j] = B[j] - D[idx2]
                end
            end
        end
    end
    #check whether the drive time is respected for all dropoff vertices after i, if so the route is feasible 
    if check_routeduration(route, B, TH) && check_ridetime(route, RT, Li)
        return true
    else
        return false
    end

end

# tailed version to check feasibility of a portion of a route
# t_start is the beginning time of service at the first node of the route
function eight_step_light(route::Vector{Int32}, t_start, ri, ei, li, s_t, dist, Q, Li, n_cus, TH)

    N_NODE = 2 * n_cus + 2
    N = length(route)
    A, B, W, D, Load, RT = calculate_ads(route, t_start, ei, li, s_t, dist, n_cus, qi)
    if !check_tw_cap(route, B, li, Load, ri, Q)
        return false
    end

    F0 = calc_Fi(1, route, W, B, RT, Li, li)
    Wp = sum(W[2:N-1])
    t_start_new = ei[route[1]] + min(F0, Wp)
    A, B, W, D, Load, RT = calculate_ads(route, t_start_new, ei, li, s_t, dist, n_cus, qi)

    if check_ridetime(route, RT, Li)
        return true
    end

    for i in 2:N-1
        if 1 < route[i] <= n_cus + 1
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
                if n_cus + 1 < route[j] < N_NODE
                    idx2 = findfirst(x -> x == route[j] - n_cus, route)
                    # RT[j] = A[j] - D[idx2]
                    RT[j] = B[j] - D[idx2]
                end
            end
        end
    end
    #check whether the drive time is respected for all dropoff vertices after i, if so the route is feasible 
    if check_ridetime(route, RT, Li)
        return true
    else
        return false
    end

end




function disp_route(path, route::Vector{Int32}, r, ei, li, s_t, dist, Q, Li, n_cus, TH, qi)

    N = length(route)
    t_start = ei[route[1]]
    A, B, W, D, Load, RT = calculate_ads(route, t_start, ei, li, s_t, dist, n_cus, qi)
    F0 = calc_Fi(1, route, W, B, RT, Li, li)
    Wp = sum(W[2:N-1])
    t_start_new = ei[route[1]] + min(F0, Wp)
    A, B, W, D, Load, RT = calculate_ads(route, t_start_new, ei, li, s_t, dist, n_cus, qi)
    str1 = path * "route_detail_r_" * string(r) * ".xlsx"
    route_detail = DataFrame(route=route, A_i=A, B_i=B, D_i=D, Load_i=Load, RT_i=RT, Li=Li[route], ei=ei[route], li=li[route])
    XLSX.writetable(str1, overwrite=true, route_detail)

end



# check tw constraint
function output_tw_ride_time(path, fast_chg, route, darp, instance, r, ei_route, li_route, s_t_route, dist_orig, qi_route, Li_route)
 
    n_cus, s_t, TH = darp.n_cus, darp.s_t, darp.TH
    set_physical, dist_all = darp.set_physical, darp.dist_all
    Q= instance.Q_max

    N_NODE = 2*n_cus+2
    N=length(route)
    t_start = ei_route[1]
    A, B, W, D, Load, RT = calculate_ads_chg(fast_chg, darp, route, t_start, ei_route, li_route, s_t_route, dist_orig, n_cus, qi_route, set_physical, dist_all)
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
                    v1,v2 = set_physical[route[j-1]-10000], set_physical[[route[j]]] 
                    A[j] = D[j-1] + dist_all[v1, v2] 
                elseif route[j] > N_NODE
                    v1,v2 = set_physical[[route[j-1]]] , set_physical[route[j]-10000] 
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
        @show("check_ridetime_chg failed !")         
        return false, -1,-1,-1,-1,-1,-1
    end

end

# function to check the final solution if the charging operations of the solution do not 
# violate the TW and ride time constraints of customers
function disp_route_chg(path, solution::Solution, darp, instance, r, chg_operation_r, fast_chg::Fast_chargers, ei, li, Li, s_t, dist_orig)

    n_cus, start_depot, end_depot, qi  = darp.n_cus, darp.start_depot, darp.end_depot, darp.qi  
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
    done, A, B, W, D, Load, RT = output_tw_ride_time(path, fast_chg, route, darp, instance, r, ei_route, li_route, s_t_route, dist_orig, qi_route, Li_route)
 

    tw_detail = DataFrame(node = route, t_arr=A, t_begin_service=B, t_wait=W, t_departure=D, ride_time=RT, num_passenger=Load,
                        ei = ei_route, li = li_route, service_or_chg_time = s_t_route, Li=Li_route )
    str1 = path * "route_detail_chg_r_" * string(r) * ".xlsx" 
    XLSX.writetable(str1, overwrite=true, tw_detail)

end

# function disp_route(route::Vector{Int32}, ri, ei, li, s_t, dist, Q, Li, n_cus, TH, nodes)
#     N_NODE = 2*n_cus+2
#     N=length(route)
#     t_start = ei[route[1]]
#     A, B, W, D, Load, RT = calculate_ads(route, t_start, ei, li, s_t, dist, n_cus, qi)
#     F0 = calc_Fi(1, route, W, B, RT, Li, li)
#     Wp = sum(W[2:N-1])
#     t_start_new = ei[route[1]] + min(F0, Wp)
#     A, B, W, D, Load, RT = calculate_ads(route, t_start_new, ei, li, s_t, dist, n_cus, qi)
#     str1 = "route_detail_r_" * string(ri) * ".csv" 

#     if check_ridetime(route, RT, Li) && check_routeduration(route, B, TH) 
#         route_detail = DataFrame(route = route, node_id = nodes[route], A_i = A, B_i = B, W_i = W, D_i = D, Load_i = Load, RT_i = RT, Li = Li[route], ei = ei[route], li = li[route])
#         CSV.write(str1, route_detail)
#         return true 
#     end
#     str2 = "route_detail_r_" * string(ri) * "_before.csv" 
#     route_detail = DataFrame(route = route, node_id = nodes[route], A_i = A, B_i = B, W_i = W, D_i = D, Load_i = Load, RT_i = RT, Li = Li[route], ei = ei[route], li = li[route])
#     CSV.write(str2, route_detail)

#     for i in 2:N-1
#         if  1 < route[i] <= n_cus+1 
#             Fi = calc_Fi(i, route, W, B, RT, Li, li)
#             Wp = sum(W[i+1:N-1])  

#             W[i] += min(Fi, Wp)
#             B[i] = A[i] + W[i]
#             D[i] = B[i] + s_t[route[i]]

#             for j in i+1:N
#                 A[j] = D[j-1] + dist[route[j-1], route[j]]
#                 B[j] = max(A[j], ei[route[j]])
#                 W[j] = B[j] - A[j]
#                 D[j] = B[j] + s_t[route[j]]

#                 if  n_cus+1 < route[j] < N_NODE
#                     idx2 = findfirst(x->x==route[j]-n_cus, route)
#                     RT[j] = B[j] - D[idx2]
#                 end
#             end
#         end
#     end

#     #check whether the drive time is respected for all dropoff vertices after i, if so the route is feasible 
#     if check_routeduration(route, B, TH) && check_ridetime(route, RT, Li)
#         @show("pos 2: ")
#         route_detail = DataFrame(route = route, node_id = nodes[route], A_i = A, B_i = B, W_i = W, D_i = D, Load_i = Load, RT_i = RT, Li = Li[route], ei = ei[route], li = li[route])
#         CSV.write(str1, route_detail)
#         return true 
#     else           
#         if ! check_routeduration(route, B,TH) 
#             @show("check_routeduration(route, B) failed!!")
#             @show(B,B[end] - B[1])
#             @show(A, B, W, D, Load, RT)
#         end
#         ! check_ridetime(route, RT, Li) && @show("check_ridetime(route, RT) failed !!")
#         @show(route, RT,Li,ei[route],li[route])
#         return false
#     end

# end

