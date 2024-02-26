 
# please be aware that dist, dist_all, dist_orig are measured in minutes, not in kilometer 
function darp_solver(path, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver) 

    n_run = global_parameter.n_run
    parameter_energy = instance.parameter_energy 
    t_start = time(); avg_k_run = 0; best_k_run=0

    best_solu = solver(path, global_parameter, fast_chg, darp, instance, demand_data, nodes_info, timeTable, output_solver)
    t_elapse = time()- t_start 
    output_solver.t_cpu_run = t_elapse/n_run
    if ! silent_darp_solver
        if n_run > 2 
             vec_tmp = output_solver.obj_values[1:n_run]
             best_k_run, avg_k_run = minimum(vec_tmp), mean(vec_tmp)
        end
        @show(instance.instance_name, n_run) 
        @show(global_parameter.percentage_full, global_parameter.parameter_sa.step_size_remove_route, t_elapse/n_run, output_solver.avg_obj, best_solu.total_cost, avg_k_run , best_k_run )
        #verify_solution(best_solu, instance, darp, fast_chg, global_parameter) 
        n_route = best_solu.n_route ; veh_ids = best_solu.RI[1:n_route,5]
        @show(n_route, veh_ids, parameter_energy.is_electric[veh_ids],
              best_solu.total_cost, best_solu.total_cost_with_penalty, best_solu.total_chg_time, best_solu.unserved_users, best_solu.penalty_unserved, 
              global_parameter.is_gasoline_fleet, global_parameter.co2_reduc_target, global_parameter.co2_threshold, best_solu.co2_emission)
        best_solu.total_chg_time > 0 &&  @show(best_solu.vec_chg_event[1:n_route,:])         
        export_route_detail(path, best_solu , instance, darp, fast_chg) 
        # for r in 1:best_solu.n_route
        #     route = get_route(best_solu, r, darp.start_depot,darp.end_depot )
        #     check_energy_simple(best_solu, darp, route, r, parameter_energy, fast_chg )
        # end
    end
    
    # n_route = best_solu.n_route
    # veh_idx_used = best_solu.RI[1:n_route,5]
    # for r in 1:best_solu.n_route  
    #     route = get_route(best_solu, r, darp.start_depot, darp.end_depot)
    #     @show(r, route)
    # end
    # show_co2_detail(best_solu, darp, instance, parameter_energy)
    # @show(verify_tw_chg_operation(best_solu, fast_chg, darp, instance))
    # @show(check_chg_occ_constr(best_solu, max_n_fast_chg, instance))
    # export_route_detail(path, best_solu , instance, darp, fast_chg)
    # @show(best_solu.penalty_unserved, best_solu.n_route, veh_idx_used, parameter_energy.is_electric[veh_idx_used], best_solu.total_cost, best_solu.vec_chg_event[1:n_route,:] )
    # @show(n_run, t_elapse/n_run, output_solver.avg_obj) 
  
    return best_solu, output_solver
end

# please be aware that dist, dist_all, dist_orig are measured in minutes, not in kilometer 
function darp_solver_batch(path_result, silent_darp_solver, global_parameter, fast_chg, darp, instance, nodes_info, demand_data, output_solver) 

    n_run = global_parameter.n_run
    parameter_energy = instance.parameter_energy 
    t_start = time(); avg_k_run = 0; best_k_run=0
    best_solu = solver(path_result, global_parameter, fast_chg, darp, instance, demand_data, nodes_info, timeTable, output_solver)
    t_elapse = time()- t_start 
    output_solver.t_cpu_run = t_elapse/n_run
    if ! silent_darp_solver
        if n_run > 2 
             vec_tmp = output_solver.obj_values[1:n_run]
             best_k_run, avg_k_run = minimum(vec_tmp), mean(vec_tmp)
        end
        n_route = best_solu.n_route ; veh_ids = best_solu.RI[1:n_route,5]
        @show(instance.instance_name, n_run, global_parameter.percentage_full,  t_elapse/n_run, 
              best_solu.total_cost, best_solu.unserved_users, veh_ids, parameter_energy.is_electric[veh_ids],
              global_parameter.co2_reduc_target, global_parameter.co2_threshold, best_solu.co2_emission)
        verify_solution(best_solu, instance, darp, fast_chg, global_parameter) 

        # @show(n_route, veh_ids, parameter_energy.is_electric[veh_ids],
        #       best_solu.total_cost, best_solu.total_cost_with_penalty, best_solu.total_chg_time, best_solu.unserved_users, 
        #       global_parameter.is_gasoline_fleet, global_parameter.co2_reduc_target, global_parameter.co2_threshold, best_solu.co2_emission)
        # best_solu.total_chg_time > 0 &&  @show(best_solu.vec_chg_event[1:n_route,:])         
        # export_route_detail(path_result, best_solu , instance, darp, fast_chg) 
        # total_cost_best_sol,  cost_with_penalty_best_sol = cost_solution_check(best_solu, global_parameter, darp, parameter_energy, fast_chg, instance)
        # @show(total_cost_best_sol,  cost_with_penalty_best_sol)
    end
    
    # n_route = best_solu.n_route
    # veh_idx_used = best_solu.RI[1:n_route,5]
    # for r in 1:best_solu.n_route  
    #     route = get_route(best_solu, r, darp.start_depot, darp.end_depot)
    #     @show(r, route)
    # end
    # show_co2_detail(best_solu, darp, instance, parameter_energy)
    # @show(verify_tw_chg_operation(best_solu, fast_chg, darp, instance))
    # @show(check_chg_occ_constr(best_solu, max_n_fast_chg, instance))
    # export_route_detail(path, best_solu , instance, darp, fast_chg)
    # @show(best_solu.penalty_unserved, best_solu.n_route, veh_idx_used, parameter_energy.is_electric[veh_idx_used], best_solu.total_cost, best_solu.vec_chg_event[1:n_route,:] )
    # @show(n_run, t_elapse/n_run, output_solver.avg_obj) 
  
    return best_solu, output_solver
end