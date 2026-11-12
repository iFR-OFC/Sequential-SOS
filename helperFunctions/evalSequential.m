function  evalSequentialStruct = evalSequential(S)

% get number of entries. Note we also counter the zeroth iteration
nelem = numel(S.stats.single_iterations);

cost                      = nan(nelem,1);
scaled_d_ldx              = nan(nelem,1);
maxConVio                 = nan(nelem,1);
alpha                     = nan(nelem,1);
counter_correction_calls  = nan(nelem,1);
counter_linesearch_calls  = nan(nelem,1);
counter_restoration_calls = nan(nelem,1);

main_sdp_time        = nan(nelem,1);
hess_approx_time     = nan(nelem,1);
total_backstep_time  = nan(nelem,1);
filter_accept_time   = nan(nelem,1);
suff_decrease_time   = nan(nelem,1);
soc_time             = nan(nelem,1);
feas_res_time        = nan(nelem,1);
eval_kkt_con_time    = nan(nelem,1);

for k = 1:nelem
    cost(k)                         = S.stats.single_iterations{k}.cost;
    maxConVio(k)                    = S.stats.single_iterations{k}.max_con_vio;
    alpha(k)                        = S.stats.single_iterations{k}.step_length;
    scaled_d_ldx(k)                 = S.stats.single_iterations{k}.scaled_d_ldx;
    counter_restoration_calls(k)    = S.stats.single_iterations{k}.counter_restoration_calls;
    counter_correction_calls(k)     = S.stats.single_iterations{k}.counter_correction_calls;
    counter_linesearch_calls(k)     = S.stats.single_iterations{k}.counter_linesearch_calls;

    % extract timing statistics
    main_sdp_time(k)        = S.stats.single_iterations{k}.time_stats.main_sdp_time;
    hess_approx_time(k)     = S.stats.single_iterations{k}.time_stats.hess_approx_time;
    total_backstep_time(k)  = S.stats.single_iterations{k}.time_stats.total_backstep_time;
    filter_accept_time(k)   = S.stats.single_iterations{k}.time_stats.filter_accept_time;
    suff_decrease_time(k)   = S.stats.single_iterations{k}.time_stats.suff_decrease_time;
    soc_time(k)             = S.stats.single_iterations{k}.time_stats.soc_time;
    feas_res_time(k)        = S.stats.single_iterations{k}.time_stats.feas_res_time;
    eval_kkt_con_time(k)    = S.stats.single_iterations{k}.time_stats.eval_kkt_con_time;


end


maxConVio(maxConVio == 0) = eps;

evalSequentialStruct.cost                      = cost;
evalSequentialStruct.scaled_d_ldx              = scaled_d_ldx ;
evalSequentialStruct.maxConVio                 = maxConVio;
evalSequentialStruct.alpha                     = alpha;
evalSequentialStruct.time_per_iteration        = S.stats.time_per_iteration;
evalSequentialStruct.total_solve_time          = S.stats.total_solve_time;
evalSequentialStruct.solver_build_time         = S.stats.solver_build_time;
evalSequentialStruct.init_time                 = S.stats.init_time;
evalSequentialStruct.counter_restoration_calls = counter_restoration_calls;
evalSequentialStruct.counter_correction_calls  = counter_correction_calls;
evalSequentialStruct.counter_linesearch_calls  = counter_linesearch_calls;
evalSequentialStruct.iterations                = S.stats.iterations;
evalSequentialStruct.options                   = S.stats.options;


evalSequentialStruct.main_sdp_time             = main_sdp_time;
evalSequentialStruct.hess_approx_time          = hess_approx_time;
evalSequentialStruct.total_backstep_time       = total_backstep_time;
evalSequentialStruct.filter_accept_time        = filter_accept_time;
evalSequentialStruct.suff_decrease_time        = suff_decrease_time;
evalSequentialStruct.soc_time                  = soc_time;
evalSequentialStruct.feas_res_time             = feas_res_time;
evalSequentialStruct.eval_kkt_con_time         =eval_kkt_con_time;

end
 