%% -----------------------------------------------------------------------
%
%  Short Description: This script generates the evaluation for desired
%                     benchmarks. It reads in the .mat file from the
%                     desired benchmark problem. The user can decide which
%                     plots are generated
%
%% -----------------------------------------------------------------------

clc
close all
clear

% load desired benchmark
load('ROA_benchmark.mat')

%% get detailed breakdown of computational effort

seq_nRuns = numel(seqSolverStructs_regularization);

% initialize arrays
seq_main_sdp_time_arr        = nan(seq_nRuns, 1);
seq_feas_res_time_arr        = nan(seq_nRuns, 1);
seq_filter_accept_time_arr   = nan(seq_nRuns, 1);
seq_hess_approx_time_arr     = nan(seq_nRuns, 1);
seq_soc_time_arr             = nan(seq_nRuns, 1);
seq_suff_decrease_time_arr   = nan(seq_nRuns, 1);
seq_total_backstep_time_arr  = nan(seq_nRuns, 1);

seq_total_solve_time_arr     = nan(seq_nRuns, 1);
seq_solver_build_time_arr    = nan(seq_nRuns, 1);
seq_solver_init_time_arr     = nan(seq_nRuns, 1);

seq_rel_time_main_sdp            = nan(seq_nRuns, 1);
seq_rel_time_feas_res            = nan(seq_nRuns, 1);
seq_rel_time_filter_accept_time   = nan(seq_nRuns, 1);
seq_rel_time_hess_approx         = nan(seq_nRuns, 1);
seq_rel_time_soc                 = nan(seq_nRuns, 1);
seq_rel_time_suff_decrease_time  = nan(seq_nRuns, 1);
seq_rel_time_total_backstep      = nan(seq_nRuns, 1);
seq_rel_time_init                = nan(seq_nRuns, 1);

seq_eval_kkt_con_time_arr        = nan(seq_nRuns, 1);
seq_rel_time_eval_kkt_con        = nan(seq_nRuns, 1);
seq_rel_time_overhead            = nan(seq_nRuns, 1);

% backstep-only breakdown
seq_rel_time_filter_accept_backstep = nan(seq_nRuns, 1);
seq_rel_time_soc_backstep           = nan(seq_nRuns, 1);
seq_rel_time_suff_decrease_backstep = nan(seq_nRuns, 1);
seq_rel_time_backstep_other         = nan(seq_nRuns, 1);

% coordinate descent arrays
cd_totalSolveTime_arr      = nan(seq_nRuns, 1);
cd_totalSolveTime_S1_arr   = nan(seq_nRuns, 1);
cd_totalSolveTime_S2_arr   = nan(seq_nRuns, 1);
cd_init_time_arr           = nan(seq_nRuns, 1);
cd_eval_kkt_time_arr       = nan(seq_nRuns, 1);

cd_rel_totalSolveTime_S1_arr = nan(seq_nRuns, 1);
cd_rel_totalSolveTime_S2_arr = nan(seq_nRuns, 1);
cd_rel_init_time_arr         = nan(seq_nRuns, 1);
cd_rel_eval_kkt_time_arr     = nan(seq_nRuns, 1);

for seq_kk = 1:seq_nRuns

    evalSequentialStruct_regularization = evalSequential(seqSolverStructs_regularization{seq_kk});

    % extract individual times
    seq_main_sdp_time_arr(seq_kk)       = sum(evalSequentialStruct_regularization.main_sdp_time);
    seq_feas_res_time_arr(seq_kk)       = sum(evalSequentialStruct_regularization.feas_res_time);
    seq_filter_accept_time_arr(seq_kk)  = sum(evalSequentialStruct_regularization.filter_accept_time);
    seq_hess_approx_time_arr(seq_kk)    = sum(evalSequentialStruct_regularization.hess_approx_time);
    seq_soc_time_arr(seq_kk)            = sum(evalSequentialStruct_regularization.soc_time);
    seq_suff_decrease_time_arr(seq_kk)  = sum(evalSequentialStruct_regularization.suff_decrease_time);
    seq_total_backstep_time_arr(seq_kk) = sum(evalSequentialStruct_regularization.total_backstep_time);
    seq_eval_kkt_con_time_arr(seq_kk)   = sum(evalSequentialStruct_regularization.eval_kkt_con_time);

    % total solve and build times
    seq_total_solve_time_arr(seq_kk)    = evalSequentialStruct_regularization.total_solve_time;
    seq_solver_build_time_arr(seq_kk)    = evalSequentialStruct_regularization.solver_build_time;
    seq_solver_init_time_arr(seq_kk)     = evalSequentialStruct_regularization.init_time;

    % relative times w.r.t. total solve time
    if seq_total_solve_time_arr(seq_kk) > 0
        seq_rel_time_main_sdp(seq_kk)           = seq_main_sdp_time_arr(seq_kk) / seq_total_solve_time_arr(seq_kk);
        seq_rel_time_feas_res(seq_kk)           = seq_feas_res_time_arr(seq_kk) / seq_total_solve_time_arr(seq_kk);
        seq_rel_time_filter_accept_time(seq_kk) = seq_filter_accept_time_arr(seq_kk) / seq_total_solve_time_arr(seq_kk);
        seq_rel_time_hess_approx(seq_kk)        = seq_hess_approx_time_arr(seq_kk) / seq_total_solve_time_arr(seq_kk);
        seq_rel_time_soc(seq_kk)                = seq_soc_time_arr(seq_kk) / seq_total_solve_time_arr(seq_kk);
        seq_rel_time_suff_decrease_time(seq_kk)  = seq_suff_decrease_time_arr(seq_kk) / seq_total_solve_time_arr(seq_kk);
        seq_rel_time_total_backstep(seq_kk)      = seq_total_backstep_time_arr(seq_kk) / seq_total_solve_time_arr(seq_kk);
        seq_rel_time_init(seq_kk)                = seq_solver_init_time_arr(seq_kk) / seq_total_solve_time_arr(seq_kk);
        seq_rel_time_eval_kkt_con(seq_kk)        = seq_eval_kkt_con_time_arr(seq_kk) / seq_total_solve_time_arr(seq_kk);
    else
        seq_rel_time_main_sdp(seq_kk)           = 0;
        seq_rel_time_feas_res(seq_kk)           = 0;
        seq_rel_time_filter_accept_time(seq_kk) = 0;
        seq_rel_time_hess_approx(seq_kk)        = 0;
        seq_rel_time_soc(seq_kk)                = 0;
        seq_rel_time_suff_decrease_time(seq_kk)  = 0;
        seq_rel_time_total_backstep(seq_kk)      = 0;
        seq_rel_time_init(seq_kk)                = 0;
        seq_rel_time_eval_kkt_con(seq_kk)        = 0;
    end

    seq_rel_time_overhead(seq_kk) = max(0, 1 - (seq_rel_time_main_sdp(seq_kk) + ...
        seq_rel_time_feas_res(seq_kk) + seq_rel_time_total_backstep(seq_kk) + ...
        seq_rel_time_hess_approx(seq_kk) + seq_rel_time_init(seq_kk) + ...
        seq_rel_time_eval_kkt_con(seq_kk)));

    % relative times w.r.t. total backstep time
    if seq_total_backstep_time_arr(seq_kk) > 0
        seq_rel_time_filter_accept_backstep(seq_kk) = seq_filter_accept_time_arr(seq_kk) / seq_total_backstep_time_arr(seq_kk);
        seq_rel_time_soc_backstep(seq_kk)           = seq_soc_time_arr(seq_kk) / seq_total_backstep_time_arr(seq_kk);
        seq_rel_time_suff_decrease_backstep(seq_kk)  = seq_suff_decrease_time_arr(seq_kk) / seq_total_backstep_time_arr(seq_kk);
        seq_rel_time_backstep_other(seq_kk)          = max(0, 1 - (seq_rel_time_filter_accept_backstep(seq_kk) + ...
                                                                   seq_rel_time_soc_backstep(seq_kk) + ...
                                                                   seq_rel_time_suff_decrease_backstep(seq_kk)));

        seq_backstepFractions = [seq_rel_time_filter_accept_backstep(seq_kk), ...
                                 seq_rel_time_soc_backstep(seq_kk), ...
                                 seq_rel_time_suff_decrease_backstep(seq_kk), ...
                                 seq_rel_time_backstep_other(seq_kk)];
        seq_backstepLabels = {
            sprintf('Filter accept\n%.1f\\%%', 100 * seq_rel_time_filter_accept_backstep(seq_kk))
            sprintf('SOC\n%.1f\\%%', 100 * seq_rel_time_soc_backstep(seq_kk))
            sprintf('Suff. decrease\n%.1f\\%%', 100 * seq_rel_time_suff_decrease_backstep(seq_kk))
            sprintf('Other backstep\n%.1f\\%%', 100 * seq_rel_time_backstep_other(seq_kk))
        };
    else
        seq_rel_time_filter_accept_backstep(seq_kk) = 0;
        seq_rel_time_soc_backstep(seq_kk)           = 0;
        seq_rel_time_suff_decrease_backstep(seq_kk) = 0;
        seq_rel_time_backstep_other(seq_kk)         = 1;
        seq_backstepFractions = 1;
        seq_backstepLabels = {sprintf('No backstep time\n100.0\\%%')};
    end

    % first pie: total solve time
    if seq_rel_time_feas_res(seq_kk) > 0
        seq_pieFractions = [seq_rel_time_main_sdp(seq_kk), ...
                            seq_rel_time_feas_res(seq_kk), ...
                            seq_rel_time_hess_approx(seq_kk), ...
                            seq_rel_time_total_backstep(seq_kk), ...
                            seq_rel_time_init(seq_kk), ...
                            seq_rel_time_eval_kkt_con(seq_kk)];
        seq_pieLabels = {
            sprintf('Main SDP\n%.1f\\%%', 100 * seq_rel_time_main_sdp(seq_kk))
            sprintf('Feas. res.\n%.1f\\%%', 100 * seq_rel_time_feas_res(seq_kk))
            sprintf('Hess. approx.\n%.1f\\%%', 100 * seq_rel_time_hess_approx(seq_kk))
            sprintf('Backstep\n%.1f\\%%', 100 * seq_rel_time_total_backstep(seq_kk))
            sprintf('Init\n%.1f\\%%', 100 * seq_rel_time_init(seq_kk))
            sprintf('Eval KKT\n%.1f\\%%', 100 * seq_rel_time_eval_kkt_con(seq_kk))
        };
    else
        seq_pieFractions = [seq_rel_time_main_sdp(seq_kk), ...
                            seq_rel_time_hess_approx(seq_kk), ...
                            seq_rel_time_total_backstep(seq_kk), ...
                            seq_rel_time_init(seq_kk), ...
                            seq_rel_time_eval_kkt_con(seq_kk)];
        seq_pieLabels = {
            sprintf('Main SDP\n%.1f\\%%', 100 * seq_rel_time_main_sdp(seq_kk))
            sprintf('Hess. approx.\n%.1f\\%%', 100 * seq_rel_time_hess_approx(seq_kk))
            sprintf('Backstep\n%.1f\\%%', 100 * seq_rel_time_total_backstep(seq_kk))
            sprintf('Init\n%.1f\\%%', 100 * seq_rel_time_init(seq_kk))
            sprintf('Eval KKT\n%.1f\\%%', 100 * seq_rel_time_eval_kkt_con(seq_kk))
        };
    end

    seq_rel_time_overhead(seq_kk) = max(0, 1 - sum(seq_pieFractions));
    seq_pieFractions(end + 1) = seq_rel_time_overhead(seq_kk);
    seq_pieLabels{end + 1} = sprintf('Overhead\n%.1f\\%%', 100 * seq_rel_time_overhead(seq_kk));

        figure('Name', probNames{seq_kk});
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax1 = nexttile;
    pie(ax1, seq_pieFractions, seq_pieLabels);
    t1 = title(ax1, 'Total solve time');
    t1.Units = 'normalized';
    t1.Position(2) = 1.12;

    ax2 = nexttile;
    pie(ax2, seq_backstepFractions, seq_backstepLabels);
    t2 = title(ax2, 'Total backstep time');
    t2.Units = 'normalized';
    t2.Position(2) = 1.12;


    % sgtitle(probNames{seq_kk});

    %% evaluate coordinate descent
    cd_totalSolveTime_arr(seq_kk)    = cdSolverStructs{seq_kk}.stats.totalSolveTime;
    cd_totalSolveTime_S1_arr(seq_kk) = cdSolverStructs{seq_kk}.stats.totalSolveTime_S1;
    cd_totalSolveTime_S2_arr(seq_kk) = cdSolverStructs{seq_kk}.stats.totalSolveTime_S2;
    cd_init_time_arr(seq_kk)        = cdSolverStructs{seq_kk}.stats.init_time;
    cd_eval_kkt_time_arr(seq_kk)    = cdSolverStructs{seq_kk}.stats.eval_kkt_time;

    if cd_totalSolveTime_arr(seq_kk) > 0
        cd_rel_totalSolveTime_S1_arr(seq_kk) = cd_totalSolveTime_S1_arr(seq_kk) / cd_totalSolveTime_arr(seq_kk);
        cd_rel_totalSolveTime_S2_arr(seq_kk) = cd_totalSolveTime_S2_arr(seq_kk) / cd_totalSolveTime_arr(seq_kk);
        cd_rel_init_time_arr(seq_kk)         = cd_init_time_arr(seq_kk) / cd_totalSolveTime_arr(seq_kk);
        cd_rel_eval_kkt_time_arr(seq_kk)     = cd_eval_kkt_time_arr(seq_kk) / cd_totalSolveTime_arr(seq_kk);
    else
        cd_rel_totalSolveTime_S1_arr(seq_kk) = 0;
        cd_rel_totalSolveTime_S2_arr(seq_kk) = 0;
        cd_rel_init_time_arr(seq_kk)         = 0;
        cd_rel_eval_kkt_time_arr(seq_kk)     = 0;
    end

end

seq_meanNames = [
    "Main SDP"
    "Feasibility Restoration"
    "Filter acceptance"
    "Hessian approximation"
    "Second-Order-Correction"
    "Sufficient decrease"
    "Backstepping"
    "Initialization"
    "Evaluation KKT Condition"
    "Overhead"
    ];

seq_meanValues = 100 * [
    mean(seq_rel_time_main_sdp, 'omitnan')
    mean(seq_rel_time_feas_res, 'omitnan')
    mean(seq_rel_time_filter_accept_time, 'omitnan')
    mean(seq_rel_time_hess_approx, 'omitnan')
    mean(seq_rel_time_soc, 'omitnan')
    mean(seq_rel_time_suff_decrease_time, 'omitnan')
    mean(seq_rel_time_total_backstep, 'omitnan')
    mean(seq_rel_time_init, 'omitnan')
    mean(seq_rel_time_eval_kkt_con, 'omitnan')
    mean(seq_rel_time_overhead, 'omitnan')
    ];

seq_meanTable = table(seq_meanNames, seq_meanValues, ...
    'VariableNames', {'Category', 'MeanRelativeTimePercent'});

seq_meanTable = sortrows(seq_meanTable, 'MeanRelativeTimePercent', 'descend');
disp(seq_meanTable);

cd_meanNames = [
    "Solve time S1"
    "Solve time S2"
    "Initialization"
    "Evaluation KKT time"
    ];

cd_meanValues = 100 * [
    mean(cd_rel_totalSolveTime_S1_arr, 'omitnan')
    mean(cd_rel_totalSolveTime_S2_arr, 'omitnan')
    mean(cd_rel_init_time_arr, 'omitnan')
    mean(cd_rel_eval_kkt_time_arr, 'omitnan')
    ];

cd_meanTable = table(cd_meanNames, cd_meanValues, ...
    'VariableNames', {'Category', 'MeanRelativeTimePercent'});

cd_meanTable = sortrows(cd_meanTable, 'MeanRelativeTimePercent', 'descend');
disp(cd_meanTable);
