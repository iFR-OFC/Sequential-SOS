%% ***********************************************************************
%
%  Region-of-attraction Estimation Benchmark Test Runner (N-link robot arm)
%
%  The benchmark solves the ROA problem for the N-link robot arm for
%  N = Nmin:Nmax links (i.e., 2*N states) with polynomial degree deg. In
%  contrast to the other ROA examples, the problem instances vary in the
%  number of states, not in the polynomial degrees.
%
%  Note 1: The N-link SOS problem differs slightly from the other ROA
%  examples (the fixed level-set constant enters the constraint but not
%  the cost), hence the problem-specific scripts
%  solve_ROA_Nlink_seq_bench_fun.m and solve_ROA_Nlink_CD_bench_fun.m are
%  used; both follow the structure, inputs, and outputs of the generic
%  solve_ROA_seq_bench_fun.m and solve_ROA_CD_bench_fun.m.
%
%  Note 2: The user might select a second hessian approximation method for
%  sequential. This is optional, but interesting in terms of convergence
%  speed. Additionally, a clearly infeasible initial guess for sequential
%  can be enabled.
%
%  License: GNU GENERAL PUBLIC LICENSE Version 3
% ***********************************************************************
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');

close all
clear
clc

% add the problem specific folders and the shared helper functions
thisDir  = fileparts(mfilename('fullpath'));
addPaths = { ...
    fullfile(fileparts(fileparts(thisDir)),'helperFunctions'), ... % appendResultRow, plotConvergence, perf_profile_cdf, evalSequential
    thisDir, ...                                                 % Nlink problem and solver scripts
    fullfile(thisDir,'robot')};                                  % dynamics and LQR data of the N-link robot arm
for ii = 1:numel(addPaths)
    addpath(addPaths{ii})
end

% enable plotting of convergence for each run
enable_plot_convergence = true;
enable_second_sequential = false;

% run sequential additionally with a clearly infeasible initial guess
enable_seq_badInit = false;

% set to true to generate tikz figure using MATLAB2TIKZ
genTikz = false;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization';      % our default (see also sequentialCommon.m)
hess_est_method_2 = 'levenberg';           % 'BFGS', 'levenberg', 'mirroring'

% common options shared by the methods
commonOpts.maxiter       = 100;
commonOpts.tolerance_opt = 1e-4;
commonOpts.maxSolveTime  = 3600;
commonOpts.solver        = 'mosek';
commonOpts.miniter       = 5; % minimum iterations before an aborted CD run still counts as successful

% MOSEK settings


commonOpts.sdp_solver_opts.mosek_param = struct( ...
    'MSK_IPAR_NUM_THREADS',            1);

commonOpts.sdp_solver_opts.mosek_echo = 0;
    % 'MSK_DPAR_INTPNT_CO_TOL_PFEAS',    1e-7, ...
    % 'MSK_DPAR_INTPNT_CO_TOL_DFEAS',    1e-7, ...
    % 'MSK_DPAR_INTPNT_CO_TOL_REL_GAP',  1e-7, ...
    % 'MSK_DPAR_INTPNT_CO_TOL_MU_RED',   1e-7);

% only degree 2 to keep computational effort reasonable for larger number
% of states
deg    = 2;

% number of links: the robot arm has 2*N states
Nmin   = 2;
Nmax   = 7;  % n = 2*N

% variable names of the results table (one row per run)
tableVariableNames = {'Problem Name','no. con. (SDP)','no. dec. (SDP)','Size Gram matrices', ...
    'Build Time (CD)','Build Time (SQ)','Solve Time (CD)','Solve Time (SQ)', ...
    'Iterations (CD)','Iterations (SQ)','Final Cost (CD)','Final Cost (SQ)', ...
    'Solved (CD)','Solved (SQ)'};

% extend the results table by the sequential run with the second Hessian
% estimation method (cf. enable_second_sequential)
if enable_second_sequential
    tableVariableNames = [tableVariableNames, ...
        ['Build Time (SQ ' hess_est_method_2 ')'], ...
        ['Solve Time (SQ ' hess_est_method_2 ')'], ...
        ['Iterations (SQ ' hess_est_method_2 ')'], ...
        ['Final Cost (SQ ' hess_est_method_2 ')'], ...
        ['Solved (SQ ' hess_est_method_2 ')']];
end

% initialize the cell array that accumulates one row per run of the
% results table (one column per table variable)
tableRows = cell(0, numel(tableVariableNames));

% initialize the cell arrays of stored runs (one entry per run)
seqSolverStructs_regularization = {};
seqSolverStructs_hessApprox2 = {};
seqSolverStructs_badInit = {};
cdSolverStructs  = {};
probNames         = {}; % problem name per run (for the evaluation output)

% preallocate (only used if the corresponding run is disabled)
evalSequentialStruct_Nlink_hessApprox2 = [];

for kk = 1:(Nmax-Nmin+1)

    % current number of links
    n = Nmin + kk - 1;

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Nlink = generateNlink_probStruct(n,deg);

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_Nlink.Name, kk, Nmax-Nmin+1)
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_Nlink(probStruct_Nlink);

    % solve Nlink example using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Nlink_regularization = solve_ROA_Nlink_seq_bench_fun(probStruct_Nlink,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_Nlink_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_Nlink_regularization                 = evalSequential(solStruct_seq_Nlink_regularization);
    evalSequentialStruct_Nlink_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Nlink_regularization.probNumber      = kk;


    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Nlink_hessApprox2 = solve_ROA_Nlink_seq_bench_fun(probStruct_Nlink,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Nlink_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_Nlink_hessApprox2 = evalSequential(solStruct_seq_Nlink_hessApprox2);
        evalSequentialStruct_Nlink_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Nlink_hessApprox2.probNumber = kk;
    end

    if enable_seq_badInit
        % solve Nlink example with a clearly infeasible initial guess
        x0_bad = generate_initGuess_Nlink_badInit(probStruct_Nlink);
        solStruct_seq_Nlink_badInit = solve_ROA_Nlink_seq_bench_fun(probStruct_Nlink,commonOpts,x0_bad,hess_est_method_1);

        seqSolverStructs_badInit{end+1} = solStruct_seq_Nlink_badInit;
    end

    % solve Nlink example using custom coordinate-descent
    solStruct_CD_Nlink = solve_ROA_Nlink_CD_bench_fun(probStruct_Nlink,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_Nlink;
    probNames{end+1}       = sprintf('%s (instance %d of %d)', probStruct_Nlink.Name, kk, Nmax-Nmin+1);


    % append the current run to the accumulating results table
    if enable_second_sequential
        tableRows = appendResultRow(tableRows,probStruct_Nlink,solStruct_CD_Nlink, ...
            solStruct_seq_Nlink_regularization,solStruct_seq_Nlink_hessApprox2,true);
    else
        tableRows = appendResultRow(tableRows,probStruct_Nlink,solStruct_CD_Nlink, ...
            solStruct_seq_Nlink_regularization,[],true);
    end

    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_Nlink, ...
                        enable_second_sequential,...
                        evalSequentialStruct_Nlink_regularization,...
                        evalSequentialStruct_Nlink_hessApprox2,...
                        solStruct_CD_Nlink)
    end
end

% display the results table accumulated over all runs
resultsTable = cell2table(tableRows,'VariableNames',tableVariableNames);
disp(resultsTable)

% store the results of all runs for later evaluation
save('Nlink_benchmark_results.mat', ...
    'seqSolverStructs_regularization','seqSolverStructs_hessApprox2', ...
    'seqSolverStructs_badInit','cdSolverStructs','tableRows','tableVariableNames')


%% Evaluation 

runs = Nmin:Nmax;

% extract solve times and iterations of the sequential runs
solverTimes_total_Seq = cellfun(@(s) s.stats.total_solve_time, seqSolverStructs_regularization);
iteration_array_Seq   = cellfun(@(s) s.stats.iterations,      seqSolverStructs_regularization);

% extract solve times and iterations of the coordinate-descent runs
solverTimes_total_CD  = cellfun(@(s) s.stats.totalSolveTime, cdSolverStructs);
iteration_array_CD    = cellfun(@(s) s.stats.iterations,     cdSolverStructs);

% post-processing of coordinate-descent, because some runs are infeasible,
% and shall not appear in the computational statistics; infeasible runs
% have a non-success return status
idx_infes_cd = cellfun(@(s) strcmp(s.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS'), cdSolverStructs);

% extract solve times and iterations of the sequential runs with bad initial guess
if enable_seq_badInit
    solverTimes_total_Seq_badInit = cellfun(@(s) s.stats.total_solve_time, seqSolverStructs_badInit);
    iteration_array_Seq_badInit   = cellfun(@(s) s.stats.iterations,      seqSolverStructs_badInit);
end

% plot solve times
figure('Name', 'Realtive improvement regarding solve times and iterations')
subplot(211)
hold on
rel_solveTime = (solverTimes_total_CD(idx_infes_cd)'-solverTimes_total_Seq(idx_infes_cd)')./solverTimes_total_CD(idx_infes_cd)'*100;
mean_solveTimeRel = mean(rel_solveTime);
yline(mean_solveTimeRel)
hold on
bar(2*runs(idx_infes_cd),rel_solveTime);
% xlabel('Number of states')
legend('Mean value','rel. Improvement','Location','best')
ylabel('Relative improvment in solve time in percent')
xticks((2:Nmax)*2)

% plot iterations (no subiterations of feasibiltity restoration)
subplot(212)
hold on
rel_iteration = (iteration_array_CD(idx_infes_cd)'-iteration_array_Seq(idx_infes_cd)')./iteration_array_CD(idx_infes_cd)'*100;
mean_rel_iteration = mean(rel_iteration);
yline(mean_rel_iteration)
hold on
bar(2*runs(idx_infes_cd),rel_iteration);
xlabel('Number of states')
ylabel('Relative improvment in iterations in percent')
legend('Mean value','rel. Improvement','Location','best')
xticks((2:Nmax)*2)


if genTikz
cleanfigure();
matlab2tikz('NlinkRobot.tex','width','\figW','height','\figH')
end


% plot solve times
figure('Name', 'Solve times and iterations')
subplot(211)
semilogy(runs*2,solverTimes_total_Seq,'-o')
hold on
semilogy(runs(idx_infes_cd)*2,solverTimes_total_CD(idx_infes_cd),'-^')
xlabel('Number of states')
ylabel('Solver time [s]')
legend('sequential','coordinate-descent','Location','northwest')
xticks((2:Nmax)*2)

% plot iterations (no subiterations of feasibiltity restoration)
subplot(212)
 plot(runs*2,iteration_array_Seq,'-o')
hold on
plot(runs(idx_infes_cd)*2,iteration_array_CD(idx_infes_cd),'-^')
xlabel('Number of states')
ylabel('Iterations [-]')
legend('sequential','coordinate-descent','Location','northwest')
xticks((2:Nmax)*2)


% time per iteration
figure('Name', 'Time per iteration')
plot(runs*2,solverTimes_total_Seq./iteration_array_Seq,'-o')
hold on
%plot(runs*2,solverTimes_total_Seq_badInit./iteration_array_Seq_badInit,'-o')
plot(runs(idx_infes_cd)*2,solverTimes_total_CD(idx_infes_cd)./iteration_array_CD(idx_infes_cd),'-^')

xlabel('Number of states')
ylabel('Time/iteration [s]')
legend('sequential','coordinate-descent','Location','northwest')
xticks((2:Nmax)*2)

% get number of linear constraints and decision variables (sequential)
n_con = nan(numel(runs),1);
n_dec = nan(numel(runs),1);
for j = 1:numel(runs)
   n_con(j) = seqSolverStructs_regularization{j}.info.sossol_convex.sdpsol.conic.size_a(1);
   n_dec(j) = seqSolverStructs_regularization{j}.info.sossol_convex.sdpsol.conic.size_a(2);
end


%% Performance Profile (Dolan & More, Math. Program. 91, 2002)

% sequential runs stored per Hessian estimation method (cf. hess_est_method)
seqMethodStructs = {seqSolverStructs_regularization};
seqMethodNames   = {hess_est_method_1};
if enable_second_sequential
    seqMethodStructs{end+1} = seqSolverStructs_hessApprox2;
    seqMethodNames{end+1}   = hess_est_method_2;
end

% solver names: one sequential solver per Hessian estimation method plus coordinate descent
% (explicit loop: arrayfun/concat can produce nested cells, which fprintf rejects)
solverNames = cell(1, numel(seqMethodNames) + 1);
for s = 1:numel(seqMethodNames)
    solverNames{s} = ['Sequential (' seqMethodNames{s} ')'];
end
solverNames{end} = 'Coordinate-descent';

% evaluate the stored test data: shifted geometric means and the
% performance profile plots (legend labels passed as inputs)
plotPerformanceProfiles(seqMethodStructs, cdSolverStructs, solverNames, commonOpts.maxSolveTime, probNames)

% remove the added folders from the path again now that the benchmarks
% are done (exactly the folders added at the top of this script)
for ii = 1:numel(addPaths)
    rmpath(addPaths{ii})
end
