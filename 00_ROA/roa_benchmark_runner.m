%% ***********************************************************************
%
%  Region-of-attraction Estimation Benchmark Test Runner
%
%  Benchmark 1: Van-der-Pol Oscillator
%  Benchmark 2: Multi-Input-Multi-Output
%  Benchmark 3: Double Pendulum
%  Benchmark 4: Aircraft
%  Benchmark 5: F-18
%  Benchmark 6: Network
%
%  Note 1: Each runs uses various combinations of polynomial degrees. For SOS
%  problem to be feasible or to improve solution quality higher degree are
%  needed. Hence, with a higher degree more decision variables occur and
%  the Gram matrices become larger. Therefore, we can check the behavior
%  for over several benchmark problems, i.e., we have sub-benchmarks so to
%  say that increase our benchmark set.
%
%  Note 2: The user might select a second hessian approximation method for
%  sequential. This is optional, but interesting in terms of convergence
%  speed.
%  
%
%
% ***********************************************************************
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');

close all
clear
clc

% add the problem specific folders
addpath("F18/")
addpath("VDP/")
addpath("MultiInput/")
addpath("DoublePendulum/")
addpath("Aircraft/")
addpath("Network/")

% enable plotting of convergence for each run
enable_plot_convergence = true;
enable_second_sequential = true;

% Hessian estimation methods of the sequential solver to compare
hess_est_method_1 = 'regularization'; % our default (see also sequentialCommon.m)
hess_est_method_2 = 'levenberg';           % 'BFGS', 'levenberg', 'mirroring'

% common options shared by the methods
commonOpts.maxiter       = 350;
commonOpts.tolerance_opt = 1e-4;
commonOpts.maxSolveTime  = 3600;
commonOpts.solver        = 'mosek';


% initialize the cell array that accumulates one row per run of the results table
tableRows = cell(0, 12);

% variable names of the results table (one row per run)
tableVariableNames = {'Problem Name','no. con. (SDP)','no. dec. (SDP)','Size Gram matrices', ...
    'Build Time (CD)','Build Time (SQ)','Solve Time (CD)','Solve Time (SQ)', ...
    'Iterations (CD)','Iterations (SQ)','Final Cost (CD)','Final Cost (SQ)'};


%% VDP

% define the problem instances based on different polynomial degrees
degV = {2,2:4,1:4,2:6};
degS = {2,2,2:4,2:4};

% initialize the cell arrays of stored runs (one entry per run)
seqSolverStructs_regularization = {};
seqSolverStructs_hessApprox2 = {};
cdSolverStructs  = {};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_VDP = generateVDP_probStruct(degV{kk},degS{kk});
    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_VDP.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_VDP(probStruct_VDP);

    % solve VDP example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_VDP_regularization = solve_ROA_seq_bench_fun(probStruct_VDP,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_VDP_regularization;

    % extract size and timing statistic from sequential and add some meta
    % data for plotting
    evalSequentialStruct_VDP_regularization                 = evalSequential(solStruct_seq_VDP_regularization);
    evalSequentialStruct_VDP_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_VDP_regularization.probNumber      = kk;
    

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_VDP_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_VDP,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_VDP_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_VDP_hessApprox2 = evalSequential(solStruct_seq_VDP_hessApprox2);
        evalSequentialStruct_VDP_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_VDP_hessApprox2.probNumber = kk;
    end

    % solve VDP example with two states using custom coordinate-descent
    solStruct_CD_VDP = solve_ROA_CD_bench_fun(probStruct_VDP,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_VDP;


    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows,probStruct_VDP,solStruct_CD_VDP,solStruct_seq_VDP_regularization);

    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_VDP, ...
                        enable_second_sequential,...
                        evalSequentialStruct_VDP_regularization,...
                        evalSequentialStruct_VDP_hessApprox2,...
                        solStruct_CD_VDP)
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Multi-Input-Multi-Output

% define the problem instances based on different polynomial degrees
degV = {2, 2:4, 2:4, 2:6};
degS = {2,   2, 2:4, 2:4};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_MiMo = generateMiMo_probStruct(degV{kk},degS{kk});
    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_MiMo.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_MiMo(probStruct_MiMo);

    % solve MiMo example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_MiMo_regularization = solve_ROA_seq_bench_fun(probStruct_MiMo,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_MiMo_regularization;

    % extract size and timing statistic from sequential
    evalSequentialStruct_MiMo_regularization = evalSequential(solStruct_seq_MiMo_regularization);
    evalSequentialStruct_MiMo_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_MiMo_regularization.probNumber = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_MiMo_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_MiMo,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_MiMo_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_MiMo_hessApprox2 = evalSequential(solStruct_seq_MiMo_hessApprox2);
        evalSequentialStruct_MiMo_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_MiMo_hessApprox2.probNumber = kk;
    end

    % solve MiMo example with two states using custom coordinate-descent
    solStruct_CD_MiMo = solve_ROA_CD_bench_fun(probStruct_MiMo,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_MiMo;


    % Generate table

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows,probStruct_MiMo,solStruct_CD_MiMo,solStruct_seq_MiMo_regularization);

    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_MiMo, ...
                        enable_second_sequential,...
                        evalSequentialStruct_MiMo_regularization,...
                        evalSequentialStruct_MiMo_hessApprox2,...
                        solStruct_CD_MiMo)
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Double Pendulum

% define the problem instances based on different polynomial degrees
degV = {2;2:4;2:4};
degS = {2;2;2:4};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_DoublePendulum = generateDoublePendulum_probStruct(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_DoublePendulum.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_DoublePendulum(probStruct_DoublePendulum);

    % solve DoublePendulum example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_DoublePendulum_regularization = solve_ROA_seq_bench_fun(probStruct_DoublePendulum,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_DoublePendulum_regularization;

    % extract size and timing statistic from sequential
    evalSequentialStruct_DoublePendulum_regularization = evalSequential(solStruct_seq_DoublePendulum_regularization);
    evalSequentialStruct_DoublePendulum_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_DoublePendulum_regularization.probNumber = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_DoublePendulum_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_DoublePendulum,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_DoublePendulum_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_DoublePendulum_hessApprox2 = evalSequential(solStruct_seq_DoublePendulum_hessApprox2);
        evalSequentialStruct_DoublePendulum_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_DoublePendulum_hessApprox2.probNumber = kk;
    end

    % solve DoublePendulum example with two states using custom coordinate-descent
    solStruct_CD_DoublePendulum = solve_ROA_CD_bench_fun(probStruct_DoublePendulum,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_DoublePendulum;

    %% Generate table

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows,probStruct_DoublePendulum,solStruct_CD_DoublePendulum,solStruct_seq_DoublePendulum_regularization);


    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_DoublePendulum, ...
                        enable_second_sequential,...
                        evalSequentialStruct_DoublePendulum_regularization,...
                        evalSequentialStruct_DoublePendulum_hessApprox2,...
                        solStruct_CD_DoublePendulum)
    end
end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


%% Aircraft 4D


% define the problem instances based on different polynomial degrees
degV = {2,2:4,1:4};
degS = {2:4,2:4,2:4};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_Aircraft = generateAircraft_probStruct(degV{kk},degS{kk});
    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_Aircraft.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_Aircraft(probStruct_Aircraft);

    % solve Aircraft example with two states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_Aircraft_regularization = solve_ROA_seq_bench_fun(probStruct_Aircraft,commonOpts,x0,hess_est_method);

    seqSolverStructs_regularization{end+1} = solStruct_seq_Aircraft_regularization;

    % extract size and timing statistic from sequential
    evalSequentialStruct_Aircraft_regularization = evalSequential(solStruct_seq_Aircraft_regularization);
    evalSequentialStruct_Aircraft_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_Aircraft_regularization.probNumber = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_Aircraft_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_Aircraft,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_Aircraft_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_Aircraft_hessApprox2 = evalSequential(solStruct_seq_Aircraft_hessApprox2);
        evalSequentialStruct_Aircraft_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_Aircraft_hessApprox2.probNumber = kk;
    end

    % solve Aircraft example with two states using custom coordinate-descent
    solStruct_CD_Aircraft = solve_ROA_CD_bench_fun(probStruct_Aircraft,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_Aircraft;

    %% Generate table

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows,probStruct_Aircraft,solStruct_CD_Aircraft,solStruct_seq_Aircraft_regularization);

    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_Aircraft, ...
                        enable_second_sequential,...
                        evalSequentialStruct_Aircraft_regularization,...
                        evalSequentialStruct_Aircraft_hessApprox2,...
                        solStruct_CD_Aircraft)
    end

end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))

%% F/A-18

% define the problem instances based on different polynomial degrees
degV = {2, 2:4, 2:4};
degS = {2,   2, 2:4};

for kk = 1:numel(degV)

    % generate problem structure (dynamics, polynomial dec. variables)
    probStruct_F18 = generateF18_probStruct(degV{kk},degS{kk});

    fprintf('------------------------------------------------------------------------------------------ \n')
    fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_F18.Name, kk, numel(degV))
    fprintf('------------------------------------------------------------------------------------------\n')

    % generate an initial guess
    x0 = generate_initGuess_F18(probStruct_F18);

    % solve F-18 example with seven states using sequential SOS
    hess_est_method = hess_est_method_1;
    solStruct_seq_7D_regularization = solve_ROA_seq_bench_fun(probStruct_F18,commonOpts,x0,hess_est_method);

    % store in a common struct for all sequential instances and problem
    seqSolverStructs_regularization{end+1} = solStruct_seq_7D_regularization;

    % extract size and timing statistic from sequential
    evalSequentialStruct_7D_regularization = evalSequential(solStruct_seq_7D_regularization);
    evalSequentialStruct_7D_regularization.hess_est_method = hess_est_method_1;
    evalSequentialStruct_7D_regularization.probNumber = kk;

    if enable_second_sequential

        hess_est_method = hess_est_method_2;
        solStruct_seq_7D_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_F18,commonOpts,x0,hess_est_method);

        seqSolverStructs_hessApprox2{end+1} = solStruct_seq_7D_hessApprox2;

        % extract size and timing statistic from sequential
        evalSequentialStruct_7D_hessApprox2 = evalSequential(solStruct_seq_7D_hessApprox2);
        evalSequentialStruct_7D_hessApprox2.hess_est_method = hess_est_method_2;
        evalSequentialStruct_7D_hessApprox2.probNumber = kk;
    end

    % solve F-18 example with seven states using custom coordinate-descent
    solStruct_CD_7D = solve_ROA_CD_bench_fun(probStruct_F18,commonOpts,x0);

    cdSolverStructs{end+1} = solStruct_CD_7D;

    %% Generate table

    % append the current run to the accumulating results table
    tableRows = appendResultRow(tableRows,probStruct_F18,solStruct_CD_7D,solStruct_seq_7D_regularization);


    if enable_plot_convergence
        % Evaluate
        plotConvergence(probStruct_F18, ...
                        enable_second_sequential,...
                        evalSequentialStruct_7D_regularization,...
                        evalSequentialStruct_7D_hessApprox2,...
                        solStruct_CD_7D)
    end

end

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))

%% 12-State Dynamical System

% generate problem structure (dynamics, polynomial dec. variables)
probStruct_12State= generate12State_probStruct();
fprintf('------------------------------------------------------------------------------------------ \n')
fprintf('Region-of-Attraction Estimation Benchmark for %s Instance %d of %d. \n', probStruct_12State.Name, 1, 1)
fprintf('------------------------------------------------------------------------------------------\n')

% generate an initial guess
x0 = generate_initGuess_12State(probStruct_12State);

% solve 12 state dynamical system example with 12 states using sequential SOS
hess_est_method = hess_est_method_1;
solStruct_seq_12D_regularization = solve_ROA_seq_bench_fun(probStruct_12State,commonOpts,x0,hess_est_method);

seqSolverStructs_regularization{end+1} = solStruct_seq_12D_regularization;

% extract size and timing statistic from sequential
evalSequentialStruct_12D_regularization = evalSequential(solStruct_seq_12D_regularization);
evalSequentialStruct_12D_regularization.hess_est_method = hess_est_method_1;
evalSequentialStruct_12D_regularization.probNumber = 1;

if enable_second_sequential

    hess_est_method = hess_est_method_2;
    solStruct_seq_12D_hessApprox2 = solve_ROA_seq_bench_fun(probStruct_12State,commonOpts,x0,hess_est_method);

    seqSolverStructs_hessApprox2{end+1} = solStruct_seq_12D_hessApprox2;

    % extract size and timing statistic from sequential
    evalSequentialStruct_12D_hessApprox2 = evalSequential(solStruct_seq_12D_hessApprox2);
    evalSequentialStruct_12D_hessApprox2.hess_est_method = hess_est_method_2;
    evalSequentialStruct_12D_hessApprox2.probNumber = 1;
end

% solve 12 state dynamical system example with 12 states using custom coordinate-descent
solStruct_CD_12D = solve_ROA_CD_bench_fun(probStruct_12State,commonOpts,x0);

cdSolverStructs{end+1} = solStruct_CD_12D;

%% Generate table

% append the current run to the accumulating results table
tableRows = appendResultRow(tableRows,probStruct_12State,solStruct_CD_12D,solStruct_seq_12D_regularization);

% display the results table accumulated over all runs so far
disp(cell2table(tableRows,'VariableNames',tableVariableNames))


if enable_plot_convergence
    % Evaluate
    plotConvergence(probStruct_12State, ...
                    enable_second_sequential,...
                    evalSequentialStruct_12D_regularization,...
                    evalSequentialStruct_12D_hessApprox2,...
                    solStruct_CD_12D)
end


%% Results table (all benchmarks, all runs)

% final results table accumulated over all runs and all problems
resultsTable = cell2table(tableRows,'VariableNames',tableVariableNames);
disp(resultsTable)


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

% number of problems n_p: one stored run per cell entry
n_p = numel(cdSolverStructs);

% number of solvers n_s: one sequential solver per Hessian estimation method and coordinate descent
n_s = numel(solverNames);

% computing time t_{p,s}: one row per problem p in P, one column per solver s in S
% extract each run from the cell arrays of stored solver structs
t_ps = zeros(n_p, n_s);
for p = 1:n_p
    % sequential runs of problem p, one per Hessian estimation method
    for s = 1:numel(seqMethodStructs)
        solStruct_seq = seqMethodStructs{s}{p};
        t_ps(p, s) = solStruct_seq.stats.total_solve_time;
    end
    % coordinate-descent run of problem p
    solStruct_cd  = cdSolverStructs{p};
    % store the computing times t_{p,s} of all solvers
    t_ps(p, n_s) = solStruct_cd.stats.totalSolveTime;
end

% logical matrix: solved_ps(p,s) = true iff solver s solves problem p
% a run only counts as solved for the return status 'SOLVER_RET_SUCCESS'
solved_ps = false(n_p, n_s);
for p = 1:n_p
    % sequential runs of problem p, one per Hessian estimation method
    for s = 1:numel(seqMethodStructs)
        solStruct_seq = seqMethodStructs{s}{p};
        solved_ps(p, s) = strcmp(solStruct_seq.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS');
    end
    % coordinate-descent run of problem p
    solStruct_cd  = cdSolverStructs{p};
    % solved_ps(p,n_s): coordinate descent
    solved_ps(p, n_s) = strcmp(solStruct_cd.stats.UNIFIED_RETURN_STATUS,  'SOLVER_RET_SUCCESS');
end

% drop problems that no solver could solve (cf. drop_unsolvable! in
% ClarabelBenchmarks): all-solvers-failed runs must not distort the
% profile denominators nor the shifted geometric means
keep_p = any(solved_ps, 2);
t_ps(~keep_p, :)      = [];
solved_ps(~keep_p, :) = [];

% number of problems n_p remaining in the test set
n_p = size(t_ps, 1);

% solver colors for reporting and plotting (cf. get_style)
solverColors = lines(n_s);

% Shifted geometric means (cf. shifted_geometric_means in ClarabelBenchmarks)

% shift k of the shifted geometric mean
k_shift = 1;

% enforce the solve time cap; failed solves count at the cap
t_capped             = min(t_ps, commonOpts.maxSolveTime);
t_capped(~solved_ps) = commonOpts.maxSolveTime;

% shifted geometric mean g_s = exp( mean_p log(t_{p,s} + k) ) - k
g_s = exp(mean(log(t_capped + k_shift), 1)) - k_shift;

% normalized shifted geometric mean r_s: best solver scores 1
r_s = g_s / min(g_s);

% failure rate of each solver in percent
failrate_s = 100 * sum(~solved_ps, 1) / n_p;

% report the score table
fprintf('Shifted geometric mean (shift k = %g, time cap %g s):\n', k_shift, commonOpts.maxSolveTime)
for s = 1:n_s
    fprintf('  %-18s: g_s = %8.2f s, normalized score r_s = %5.2f, failure rate = %5.1f %%\n', ...
        solverNames{s}, g_s(s), r_s(s), failrate_s(s))
end

% Relative performance profile (cf. performance_profile in ClarabelBenchmarks)

% failed solves have no valid time: their ratio is Inf, i.e. they never
% count within the plotted range (this is Dolan & More's r_M convention)
t_valid             = t_ps;
t_valid(~solved_ps) = Inf;

% relative performance ratio u_{p,s} = t_{p,s} / min{ t_{p,s} : s in S }
% (u_{p,s} is exactly the performance ratio r_{p,s} of Dolan & More;
%  the baseline minimum runs over successful solves only)
u_ps = t_valid ./ min(t_valid, [], 2);

% performance levels: tau = 10^k, k = 0 : 0.01 : 3 (cf. ClarabelBenchmarks)
tau   = 10.^(0:0.01:3);
f_s_r = perf_profile_cdf(u_ps, tau);

% Absolute performance profile (cf. time_profile in ClarabelBenchmarks)

% axis range: powers of ten bracketing the successful solve times
min_time = 10^floor(log10(min(t_ps(solved_ps))));
max_time = 10^ceil(log10(max(t_ps(solved_ps))));

% exact step profile per solver at the sorted successful solve times;
% failed solves raise no step, so the curve saturates below one
for s = 1:n_s
    % sorted solve times of the successful runs of solver s
    t_s = sort(t_ps(solved_ps(:, s), s));
    % fraction of problems solved within each time (denominator: all n_p)
    y_s = (1:numel(t_s)).' / n_p;
    % extend the step function to the axis corners
    if isempty(t_s)
        t_prof_s{s} = [min_time, max_time];
        y_prof_s{s} = [0, 0];
    else
        t_prof_s{s} = [min_time; t_s; max_time];
        y_prof_s{s} = [0; y_s; y_s(end)];
    end
end

% Plot: two-panel performance profiles (cf. Fig. 1 in Clarabel)
figure('Name', 'Performance Profiles')
tiledlayout(1,2)

% (a) relative performance profile
ax_r = nexttile;
for s = 1:n_s
    semilogx(ax_r, tau, f_s_r(:, s), 'LineWidth', 2, 'Color', solverColors(s, :)), hold(ax_r, 'on')
end
grid(ax_r, 'on'), grid(ax_r, 'minor')
xlim(ax_r, [1, 100]), ylim(ax_r, [0, 1])
xticks(ax_r, [1, 10, 100])
yticks(ax_r, 0:0.2:1)
xlabel(ax_r, 'performance ratio $\tau$')
ylabel(ax_r, 'ratio of problems solved')
title(ax_r, 'Performance ratio')
legend(ax_r, solverNames, 'Location', 'southeast')

% (b) absolute performance profile
ax_a = nexttile;
for s = 1:n_s
    semilogx(ax_a, t_prof_s{s}, y_prof_s{s}, 'LineWidth', 2, 'Color', solverColors(s, :)), hold(ax_a, 'on')
end
grid(ax_a, 'on'), grid(ax_a, 'minor')
xlim(ax_a, [min_time, max_time]), ylim(ax_a, [0, 1])
xticks(ax_a, 10.^(log10(min_time):log10(max_time)))
yticks(ax_a, 0:0.2:1)
xlabel(ax_a, 'solve time $t$ [s]')
ylabel(ax_a, 'fraction of problems solved within $t$')
title(ax_a, 'Solution time profile')
legend(ax_a, solverNames, 'Location', 'southeast')

