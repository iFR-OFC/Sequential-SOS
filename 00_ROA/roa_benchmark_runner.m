%% ***********************************************************************
%
%  Region-of-attraction Estimation Benchmark Test Runner
% 
%  Run 1: Van-der-Pol Oscillator 
%  Run 2: Aircraft 4D
%  Run 3: F-18
%  Run 4: Network
%
% ***********************************************************************
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');

close all
clear
clc

addpath("F18/")
addpath("VDP/")
addpath("Aircraft_4D/")
addpath("12State/")


commonOpts.maxiter       = 350;
commonOpts.tolerance_opt = 1e-4;
% maximum allowable solve time; assigned to failed solves (Clarabel convention)
commonOpts.maxSolveTime  = 3600;


%% VDP

% generate problem structure (dynamics, polynomial dec. variables)
degV = {2,2:4,1:4,2:6};
degS = {2,2,2:4,2:4};

% initialize the cell arrays of stored runs (one entry per run)
seqSolverStructs = {};
cdSolverStructs  = {};

for kk = 1:numel(degV)


probStruct_VDP = generateVDP_benchmark(degV{kk},degS{kk});
fprintf('-------------------------------------------------- \n')
fprintf('Run 1 - Region-of-Attraction for %s. \n',probStruct_VDP.Name)
fprintf('-------------------------------------------------- \n')

% generate an initial guess
x0 = generate_initGuess_VDP(probStruct_VDP);

% solve VDP example with two states using sequential SOS
solStruct_seq_VDP     = solve_ROA_seq_bench_fun(probStruct_VDP,commonOpts,x0);

seqSolverStructs{end+1} = solStruct_seq_VDP;

% extract size and timing statistic from sequential
evalSequentialStruct_VDP = evalSequential(solStruct_seq_VDP);

% solve VDP example with two states using custom coordinate-descent
solStruct_CD_VDP = solve_ROA_CD_bench_fun(probStruct_VDP,commonOpts,x0);

cdSolverStructs{end+1} = solStruct_CD_VDP;


% Evaluate

% Plot convergence diagnostics
figure('Name','VDP Benchmark')
tiledlayout(1,2)
nexttile
semilogy(abs(evalSequentialStruct_VDP.cost(2:end)-evalSequentialStruct_VDP.cost(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,
semilogy(abs(solStruct_CD_VDP.stats.cost(2:end)-solStruct_CD_VDP.stats.cost(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,
xlabel('Iterations')
ylabel('$ |f(\xi_{k+1}) - f(\xi_{k}) |$')
legend('Sequential','Coordinate-descent')

% constraint violation of sequential SOS
% ------------------------------------------------------------------------
% Note: If we are in the SOS cone, constraint violation is negative if
% signed distance is used. In that case, constraint violation is set to
% zero. For plotting on logarithmic scale zero values are set to eps!
% ------------------------------------------------------------------------
nexttile, semilogy(evalSequentialStruct_VDP.maxConVio,'LineWidth',1.2), grid on, title('Maximum constraint violation') , hold on 
plot([0 solStruct_seq_VDP.stats.iterations],[solStruct_seq_VDP.stats.options.tolerance_con solStruct_seq_VDP.stats.options.tolerance_con],'k--','LineWidth',1.2)
plot([0 solStruct_seq_VDP.stats.iterations],[eps eps],'k-.','LineWidth',1.2)
xlabel('Iterations')
ylabel('$\min(\mathrm{eps},\theta(\xi_{k}))$')
legend('$ \theta(\xi_{k})$','$\epsilon_\mathrm{con}$','eps')
end

%% Aircraft 4D


% generate problem structure (dynamics, polynomial dec. variables)
degV = {2,2:4,1:4};
degS = {2:4,2:4,2:4};

for kk = 1:numel(degV)

% generate problem structure (dynamics, polynomial dec. variables)
probStruct_Aircraft = generateAircraft_benchmark(degV{kk},degS{kk});
fprintf('-------------------------------------------------- \n')
fprintf('Run 2 - Region-of-Attraction for %s. \n',probStruct_Aircraft.Name)
fprintf('-------------------------------------------------- \n')

% generate an initial guess
x0 = generate_initGuess_Aircraft(probStruct_Aircraft);

% solve Aircraft example with two states using sequential SOS
solStruct_seq_Aircraft     = solve_ROA_seq_bench_fun(probStruct_Aircraft,commonOpts,x0);

seqSolverStructs{end+1} = solStruct_seq_Aircraft;

% extract size and timing statistic from sequential
evalSequentialStruct_Aircraft = evalSequential(solStruct_seq_Aircraft);

% solve Aircraft example with two states using custom coordinate-descent
solStruct_CD_Aircraft = solve_ROA_CD_bench_fun(probStruct_Aircraft,commonOpts,x0);

cdSolverStructs{end+1} = solStruct_CD_Aircraft;

% Evaluate

% Plot convergence diagnostics
figure('Name','Aircraft Benchmark')
tiledlayout(1,2)
nexttile
semilogy(abs(evalSequentialStruct_Aircraft.cost(2:end)-evalSequentialStruct_Aircraft.cost(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,
semilogy(abs(solStruct_CD_Aircraft.stats.cost(2:end)-solStruct_CD_Aircraft.stats.cost(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,
xlabel('Iterations')
ylabel('$ |f(\xi_{k+1}) - f(\xi_{k}) |$')
legend('Sequential','Coordinate-descent')

% constraint violation of sequential SOS
% ------------------------------------------------------------------------
% Note: If we are in the SOS cone, constraint violation is negative if
% signed distance is used. In that case, constraint violation is set to
% zero. For plotting on logarithmic scale zero values are set to eps!
% ------------------------------------------------------------------------
nexttile, semilogy(evalSequentialStruct_Aircraft.maxConVio,'LineWidth',1.2), grid on, title('Maximum constraint violation') , hold on 
plot([0 solStruct_seq_Aircraft.stats.iterations],[solStruct_seq_Aircraft.stats.options.tolerance_con solStruct_seq_Aircraft.stats.options.tolerance_con],'k--','LineWidth',1.2)
plot([0 solStruct_seq_Aircraft.stats.iterations],[eps eps],'k-.','LineWidth',1.2)
xlabel('Iterations')
ylabel('$\min(\mathrm{eps},\theta(\xi_{k}))$')
legend('$ \theta(\xi_{k})$','$\epsilon_\mathrm{con}$','eps')

end
%% F/A-18

degV = {2,2:4,2:4};
degS = {2,2,2:4};

for kk = 1:numel(degV)
% generate problem structure (dynamics, polynomial dec. variables)
probStruct_F18 = generateF18_benchmark(degV{kk},degS{kk});

fprintf('-------------------------------------------------- \n')
fprintf('Run 3 - Region-of-Attraction for %s. \n',probStruct_F18.Name)
fprintf('-------------------------------------------------- \n')

% generate an initial guess
x0 = generate_initGuess_F18(probStruct_F18);

% solve F-18 example with seven states using sequential SOS
solStruct_seq_7D     = solve_ROA_seq_bench_fun(probStruct_F18,commonOpts,x0);

seqSolverStructs{end+1} = solStruct_seq_7D;

% extract size and timing statistic from sequential
evalSequentialStruct_7D = evalSequential(solStruct_seq_7D);

% solve F-18 example with seven states using custom coordinate-descent
solStruct_CD_7D = solve_ROA_CD_bench_fun(probStruct_F18,commonOpts,x0);

cdSolverStructs{end+1} = solStruct_CD_7D;

% Plot convergence diagnostics
figure('Name','12State Benchmark')
tiledlayout(1,2)
nexttile
semilogy(abs((evalSequentialStruct_7D.cost(2:end)-evalSequentialStruct_7D.cost(1:end-1))),'LineWidth',1.2), grid on, title('Cost'), hold on,
semilogy(abs(solStruct_CD_7D.stats.cost(2:end)-solStruct_CD_7D.stats.cost(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,
xlabel('Iterations')
ylabel('$ |f(\xi_{k+1}) - f(\xi_{k}) |$')
legend('Sequential','Coordinate-descent')

% constraint violation of sequential SOS
% ------------------------------------------------------------------------
% Note: If we are in the SOS cone, constraint violation is negative if
% signed distance is used. In that case, constraint violation is set to
% zero. For plotting on logarithmic scale zero values are set to eps!
% ------------------------------------------------------------------------
nexttile, semilogy(0:1:solStruct_seq_7D.stats.iterations, evalSequentialStruct_7D.maxConVio,'LineWidth',1.2), grid on, title('Maximum constraint violation') , hold on 
plot([0 solStruct_seq_7D.stats.iterations],[solStruct_seq_7D.stats.options.tolerance_con solStruct_seq_7D.stats.options.tolerance_con],'k--','LineWidth',1.2)
plot([0 solStruct_seq_7D.stats.iterations],[eps eps],'k-.','LineWidth',1.2)
xlabel('Iterations')
ylabel('$\min(\mathrm{eps},\theta(\xi_{k}))$')
legend('$ \theta(\xi_{k})$','$\epsilon_\mathrm{con}$','eps')

end

%% 12-State Dynamical System

% generate problem structure (dynamics, polynomial dec. variables)
probStruct_12State= generate12State_benchmark();
fprintf('-------------------------------------------------- \n')
fprintf('Run 4 - Region-of-Attraction for %s. \n',probStruct_12State.Name)
fprintf('-------------------------------------------------- \n')

% generate an initial guess
x0 = generate_initGuess_12State(probStruct_12State);

% solve 12 state dynamical system example with 12 states using sequential SOS
solStruct_seq_12D     = solve_ROA_seq_bench_fun(probStruct_12State,commonOpts,x0);

seqSolverStructs{end+1} = solStruct_seq_12D;

% extract size and timing statistic from sequential
evalSequentialStruct_12D = evalSequential(solStruct_seq_12D);

% solve 12 state dynamical system example with 12 states using custom coordinate-descent
solStruct_CD_12D = solve_ROA_CD_bench_fun(probStruct_12State,commonOpts,x0);

cdSolverStructs{end+1} = solStruct_CD_12D;


% Plot convergence diagnostics
figure('Name','12StateBenchmark')
tiledlayout(1,2)
nexttile
semilogy(abs((evalSequentialStruct_12D.cost(2:end)-evalSequentialStruct_12D.cost(1:end-1))),'LineWidth',1.2), grid on, title('Cost'), hold on,
semilogy(abs(solStruct_CD_12D.stats.cost(2:end)-solStruct_CD_12D.stats.cost(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,
xlabel('Iterations')
ylabel('$ |f(\xi_{k+1}) - f(\xi_{k}) |$')
legend('Sequential','Coordinate-descent')

% constraint violation of sequential SOS
% ------------------------------------------------------------------------
% Note: If we are in the SOS cone, constraint violation is negative if
% signed distance is used. In that case, constraint violation is set to
% zero. For plotting on logarithmic scale zero values are set to eps!
% ------------------------------------------------------------------------
nexttile, semilogy(0:1:solStruct_seq_12D.stats.iterations, evalSequentialStruct_12D.maxConVio,'LineWidth',1.2), grid on, title('Maximum constraint violation') , hold on 
plot([0 solStruct_seq_12D.stats.iterations],[solStruct_seq_12D.stats.options.tolerance_con solStruct_seq_12D.stats.options.tolerance_con],'k--','LineWidth',1.2)
plot([0 solStruct_seq_12D.stats.iterations],[eps eps],'k-.','LineWidth',1.2)
xlabel('Iterations')
ylabel('$\min(\mathrm{eps},\theta(\xi_{k}))$')
legend('$ \theta(\xi_{k})$','$\epsilon_\mathrm{con}$','eps')


%% Performance Profile (Dolan & More, Math. Program. 91, 2002)

% number of problems n_p: one stored run per cell entry
n_p = numel(seqSolverStructs);

% number of solvers n_s: sequential and coordinate descent
n_s = 2;

% computing time t_{p,s}: one row per problem p in P, one column per solver s in S
% extract each run from the cell arrays of stored solver structs
t_ps = zeros(n_p, n_s);
for p = 1:n_p
    % sequential run of problem p
    solStruct_seq = seqSolverStructs{p};
    % coordinate-descent run of problem p
    solStruct_cd  = cdSolverStructs{p};
    % store the computing times t_{p,s} of both solvers
    t_ps(p, :) = [solStruct_seq.stats.totalSolveTime, solStruct_cd.stats.totalSolveTime];
end

% logical matrix: solved_ps(p,s) = true iff solver s solves problem p
% a run only counts as solved for the return status 'SOLVER_RET_SUCCESS'
solved_ps = false(n_p, n_s);
for p = 1:n_p
    % sequential run of problem p
    solStruct_seq = seqSolverStructs{p};
    % coordinate-descent run of problem p
    solStruct_cd  = cdSolverStructs{p};
    % solved_ps(p,1): sequential, solved_ps(p,2): coordinate descent
    solved_ps(p, 1) = strcmp(solStruct_seq.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS');
    solved_ps(p, 2) = strcmp(solStruct_cd.stats.UNIFIED_RETURN_STATUS,  'SOLVER_RET_SUCCESS');
end

% drop problems that no solver could solve (cf. drop_unsolvable! in
% ClarabelBenchmarks): all-solvers-failed runs must not distort the
% profile denominators nor the shifted geometric means
keep_p = any(solved_ps, 2);
t_ps(~keep_p, :)      = [];
solved_ps(~keep_p, :) = [];

% number of problems n_p remaining in the test set
n_p = size(t_ps, 1);

% solver names and colors for reporting and plotting (cf. get_style)
solverNames  = {'Sequential', 'Coordinate-descent'};
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
semilogx(ax_r, tau, f_s_r(:, 1), 'LineWidth', 2, 'Color', solverColors(1, :)), hold(ax_r, 'on')
semilogx(ax_r, tau, f_s_r(:, 2), 'LineWidth', 2, 'Color', solverColors(2, :))
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

%% ------------------------------------------------------------------------
%  local functions
% ------------------------------------------------------------------------

% relative performance profile: fraction of problems per solver column
% whose ratio lies below each grid level,
% f_sx(k, s) = #{ p in P : data(p, s) < x_grid(k) } / n_p
% (strict comparison, as in ClarabelBenchmarks/src/performance_profile.jl;
%  Inf entries of failed solves never count)
function f_sx = perf_profile_cdf(data, x_grid)
    % logical 3-D array with entry (p, s, k) = (data(p, s) < x_grid(k)),
    % vectorized over all solvers via implicit expansion
    f_sx = squeeze(sum(data < reshape(x_grid, 1, 1, []), 1) / size(data, 1)).';
end