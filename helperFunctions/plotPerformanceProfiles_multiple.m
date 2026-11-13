function plotPerformanceProfiles_multiple(seqMethodStructsAll, cdSolverStructsAll, solverNamesAll, benchmarkNamesAll, maxSolveTimeAll, probNamesAll)

fprintf('-----------------------------------------------------------------------------------------------------\n')
fprintf('-------------------------------------- Benchmark Evalaution -----------------------------------------\n')
fprintf('-----------------------------------------------------------------------------------------------------\n')
fprintf('\n')

numBenchmarks = numel(seqMethodStructsAll);

n_sAll        = zeros(1, numBenchmarks);
t_psAll       = cell(1, numBenchmarks);
solved_psAll  = cell(1, numBenchmarks);

for mm = 1:numBenchmarks

% number of problems n_p: one stored run per cell entry
% (if no coordinate-descent runs are given, e.g. for CBF-CLF, take the
% number of sequential runs instead)
if isempty(cdSolverStructsAll{mm})
    n_p = numel(seqMethodStructsAll{mm}{1});
else
    n_p = numel(cdSolverStructsAll{mm});
end

% whether coordinate-descent runs are part of the comparison
hasCD = ~isempty(cdSolverStructsAll{mm});

% number of solvers n_s: one sequential solver per Hessian estimation method and coordinate descent
n_s = numel(solverNamesAll{mm});

% computing time t_{p,s}: one row per problem p in P, one column per solver s in S
% extract each run from the cell arrays of stored solver structs
t_ps = zeros(n_p, n_s);
for p = 1:n_p
    % sequential runs of problem p, one per Hessian estimation method
    for s = 1:numel(seqMethodStructsAll{mm})
        solStruct_seq = seqMethodStructsAll{mm}{s}{p};
        t_ps(p, s) = solStruct_seq.stats.total_solve_time;
    end
    % coordinate-descent run of problem p
    if hasCD
        solStruct_cd  = cdSolverStructsAll{mm}{p};
        % store the computing times t_{p,s} of all solvers
        t_ps(p, n_s) = solStruct_cd.stats.totalSolveTime;
    end
end

% logical matrix: solved_ps(p,s) = true iff solver s solves problem p
% a run only counts as solved for the return status 'SOLVER_RET_SUCCESS'
solved_ps = false(n_p, n_s);
for p = 1:n_p
    % sequential runs of problem p, one per Hessian estimation method
    for s = 1:numel(seqMethodStructsAll{mm})
        solStruct_seq = seqMethodStructsAll{mm}{s}{p};
        solved_ps(p, s) = strcmp(solStruct_seq.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS');
    end
    % coordinate-descent run of problem p
    if hasCD
        solStruct_cd  = cdSolverStructsAll{mm}{p};
        % solved_ps(p,n_s): coordinate descent
        solved_ps(p, n_s) = strcmp(solStruct_cd.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS');
    end
end

% drop problems that no solver could solve
keep_p = any(solved_ps, 2);

% additionally drop problems for which the solvers converged to
% different local solutions (cf. Eq. (38) in Waechter & Biegler,
% Math. Program. 106, 2004, p. 49): the final objective values of the
% unsuccessful runs are omitted and the problem is discarded if
% (f_*^max - f_*^min) / (1 + max(|f_*^min|, |f_*^max|)) > 1e-1
keep_cost = true(n_p, 1);
for p = 1:n_p
    % problems already dropped above cannot enter the comparison
    if ~keep_p(p)
        continue
    end
    % final objective values f(x_*^s) of the successful runs of problem p
    f_star = nan(n_s, 1);

    % sequential runs of problem p, one per Hessian estimation method
    for s = 1:numel(seqMethodStructsAll{mm})
        if solved_ps(p, s)
            solStruct_seq = seqMethodStructsAll{mm}{s}{p};
            f_star(s) = solStruct_seq.stats.single_iterations{end}.cost;
        end
    end

    % coordinate-descent run of problem p
    if hasCD && solved_ps(p, n_s)
        solStruct_cd = cdSolverStructsAll{mm}{p};
        f_star(n_s) = solStruct_cd.stats.cost(end);
    end

    % f_*^max and f_*^min over the final objective values of the
    % successful runs only (unsuccessful runs omitted)
    f_star_max = max(f_star, [], 'omitnan');
    f_star_min = min(f_star, [], 'omitnan');

    % discard the problem if the final objective values are not close
    if (f_star_max - f_star_min) / (1 + max(abs(f_star_min), abs(f_star_max))) > 1e-1
        keep_cost(p) = false;
    end
end

keep_p = keep_p & keep_cost;

% report which benchmarks are considered in the performance profile:
% sorted out are problems that are infeasible for all solvers as well as
% problems whose successful runs converged to different local solutions
fprintf('Benchmark selection for the performance profile (%s, %d of %d considered) \n', benchmarkNamesAll{mm}, nnz(keep_p), n_p)

for p = 1:n_p
    if ~any(solved_ps(p, :))
        fprintf('  sorted out: %s (infeasible for all solvers)\n', probNamesAll{mm}{p})
    elseif ~keep_p(p)
        fprintf('  sorted out: %s (diverging local costs, Eq. 38)\n', probNamesAll{mm}{p})
    else
        fprintf('  considered: %s\n', probNamesAll{mm}{p})
    end
end
fprintf('\n')
% remove those runs
t_ps(~keep_p, :)      = [];
solved_ps(~keep_p, :) = [];

n_sAll(mm)      = n_s;
t_psAll{mm}     = t_ps;
solved_psAll{mm} = solved_ps;

end

n_p = 0;
for mm = 1:numBenchmarks
    n_p = n_p + size(t_psAll{mm}, 1);
end
n_s = sum(n_sAll);

t_ps          = nan(n_p, n_s);
solved_ps     = false(n_p, n_s);
maxSolveTime_ps = nan(n_p, n_s);
legendNames   = cell(1, n_s);

rowOffset = 0;
colOffset = 0;
for mm = 1:numBenchmarks
    rows_mm = size(t_psAll{mm}, 1);
    t_ps(rowOffset+1:rowOffset+rows_mm, colOffset+1:colOffset+n_sAll(mm))                   = t_psAll{mm};
    solved_ps(rowOffset+1:rowOffset+rows_mm, colOffset+1:colOffset+n_sAll(mm))              = solved_psAll{mm};
    maxSolveTime_ps(rowOffset+1:rowOffset+rows_mm, colOffset+1:colOffset+n_sAll(mm))         = maxSolveTimeAll{mm};
    for s = 1:n_sAll(mm)
        legendNames{colOffset+s} = [benchmarkNamesAll{mm} ': ' solverNamesAll{mm}{s}];
    end
    rowOffset = rowOffset + rows_mm;
    colOffset = colOffset + n_sAll(mm);
end
applicable_ps = ~isnan(t_ps);

% number of problems n_p remaining in the test set
n_p = size(t_ps, 1);

% solver colors for reporting and plotting
solverColors = lines(n_s);

% shift k of the shifted geometric mean
k_shift = 1;

% enforce the solve time cap; failed solves count at the cap
t_capped                                  = min(t_ps, maxSolveTime_ps);
t_capped(~solved_ps & applicable_ps)      = maxSolveTime_ps(~solved_ps & applicable_ps);

% shifted geometric mean g_s = exp( mean_p log(t_{p,s} + k) ) - k
g_s = exp(mean(log(t_capped + k_shift), 1, 'omitnan')) - k_shift;

% normalized shifted geometric mean r_s: best solver scores 1
r_s = g_s / min(g_s);

% failure rate of each solver in percent
failrate_s = 100 * sum(~solved_ps & applicable_ps, 1) ./ sum(applicable_ps, 1);

% report the score table
fprintf('-----------------------------------------------------------------------------------------------------\n')
fprintf('\n')
fprintf('Shifted geometric mean (shift k = %g):\n', k_shift)
fprintf('\n')
for s = 1:n_s
    fprintf('  %-60s: g_s = %8.2f s, normalized score r_s = %5.2f, failure rate = %5.1f %%\n', ...
        legendNames{s}, g_s(s), r_s(s), failrate_s(s))
end
fprintf('-----------------------------------------------------------------------------------------------------\n')
fprintf('\n')
% Relative performance profile 

% failed solves have no valid time: their ratio is Inf, i.e. they never
% count within the plotted range (this is Dolan & More's r_M convention)
t_valid                                  = t_ps;
t_valid(~solved_ps & applicable_ps)      = maxSolveTime_ps(~solved_ps & applicable_ps);

% relative performance ratio u_{p,s} = t_{p,s} / min{ t_{p,s} : s in S }
% (u_{p,s} is exactly the performance ratio r_{p,s} of Dolan & More;
%  the baseline minimum runs over successful solves only)
u_ps = t_valid ./ min(t_valid, [], 2, 'omitnan');

% performance levels: tau = 10^k, k = 0 : 0.01 : 3 (cf. ClarabelBenchmarks.jl)
tau   = 10.^(0:0.01:3);
f_s_r = zeros(numel(tau), n_s);
for s = 1:n_s
    f_s_r(:, s) = perf_profile_cdf(u_ps(~isnan(u_ps(:, s)), s), tau);
end


% axis range: powers of ten bracketing the successful solve times
t_success = t_ps(solved_ps & applicable_ps);
min_time = 10^floor(log10(min(t_success)));
max_time = 10^ceil(log10(max(t_success)));

% exact step profile per solver at the sorted successful solve times;
% failed solves raise no step, so the curve saturates below one
for s = 1:n_s

    % sorted solve times of the successful runs of solver s
    t_s = sort(t_ps(solved_ps(:, s), s));

    % fraction of problems solved within each time (denominator: all n_p)
    y_s = (1:numel(t_s)).' / sum(applicable_ps(:, s));

    % extend the step function to the axis corners
    if isempty(t_s)
        t_prof_s{s} = [min_time, max_time];
        y_prof_s{s} = [0, 0];
    else
        t_prof_s{s} = [min_time; t_s; max_time];
        y_prof_s{s} = [0; y_s; y_s(end)];
    end
end

%% Plot performance profiles

figure('Name', 'Performance Profiles')
tiledlayout(1,2) %  relative performance profile
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
legend(ax_r, legendNames, 'Location', 'southeast')


ax_a = nexttile; % absolute performance profile
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
legend(ax_a, legendNames, 'Location', 'southeast')


fprintf('Detailed results in tabular form:\n')
fprintf('\n')

end
