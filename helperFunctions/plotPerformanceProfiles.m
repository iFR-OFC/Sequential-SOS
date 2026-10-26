function plotPerformanceProfiles(seqMethodStructs, cdSolverStructs, solverNames, maxSolveTime)

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
    for s = 1:numel(seqMethodStructs)
        if solved_ps(p, s)
            solStruct_seq = seqMethodStructs{s}{p};
            f_star(s) = solStruct_seq.stats.single_iterations{end}.cost;
        end
    end
    % coordinate-descent run of problem p
    if solved_ps(p, n_s)
        solStruct_cd = cdSolverStructs{p};
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
if any(keep_p & ~keep_cost)
    fprintf('Cost exclusion (Eq. 38): discarded %d problem(s) with diverging final objective values\n', ...
        sum(keep_p & ~keep_cost))
end
keep_p = keep_p & keep_cost;

t_ps(~keep_p, :)      = [];
solved_ps(~keep_p, :) = [];

% number of problems n_p remaining in the test set
n_p = size(t_ps, 1);

% solver colors for reporting and plotting (cf. get_style)
solverColors = lines(n_s);

% shift k of the shifted geometric mean
k_shift = 1;

% enforce the solve time cap; failed solves count at the cap
t_capped             = min(t_ps, maxSolveTime);
t_capped(~solved_ps) = maxSolveTime;

% shifted geometric mean g_s = exp( mean_p log(t_{p,s} + k) ) - k
g_s = exp(mean(log(t_capped + k_shift), 1)) - k_shift;

% normalized shifted geometric mean r_s: best solver scores 1
r_s = g_s / min(g_s);

% failure rate of each solver in percent
failrate_s = 100 * sum(~solved_ps, 1) / n_p;

% report the score table
fprintf('Shifted geometric mean (shift k = %g, time cap %g s):\n', k_shift, maxSolveTime)
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

end
