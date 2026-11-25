function cost_spread_calculation(seqMethodStructs, cdSolverStructs, solverNames, maxSolveTime, probNames)

fprintf('-----------------------------------------------------------------------------------------------------\n')
fprintf('-------------------------------------- Benchmark Evalaution -----------------------------------------\n')
fprintf('-----------------------------------------------------------------------------------------------------\n')
fprintf('\n')

% number of problems n_p: one stored run per cell entry
if isempty(cdSolverStructs)
    n_p = numel(seqMethodStructs{1});
else
    n_p = numel(cdSolverStructs);
end

% whether coordinate-descent runs are part of the comparison
hasCD = ~isempty(cdSolverStructs);

% optional: problem names for the evaluation output (one per run);
% if not provided, the runs are labeled by their index
if nargin < 5
    probNames = arrayfun(@(p) sprintf('Problem %d', p), (1:n_p).', 'UniformOutput', false);
end

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
    if hasCD
        solStruct_cd  = cdSolverStructs{p};
        % store the computing times t_{p,s} of all solvers
        t_ps(p, n_s) = solStruct_cd.stats.totalSolveTime;
    end
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
    if hasCD
        solStruct_cd  = cdSolverStructs{p};
        % solved_ps(p,n_s): coordinate descent
        solved_ps(p, n_s) = strcmp(solStruct_cd.stats.UNIFIED_RETURN_STATUS,  'SOLVER_RET_SUCCESS');
    end
end

% drop problems that no solver could solve
keep_p = any(solved_ps, 2);

% additionally drop problems for which the solvers converged to
% different local solutions 
% 
% see in Waechter & Biegler, Math. Program. 106, 2004, p. 49)
% if (f_*^max - f_*^min) / (1 + max(|f_*^min|, |f_*^max|)) > 1e-1, then
% discard

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
    if hasCD && solved_ps(p, n_s)
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


for k = 1:size(solved_ps,1)
    
    if solved_ps(k,3)
        cost_CD(k) = cdSolverStructs{k}.stats.cost(end);
    else
        cost_CD(k) = -1;
    end

     if solved_ps(k,1)
      cost_seq_hess_1(k) = seqMethodStructs{1}{k}.stats.single_iterations{end}.cost;
     else
         cost_seq_hess_1(k) = -1;
     end


     if solved_ps(k,2)
         cost_seq_hess_2(k) = seqMethodStructs{2}{k}.stats.single_iterations{end}.cost;
        else
         cost_seq_hess_2(k) = -1;
     end
        
    end




% minimum solve times
cost_CD(cost_CD == -1) =[];
cost_seq_hess_1(cost_seq_hess_1 == -1) =[];
cost_seq_hess_2(cost_seq_hess_2 == -1) =[];

min(cost_CD)
min(cost_seq_hess_1)
min(cost_seq_hess_2)


max(cost_CD)
max(cost_seq_hess_1)
max(cost_seq_hess_2)



% "global" minimum
min([min(cost_CD), min(cost_seq_hess_1),min(cost_seq_hess_2)])


% table: Method, Min cost, Max cost

% printf: "global minimum, total number of successfully solved problems

end
