%--------------------------------------------------------------------------
%
%   Short Description:  Calculate an inner-estimate of the
%                       region-of-attraction for the N-link robot arm
%                       using coordinate descent.
%
%                       Note: In contrast to the other ROA examples (see
%                       solve_ROA_CD_bench_fun.m), the squared distance of
%                       V itself to the safe set is minimized; the fixed
%                       level-set constant b enters the SOS constraints
%                       only. The structure, inputs, and outputs are
%                       identical to solve_ROA_CD_bench_fun.m.
%
%   Reference:  T. Cunis and J. Olucak,
%   "CaΣoS: A nonlinear sum-of-squares optimization suite⋆,"
%    2025 American Control Conference (ACC), Denver, CO,
%    USA, 2025, pp. 1659-1666, doi: 10.23919/ACC63710.2025.11107794.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%--------------------------------------------------------------------------

function solStruct = solve_ROA_Nlink_CD_bench_fun(probStruct,commonOpts,x0)

x = probStruct.x;
f = probStruct.f;
g = probStruct.g;
s1 = probStruct.s1;
V = probStruct.V;
b0 = probStruct.b;


% enforce positivity
l = 1e-6*(x'*x);

% start to measure build time of all parameterized solver
buildTimes_start = tic;


%% setup solver

% solver 1: s1-step
sos1 = struct('x',s1,'p',V);     % parameter

% constraint
sos1.('g') = [s1; s1*(V-b0)-nabla(V,x)*f-l];

% states + constraint are SOS cones
opts.Kx = struct('lin', 1);
opts.Kc = struct('sos', 2);
opts.newton_solver =[];
opts.error_on_fail = 0; % we do our own check in the main loop
opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

% build first solver
S1 = casos.sossol('S1',commonOpts.solver ,sos1,opts);

% solver 2: V-step
sos2 = struct('x',V, ...        % dec.var
    'f',dot(g-V,g-V),...        % squared distance of V to the safe set (no level-set shift)
    'p',s1);                % parameter

% constraints
sos2.('g') = [V-l;
    s1*(V-b0)-nabla(V,x)*f-l];

opts     = struct;
opts.Kx = struct('sos', 0, 'lin', 1);
opts.Kc = struct('sos', 2);
opts.newton_solver =[];
opts.error_on_fail = 0; % we do our own check in the main loop
opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

% build third solver
S2 = casos.sossol('S',commonOpts.solver ,sos2,opts);

totalBuildTime = toc(buildTimes_start);

%% Initial guess Lyapunov function

fval_old = [];


maxIter = commonOpts.maxiter;


Vval = x0(1);

fprintf('-------------------------------------------------- \n')
fprintf('Custom Coordinate Descent for Region-of-Attraction \n')
fprintf('-------------------------------------------------- \n')



UNIFIED_RETURN_STATUS = 'SOLVER_RET_SUCCESS';

% number of initial-guess restarts needed until the problem became feasible
% (0: first guess feasible, 1: second starter needed, 2: third starter needed)
initial_guess_trials = 0;

%% V-s-iteration
for iter = 1:maxIter

    %% s1-step
    s1start = tic;

        sol1 = S1('p',Vval);

        if strcmp(S1.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
            s1val  = sol1.x;
            % old_UNIFIED_RETURN_STATUS1 = S1.stats.UNIFIED_RETURN_STATUS;
        else
            % if we just made a few iterations and become infeasible we
            % assume
            if iter <= 1

                s1val = x0(2);

                initial_guess_trials = 1;

                fprintf('Infeasible in solver 1 in iteration %d. We provide another initial guess. \n',iter)

            else
                fprintf('Infeasible in solver 1 in iteration %d; retaining previous iterate\n', iter);
                UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
                break
            end
        end


    solvetime_all1(iter) =  toc(s1start);


    %% V-step
    S2start = tic;

    sol2 = S2('p',s1val);

    cost_CD(iter) = full(sol2.f);

    solvetime_all2(iter) = toc(S2start);

    if strcmp(S2.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
        % extract solution
        Vval = sol2.x(1);
        old_UNIFIED_RETURN_STATUS2 = S2.stats.UNIFIED_RETURN_STATUS;
    else
        if iter <= 1
            fprintf('Infeasible in solver 2 in iteration %d \n',iter)
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;

        else
            fprintf('Infeasible in solver 2 in iteration %d; retaining previous iterate\n', iter);
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
        end
        break
    end


    % show progress
    fprintf('Iteration %d: f = %g\n',iter,full(sol2.f));


    if ~isempty(fval_old)
        if abs(full(sol2.f - fval_old)) <= commonOpts.tolerance_opt
            break
        else
            fval_old = sol2.f;
        end
    else
        fval_old = sol2.f;
    end



end % end for-loop (inner)

totalSolveTime = sum(solvetime_all1) + sum(solvetime_all2);


if iter == maxIter
    UNIFIED_RETURN_STATUS = 'RET_LIMITED';
    stats.totalSolveTime = totalSolveTime;
end

% ************************************************************************
% check if coordinate-descent was infeasible after a few iterations
% Note: in this case we did not reach the termination criterion, but have a
% feasible solution

% Meaning: the run only counts as successful if a subsolver became
% infeasible after a minimum number of iterations (cf. commonOpts.miniter)
% and the cost function made some progress over these iterations;
% otherwise the run is marked as RET_LIMITED
if ~strcmp(UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS') && iter < maxIter
    if iter >= commonOpts.miniter && abs(cost_CD(1) - cost_CD(end)) > commonOpts.tolerance_opt
        UNIFIED_RETURN_STATUS = 'SOLVER_RET_SUCCESS';
    else
        UNIFIED_RETURN_STATUS = 'RET_LIMITED';
    end
end
% ************************************************************************

% summary output (cf. sequential solver)
fprintf('------------------------------------------------------------------------------------------ \n')
if strcmp(UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS')
    fprintf('Solution status: Optimal solution found\n')
elseif iter == maxIter
    fprintf('Solution status: Maximum number of iterations reached\n')
else
    fprintf('Solution status: Problem infeasible.\n')
end
fprintf('Solution time: %.4f s\n',totalSolveTime)
fprintf('Build time: %.4f s\n',totalBuildTime)
fprintf('Total time: %.4f s\n',totalBuildTime+totalSolveTime)
fprintf('------------------------------------------------------------------------------------------ \n')



stats.iterations     = iter;
stats.UNIFIED_RETURN_STATUS        = UNIFIED_RETURN_STATUS;
stats.cost           = cost_CD;
stats.initial_guess_trials = initial_guess_trials;

if strcmp(UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS') && iter < maxIter
    stats.totalSolveTime = totalSolveTime;
else
    if iter < maxIter
        stats.totalSolveTime = inf;
    end
end

stats.totalBuildTime = totalBuildTime;
stats.solver1_stats  = S1.stats;
stats.solver2_stats  = S2.stats;

% store final solution
sol.x.V = Vval;
sol.x.s1 = s1val;

solStruct.sol          = sol;
solStruct.stats        = stats;

% store information about subsolver e.g. for problem size
solStruct.inf0.solver1 =  S1.info;
solStruct.inf0.solver2 =  S2.info;


end
