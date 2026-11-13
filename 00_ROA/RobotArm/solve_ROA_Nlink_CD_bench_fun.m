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

% de-pack for better readability
x       = probStruct.x;
f       = probStruct.f;
g       = probStruct.g;
s1      = probStruct.s1;
V       = probStruct.V;
b0      = probStruct.b;

% enforce positivity
l = 1e-6*(x'*x);


%% setup solver
startBuildTime = tic;

% solver 1: s1-step
sos1 = struct('x',s1,...            % decision variables
             'p',V);            % parameter

% constraint
sos1.('g') = [
              s1;
              s1*(V-b0)-nabla(V,x)*f-l
             ];

% states + constraint are SOS cones
opts.Kx           = struct('lin', 1);
opts.Kc           = struct('sos', 2);
opts.newton_solver = [];
opts.error_on_fail = 0; % we do our own check in the main loop

if  strcmpi(commonOpts.solver,'mosek')

    opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
    opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

elseif strcmpi(commonOpts.solver,'clarabel')

    opts.sdpsol_options.clarabel = commonOpts.sdp_solver_opts;

elseif strcmpi(commonOpts.solver,'scs')

    opts.sdpsol_options.scs = commonOpts.sdp_solver_opts;

end

% build first solver
S1 = casos.sossol('S1',commonOpts.solver,sos1,opts);

% solver 2: V-step
cost = dot(g-V,g-V);            % squared distance of V to the safe set (no level-set shift)

sos2 = struct('x',V,...             % decision variables
    'f',cost,...                   % cost
    'p',s1);                       % parameter

% constraints
sos2.('g') = [
              V-l;
              s1*(V-b0)-nabla(V,x)*f-l
             ];

% states + constraint cones
opts      = struct;
opts.Kx   = struct('sos', 0, 'lin', 1);
opts.Kc   = struct('sos', 2);
opts.newton_solver = [];
opts.error_on_fail = 0; % we do our own check in the main loop

if  strcmpi(commonOpts.solver,'mosek')

    opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
    opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

elseif strcmpi(commonOpts.solver,'clarabel')

    opts.sdpsol_options.clarabel = commonOpts.sdp_solver_opts;

elseif strcmpi(commonOpts.solver,'scs')

    opts.sdpsol_options.scs = commonOpts.sdp_solver_opts;

end

% build third solver
S2 = casos.sossol('S',commonOpts.solver,sos2,opts);

% functions for KKT condition

% common constraints for Lagrangian
g_common  = [
             s1;
             V-l;
             s1*(V-b0)-nabla(V,x)*f-l
            ];

% same order as initial guess for decision variables
sos_x     = [V;s1];

% get lagrange multiplier (sparsity)
lam_gs    = casos.PS.sym('lam_gs', sparsity(g_common));

% lagrangian
lagrangian = casos.PS(cost + dot(lam_gs,g_common));

% derivative lagrangian
eval_gradLag  =  casos.Function('f',{sos_x,lam_gs}, { norm(jacobian( poly2basis(lagrangian), poly2basis(sos_x) ),inf) });

% inf. norm of constraints (for scaled termination criterion)
eval_gradLag_inf_norm  =  casos.Function('f',{sos_x,lam_gs}, { norm(poly2basis(casos.PS(dot(lam_gs,g_common))),inf) });

% cost function of the V-step (needed to evaluate the cost at the zeroth
% iterate/initial guess, cf. SequentialCommon: eval_cost of full(nlsos.f))
eval_cost = casos.Function('f',{sos_x}, { cost });

totalBuildTime = toc(startBuildTime);

%% Preparation

% initial guess
Vval = x0(1);

fprintf('------------------------------------------------------------------------------------------\n')
fprintf('Custom Coordinate Descent for Region-of-Attraction\n')
fprintf('------------------------------------------------------------------------------------------\n')
fprintf('-------------------------------------------------------\n')
fprintf('Settings \n ')
fprintf('Solver\t\t\t\t- %s\n ', commonOpts.solver);
fprintf('Maximum iterations\t\t- %d\n ', commonOpts.maxiter);
fprintf(['' char(949) '_opt\t\t\t\t- %d\n'],commonOpts.tolerance_opt);
fprintf('-------------------------------------------------------\n');
fprintf('\n');
fprintf('Solver build time \t\t:  %0.2d s\n', totalBuildTime);
fprintf('Start solving ... \n');
fprintf('\n');

% x_prev for termination criterion; take initial guess as first "previous"
% iteration
x_prev = x0;

% fallback final solution (same order as the sequential solver) in case the
% loop breaks before the first full iterate
sol_x_it = x0;

UNIFIED_RETURN_STATUS = 'SOLVER_RET_SUCCESS';

% number of initial-guess restarts needed until the problem became feasible
% (0: first guess feasible, 1: second starter needed, 2: third starter needed)
initial_guess_trials = 0;

% first-order KKT residual per iteration; entries beyond the last completed
% iteration stay inf (e.g. if a subsolver was infeasible)
nabla_L        = inf;
nabla_L_scaled = inf;

% no subsolver time at the zeroth iterate
solvetime_all1(1) = 0;
solvetime_all2(1) = 0;

% zeroth iterate (initial guess)
cost_CD(1) = full(eval_cost(x0));

% no multipliers available at the initial guess -> use zeros (cf. sequential)
lam_gs_0 = zeros(numel(lam_gs),1);

% no previous iterate at the zeroth iterate (x_prev = x0)
inf_delta_x_0 = norm(full(poly2basis(x_prev)),inf);

nabla_L(1)        = full(eval_gradLag(x0,lam_gs_0));
inf_g_0           = full(eval_gradLag_inf_norm(x0,lam_gs_0));
nabla_L_scaled(1) = nabla_L(1)*inf_delta_x_0/(max(1,abs(cost_CD(1)) + inf_g_0));

% display output header
fprintf('------------------------------------------------------------------------------------------\n');
fprintf('%-8s%-15s%-15s\n','iter','obj','||dLdx||');
fprintf('------------------------------------------------------------------------------------------\n');

% display output current iterate (zeroth iteration = initial guess)
fprintf('%-8d%-15g%-15g\n',0,cost_CD(1),nabla_L(1));

buildBisecTime = 0;
solveBisecTime  = 0;
solved = 1;

% counter until solution is marked almost optimal
counter_almost_opt_iter    = 0;

% flag to indicate that the maximum solve time was reached
time_limit_reached = false;

%% Main-loop (V-s-iteration)
% real k-th iteration runs at iter = k+1 (cf. sequential)

maxIter = commonOpts.maxiter;
maxTime = commonOpts.maxSolveTime;

if ~isfield(commonOpts,'almost_opt_count'), commonOpts.almost_opt_count = 10; end

startSolve_mainLoop = tic;
for iter = 2:maxIter+1


    %% s1-step
    s1start = tic;

    sol1 = S1('p',Vval);

    s1_ok = strcmp(S1.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS');
    if s1_ok
        s1val = sol1.x;
        UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
    else
        % if we just made a few iterations and become infeasible we
        % assume the initial guess is bad
        if iter <= 2

            s1val = x0(2);

            initial_guess_trials = 1;

            fprintf('Infeasible in solver 1 in iteration %d. We provide another initial guess. \n',iter-1)

        else
            fprintf('Infeasible in solver 1 in iteration %d; retaining previous iterate\n', iter-1);
            UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
            break
        end
    end

    solvetime_all1(iter) =  toc(s1start);

    % V-step
    S2start = tic;

    sol2 = S2('p',s1val);

    cost_CD(iter) = full(sol2.f);

    solvetime_all2(iter) = toc(S2start);

    if strcmp(S2.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
        % extract solution
        Vval = sol2.x(1);
    else
        if iter <= 2
            fprintf('Infeasible in solver 2 in iteration %d \n',iter-1)
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
        else
            fprintf('Infeasible in solver 2 in iteration %d; retaining previous iterate\n', iter-1);
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
        end
        break
    end

    % extract lagrange multiplier from each subproblem
    if s1_ok
        % in same order as constraints
        lam_gs_it = [sol1.lam_g(1);sol2.lam_g];
    else
        % in case S1 was infeasible, set multiplier to zero
        lam_gs_it = [0;sol2.lam_g];
    end

    % get complete solution vector, i.e, same order
    sol_x_it  = [sol2.x;sol1.x];

    % for scaled termination
    delta_x     = sol_x_it - x_prev;
    inf_delta_x = norm(full(poly2basis(delta_x)),inf);

    % evaluate first order KKT at current iterate
    nabla_L(iter) = full(eval_gradLag( sol_x_it,lam_gs_it));

    % inf. norm of lambda*constraint(s)
    inf_g   = full(eval_gradLag_inf_norm(sol_x_it,lam_gs_it));

    % display output header every 10th iteration
    if ~mod(iter-1,10) && iter-1 > 0
        fprintf('%-8s%-15s%-15s\n','iter','obj','||dLdx||');
        fprintf('------------------------------------------------------------------------------------------\n');
    end

    % display output current iterate
    fprintf('%-8d%-15g%-15g\n',iter-1,cost_CD(iter),nabla_L(iter));

    % just for plotting
    nabla_L_scaled(iter) = nabla_L(iter)*inf_delta_x/(max(1,abs(full(sol2.f)) + inf_g ));

    % scaled termination criterion
    if nabla_L(iter) <= commonOpts.tolerance_opt*(max(1,abs(full(sol2.f)) + inf_g ))/inf_delta_x
        break
    else
        % set current iterate to previous iterate
        x_prev =  sol_x_it;
    end

    % check if current iterate is solved to an acceptable level
    if nabla_L(iter) <= 10*commonOpts.tolerance_opt*(max(1,abs(full(sol2.f)) + inf_g ))/inf_delta_x

        % increase counter
        counter_almost_opt_iter = counter_almost_opt_iter +1;

        % we have 5-10 in a row, seems we solved it to an acceptable level/almost optimal
        if counter_almost_opt_iter >= commonOpts.almost_opt_count
            solved = 4;
            break
        end
    else
        % if current iterate is not solved to an acceptable level, reset counter
        counter_almost_opt_iter = 0;
    end

    % at the end of the loop, check solve time
    if toc(startSolve_mainLoop) >= maxTime && iter <= maxIter
        % abort: maximum solve time reached
        time_limit_reached = true;
        break
    end

end % end for-loop (inner)

totalSolveTime = toc(startSolve_mainLoop) - buildBisecTime;

% number of real CD iterations; the zeroth iterate (initial guess) is stored
% at index 1 and is not counted (cf. sequential solver: stats.iterations = iter-1)
numIter = iter - 1;

%% Check and prepare output

% ************************************************************************
% check if coordinate-descent was infeasible after a few iterations
% Note: in this case we did not reach the termination criterion, but have a
% feasible solution

% Meaning: the run only counts as successful if a subsolver became
% infeasible after a minimum number of iterations
% and the cost function made some progress over these iterations;
% otherwise the run is marked as RET_LIMITED
if ~strcmp(UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS') && iter <= maxIter

    if numIter >= commonOpts.miniter && abs(cost_CD(1) - cost_CD(end)) > commonOpts.tolerance_opt

        UNIFIED_RETURN_STATUS = 'SOLVER_RET_SUCCESS';
        solved = 2;
    else
        UNIFIED_RETURN_STATUS = 'RET_LIMITED';
    end
end


% check the case with max. iter i.e. it might be CD just needs many
% iterations to converge according to KKT conditions, yet, the solution
% might be already close

% we are feasible, but are not converged until max. iter or max. time
if (strcmp(UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS') && iter >= maxIter) || time_limit_reached
    UNIFIED_RETURN_STATUS = 'RET_LIMITED';
end

% summary output
fprintf('------------------------------------------------------------------------------------------ \n')
if strcmp(UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS') && solved == 1
    fprintf('Solution status: Optimal solution found\n')
elseif solved == 2
    fprintf('Solution status: Abort midway. Solved to feasible point, but did not reach KKT point\n')
elseif solved == 4
    fprintf('Solution status: Almost optimal.\n')
elseif time_limit_reached
    fprintf('Solution status: Maximum solve time reached\n')
elseif strcmp(UNIFIED_RETURN_STATUS,'RET_LIMITED') && iter > maxIter
    fprintf('Solution status: Maximum number of iterations reached\n')
else
    fprintf('Solution status: Problem infeasible.\n')
end

fprintf('Solution time: %.4f s\n',totalSolveTime)
fprintf('Build time: %.4f s\n',totalBuildTime + buildBisecTime)
fprintf('Total time: %.4f s\n',totalBuildTime + buildBisecTime + totalSolveTime)
fprintf('------------------------------------------------------------------------------------------ \n')

% pad diagnostic arrays to numIter+1 entries (index 1 = zeroth iterate) so
% that they line up with 0:stats.iterations in the convergence plots
cost_CD(end+1:numIter+1)        = cost_CD(end);
nabla_L(end+1:numIter+1)        = inf;
nabla_L_scaled(end+1:numIter+1) = inf;

stats.iterations            = numIter;
stats.UNIFIED_RETURN_STATUS = UNIFIED_RETURN_STATUS;
stats.cost                  = cost_CD;
stats.nabla_L               = nabla_L;
stats.nabla_L_scaled        = nabla_L_scaled;

stats.initial_guess_trials  = initial_guess_trials;

stats.totalSolveTime     = totalSolveTime;
stats.totalBuildTime     = totalBuildTime + buildBisecTime;
stats.time_per_iteration = stats.totalSolveTime / stats.iterations;

% write to struct
stats.solver1_stats      = S1.stats;
stats.solver2_stats      = S2.stats;

% individual solver times
stats.totalSolveTime_S1     =  sum(solvetime_all1);
stats.totalSolveTime_S2     =  sum(solvetime_all2);
stats.totalSolveTime_bisect = solveBisecTime;


% store final solution (same order as the sequential solver: [V;s1])
sol.x = sol_x_it;

solStruct.sol          = sol;
solStruct.stats        = stats;

% store information about sub-solver e.g. for problem size
solStruct.info.solver1 =  S1.info;
solStruct.info.solver2 =  S2.info;

% store problem size explicitly
solStruct.info.probSize_S1.A   = S1.info.sdpsol.conic.size_a;
solStruct.info.probSize_S1.gram = S1.info.gram.Kc;

solStruct.info.probSize_S2.A   = S2.info.sdpsol.conic.size_a;
solStruct.info.probSize_S2.gram = S2.info.gram.Kc;

end
