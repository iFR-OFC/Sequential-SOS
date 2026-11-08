function solStruct = solve_CbfClf_CD_bench_fun(probStruct,commonOpts,x0)


x  = probStruct.x;
u  = probStruct.u;
f  = probStruct.f;
g  = probStruct.g;
V  = probStruct.V;
h  = probStruct.h;
kappa = probStruct.kappa;
s1 = probStruct.s1;
s2 = probStruct.s2;
s3 = probStruct.s3;
s4 = probStruct.s4;
s5 = probStruct.s5;
gammaV = probStruct.gammaV;
gammaW = probStruct.gammaW;
b  = probStruct.b;

umin = probStruct.umin;
umax = probStruct.umax;

l = 1e-6*(x'*x);

% b = 1;
%% setup solver
startBuildTime = tic;

sos = struct('x',[s1;s2;s3;s4;s5;kappa],... % decision variables
    'p',[V;h;b]);    % parameter

% SOS constraints
sos.('g') = [
    s1;
    s2;
    s3;
    s4;
    s5;
    s1*(h-b) - g;                                         % state constraints
    s2*(h-b) - nabla(h,x)*subs(f,u,kappa) - gammaW*(h-b); % CBF dissipation
    s3*(h-b) + kappa - umin;                              % control constraints
    s4*(h-b) + umax - kappa;
    s5*(h-b) - nabla(V,x)*subs(f,u,kappa) - gammaV*V     % CLF dissipation
    ];

% states + constraint cones
opts.Kx           = struct('lin', length(sos.x));
opts.Kc           = struct('sos', length(sos.g));
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

% solver for multipliers, gammas and control law
S1 = casos.sossol('S1',commonOpts.solver,sos,opts);


% setup second solver
cost = dot(g-(h-b),g-(h-b));

sos = struct('x',[V;h;b],...                    % decision variables
    'f',cost,...                                % cost
    'p',[s1;s2;s3;s4;s5;kappa]);  % parameter

% SOS constraints
sos.('g') = [
    V - l;                                                % CLF positivity
    b - 1e-6;                                             % keep the level set positive
    s1*(h-b) - g;                                         % state constraints
    s2*(h-b) - nabla(h,x)*subs(f,u,kappa) - gammaW*(h-b); % CBF dissipation
    s3*(h-b) + kappa - umin;                              % control constraints
    s4*(h-b) + umax - kappa;
    s5*(h-b) - nabla(V,x)*subs(f,u,kappa) - gammaV*V     % CLF dissipation
    ];

% states + constraint cones
opts      = struct;
opts.Kx   = struct('lin', length(sos.x));
opts.Kc   = struct('sos', length(sos.g));
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

% solver for CLF, CBF and level set
S2 = casos.sossol('S2',commonOpts.solver,sos,opts);

% common constraints for Lagrangian (same order as sequential solver)
g_common  = [s1;
             s2;
             s3;
             s4;
             s5;
             V-l;
             b - 1e-6;
             s1*(h-b) - g;
             s2*(h-b) - nabla(h,x)*subs(f,u,kappa) - gammaW*(h-b);
             s3*(h-b) + kappa - umin;
             s4*(h-b) + umax - kappa;
             s5*(h-b) - nabla(V,x)*subs(f,u,kappa) - gammaV*V];

% same order as initial guess for decision variables
sos_x     = [V;h;s1;s2;s3;s4;s5;kappa;b];

% get lagrange multiplier
lam_gs    = casos.PS.sym('lam_gs', sparsity(g_common));

% lagrangian
lagrangian = casos.PS(cost + dot(lam_gs,g_common));

% derivative lagrangian
eval_gradLag  =  casos.Function('f',{sos_x,lam_gs}, { norm(jacobian( poly2basis(lagrangian), poly2basis(sos_x) ),inf) });

% inf. norm of constraints (for scaled termination criterion)
eval_gradLag_inf_norm  =  casos.Function('f',{sos_x,lam_gs}, { norm(poly2basis(casos.PS(dot(lam_gs,g_common))),inf) });

totalBuildTime = toc(startBuildTime);


%% Initial guess

fval_old = [];

maxIter = commonOpts.maxiter;

% initial guess for Prob 2: [V;h;b] taken from x0
% (block layout of x0: V,h,s1,s2,s3,s4,s5,kappa,b,gammaW,gammaV)
sol2_x = [x0(1);
          x0(2);
          x0(end)];

% initial guess for Prob 1: [s1;s2;s3;s4;s5;gammaW;gammaV;kappa] from x0
sol1_x0 = x0(3:end-1);

fprintf('-------------------------------------------------- \n')
fprintf('Custom Coordinate Descent for CBF-CLF \n')
fprintf('-------------------------------------------------- \n')

% x_prev for termination criterion; take initial guess as first "previous"
% iteration
x_prev = x0;

% fallback final solution (same order as the sequential solver) in case the
% loop breaks before the first full iterate
sol_x_it = x0;

UNIFIED_RETURN_STATUS = 'SOLVER_RET_SUCCESS';

% number of initial-guess restarts needed until the problem became feasible
initial_guess_trials = 1;

% first-order KKT residual per iteration; stays inf if the loop breaks
% before the first evaluation (e.g. a subsolver was infeasible)
nabla_L        = inf;
nabla_L_scaled = inf;

%% V-h-b-s-iteration
for iter = 1:maxIter


    %% s-step (Prob 1)
    s1start = tic;

    sol1 = S1('p',sol2_x);

    s1_ok = strcmp(S1.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS');
    if s1_ok
        sol1_x = sol1.x;
        UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
    else
        % if we just made a few iterations and become infeasible we
        % assume the initial guess is bad
        if iter <= 1

            sol1_x = sol1_x0;

            initial_guess_trials = 2;

            fprintf('Infeasible in solver 1 in iteration %d. We provide another initial guess. \n',iter)

        else
            fprintf('Infeasible in solver 1 in iteration %d; retaining previous iterate\n', iter);
            UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
            break
        end
    end

    solvetime_all1(iter) =  toc(s1start);


    %% V-h-b-step (Prob 2)
    S2start = tic;

    sol2 = S2('p',sol1_x);

    cost_CD(iter) = full(sol2.f);

    solvetime_all2(iter) = toc(S2start);

    s2_ok = strcmp(S2.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS');
    if s2_ok
        % extract solution
        sol2_x = sol2.x;
    else
        if iter <= 1
            fprintf('Infeasible in solver 2 in iteration %d. We try another initial guess! \n',iter)
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
            sol2_x = [x0(1);
                      x0(2)
                      x0(end)];
            initial_guess_trials = 3;
        else
            fprintf('Infeasible in solver 2 in iteration %d; retaining previous iterate\n', iter);
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
            break
        end
    end

    % show progress
    fprintf('Iteration %d: f = %g\n',iter,full(sol2.f));

    if s2_ok
    % lam of s1..s5 (s3,s4 are nu-vectors) from S1, rest from S2
    if s1_ok
        % in same order as constraints
        lam_gs_it = [sol1.lam_g(1:length([s1;s2;s3;s4;s5]));sol2.lam_g];
    else
        % in case S1 was infeasible, set multiplier to zero
        lam_gs_it = [zeros(length([s1;s2;s3;s4;s5]),1);sol2.lam_g];
    end

    % dec vars: V,h,b from S2; s1..s5,kappa from S1
    sol_x_it  = [sol2.x(1);     % V from S2
                 sol2.x(2);     % h from S2
                 sol1_x;        % s1,s2,s3,s4,s5,kappa from S1
                 sol2.x(end)];  % b from S2

    % for scaled termination
    delta_x     = sol_x_it - x_prev;
    inf_delta_x = norm(full(poly2basis(delta_x)),inf);

    % evaluate first order KKT at current iterate
    nabla_L(iter) = full(eval_gradLag( sol_x_it,lam_gs_it));

    inf_g   = full(eval_gradLag_inf_norm(sol_x_it,lam_gs_it));

    % just for plotting
    nabla_L_scaled(iter) = nabla_L(iter)*inf_delta_x/(max(1,abs(full(sol2.f)) + inf_g ));

    % scaled termination criterion
    if nabla_L(iter) <= commonOpts.tolerance_opt*(max(1,abs(full(sol2.f)) + inf_g ))/inf_delta_x
        break
    else
        % set current iterate to previous iterate
        x_prev =  sol_x_it;
    end
    end

    % if ~isempty(fval_old)
    %     if abs(full(sol2.f - fval_old)) <= commonOpts.tolerance_opt
    %         break
    %     else
    %         fval_old = sol2.f;
    %     end
    % else
    %     fval_old = sol2.f;
    % end

end % end for-loop (inner)

totalSolveTime = sum(solvetime_all1) + sum(solvetime_all2);


if iter == maxIter
    UNIFIED_RETURN_STATUS = 'RET_LIMITED';
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

stats.iterations            = iter;
stats.UNIFIED_RETURN_STATUS = UNIFIED_RETURN_STATUS;
stats.cost                  = cost_CD;
stats.nabla_L               = nabla_L;
stats.nabla_L_scaled        = nabla_L_scaled;
stats.initial_guess_trials  = initial_guess_trials;

stats.totalSolveTime = totalSolveTime;
stats.totalBuildTime = totalBuildTime;
stats.solver1_stats  = S1.stats;
stats.solver2_stats  = S2.stats;

% store final solution (same order as the sequential solver's decision
% variables, cf. sol_x_it in the main loop)
sol.x = sol_x_it;

solStruct.sol          = sol;
solStruct.stats        = stats;

% store information about subsolver e.g. for problem size
solStruct.inf0.solver1 =  S1.info;
solStruct.inf0.solver2 =  S2.info;


end
