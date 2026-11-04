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

% setup first solver (Prob 1): SOS multipliers, comparison rates and
% control law
% decision variables: [s1;s2;s3;s4;s5;gammaW;gammaV;kappa]
% parameters:         [V;h;b]
sos = struct('x',[s1;s2;s3;s4;s5;gammaW;gammaV;kappa],... % decision variables
    'p',[V;h;b]);    % parameter

% SOS constraints
sos.('g') = [
    s1;
    s2;
    s3;
    s4;
    s5;
    gammaW - 1e-6;
    gammaV - 1e-6;
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
% opts.sossol = commonOpts.solver;
% opts.tolerance_abs = 1e-6;
% opts.tolerance_rel = 1e-6;
opts.error_on_fail = 0; % we do our own check in the main loop
opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

% solver for multipliers, gammas and control law
S1 = casos.sossol('S1','mosek',sos,opts);


% setup second solver (Prob 2): CLF, CBF and level set
% decision variables: [V;h;b]
% parameters:         [s1;s2;s3;s4;s5;gammaW;gammaV;kappa]
cost = dot(g-(h-b),g-(h-b));

sos = struct('x',[V;h;b],...            % decision variables
    'f',cost,...                            % cost
    'p',[s1;s2;s3;s4;s5;gammaW;gammaV;kappa]);  % parameter

% SOS constraints
sos.('g') = [
    V - l;                                                  % CLF positivity
    b - 1e-6;                                               % keep the level set positive
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
opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

% solver for CLF, CBF and level set
S2 = casos.sossol('S2',commonOpts.solver,sos,opts);

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


UNIFIED_RETURN_STATUS = 'SOLVER_RET_SUCCESS';

% number of initial-guess restarts needed until the problem became feasible
initial_guess_trials = 1;

%% V-h-b-s-iteration
for iter = 1:maxIter


    %% s-step (Prob 1)
    s1start = tic;

    sol1 = S1('p',sol2_x);

    if strcmp(S1.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
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

    if strcmp(S2.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
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
stats.initial_guess_trials  = initial_guess_trials;

stats.totalSolveTime = totalSolveTime;
stats.totalBuildTime = totalBuildTime;
stats.solver1_stats  = S1.stats;
stats.solver2_stats  = S2.stats;

% store final solution
sol.x.V      = sol2_x(1);
sol.x.h      = sol2_x(2);
sol.x.b      = -sol1.f;
sol.x.s1     = sol1_x(1);
sol.x.s2     = sol1_x(2);
sol.x.s3     = sol1_x(3:5);
sol.x.s4     = sol1_x(6:8);
sol.x.s5     = sol1_x(9);
sol.x.gammaW = sol1_x(10);
sol.x.gammaV = sol1_x(11);
sol.x.kappa  = sol1_x(12:14);

solStruct.sol          = sol;
solStruct.stats        = stats;

% store information about subsolver e.g. for problem size
solStruct.inf0.solver1 =  S1.info;
solStruct.inf0.solver2 =  S2.info;


end
