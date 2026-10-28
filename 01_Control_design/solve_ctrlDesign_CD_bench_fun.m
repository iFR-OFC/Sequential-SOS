function solStruct = solve_ctrlDesign_CD_bench_fun(probStruct,commonOpts,x0)


x  = probStruct.x;
u  = probStruct.u;
f  = probStruct.f;
g  = probStruct.g;
s1 = probStruct.s1;
s2 = probStruct.s2;
s3 = probStruct.s3;
s4 = probStruct.s4;
V  = probStruct.V;
kappa = probStruct.kappa;
b     = probStruct.b;


l = 1e-6*(x'*x);

umin = probStruct.umin;
umax = probStruct.umax;


%% setup solver
startBuildTime = tic;

% setup first solver
sos = struct('x',[s1;s2;s3;s4;kappa],... % decision variables
    'p',[V;b]);    % parameter

% SOS constraints
sos.('g') = [
    s1;
    s2;
    s3;
    s4;
    s1*(V-b) - nabla(V,x)*subs(f,u,kappa);
    s2*(V-b) + kappa - umin;
    s3*(V-b) + umax - kappa;
    s4*(V-b) - g
    ];

% states + constraint cones
opts.Kx           = struct('lin', length(sos.x));
opts.Kc           = struct('sos', length(sos.g));
opts.newton_solver = [];
opts.error_on_fail = 0; % we do our own check in the main loop
opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

% solver for multiplier
S1 = casos.sossol('S1',commonOpts.solver,sos,opts);


% setup second solver
cost = dot(g-(V-b),g-(V-b));

sos = struct('x',[V;b],...            % decision variables
    'f',cost,...                            % cost
    'p',[s1;s2;s3;s4;kappa]);  % parameter

% SOS constraints
sos.('g') = [V-l;
            b - 1e-6;
            s1*(V-b) - nabla(V,x)*subs(f,u,kappa);
            s2*(V-b) + kappa - umin;
            s3*(V-b) + umax - kappa;
            s4*(V-b) - g
            ];

% states + constraint cones
opts      = struct;
opts.Kx   = struct('lin', length(sos.x));
opts.Kc   = struct('sos', length(sos.g));
opts.newton_solver = [];
opts.error_on_fail = 0; % we do our own check in the main loop
opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

% setup solver: solveLyapunov function
S2 = casos.sossol('S2',commonOpts.solver,sos,opts);

totalBuildTime = toc(startBuildTime);


%% Initial guess

fval_old = [];

maxIter = commonOpts.maxiter;

% initial guess [V;kappa;b] (same ordering as in x0)
sol2_x = [x0(1);
          x0(end)];

fprintf('-------------------------------------------------- \n')
fprintf('Custom Coordinate Descent for Control Design \n')
fprintf('-------------------------------------------------- \n')


UNIFIED_RETURN_STATUS = 'SOLVER_RET_SUCCESS';

useBisection = false;

% number of initial-guess restarts needed until the problem became feasible
initial_guess_trials = 1;

%% V-s-iteration
for iter = 1:maxIter


    %% s-step
    s1start = tic;
    % second fall-back
    if useBisection 

        useBisection = false;
        % solver 1: s1-step
        sos0 = struct('x',[s1;s2;s3;s4;kappa],... % decision variables
            'f',-b,...
            'p',V);    % parameter
        
        % SOS constraints
        sos0.('g') = [
            s1;
            s2;
            s3;
            s4;
            s1*(V-b) - nabla(V,x)*subs(f,u,kappa);
            s2*(V-b) + kappa - umin;
            s3*(V-b) + umax - kappa;
            s4*(V-b) - g
            ];


        % states + constraint are SOS cones
        opts = [];
        opts.Kx = struct('lin', length(sos0.x));
        opts.Kc = struct('sos', length(sos0.g));
        opts.sossol = 'mosek';
        opts.error_on_fail = 0;
        opts.sossol_options.newton_solver =[];
        opts.sossol_options.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
        opts.sossol_options.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

        % build first solver

        % solve by bisection
        S0 = casos.qcsossol('S','bisection',sos0,opts);
        % evaluate
        sol0 = S0('p',casos.PD(sol2_x(1)));

            switch(S0.stats.UNIFIED_RETURN_STATUS)
                case {'SOLVER_RET_SUCCESS'}
                   
                         sol1_x = sol0.x;
                         % break; % take the first feasible level set
                otherwise
                    UNIFIED_RETURN_STATUS = S0.stats.UNIFIED_RETURN_STATUS;
                    break
            end
 
    else
        sol1 = S1('p',sol2_x);

        if strcmp(S1.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
            sol1_x = sol1.x;
            UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
        else
            % if we just made a few iterations and become infeasible we
            % assume the initial guess is bad
            if iter <= 1

                sol1_x = x0(2:end-1);

                initial_guess_trials = 2;

                fprintf('Infeasible in solver 1 in iteration %d. We provide another initial guess. \n',iter)

            else
                fprintf('Infeasible in solver 1 in iteration %d; retaining previous iterate\n', iter);
                UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
                break
            end
        end
    end




    solvetime_all1(iter) =  toc(s1start);


    %% V-step
    S2start = tic;

    sol2 = S2('p',sol1_x);

    cost_CD(iter) = full(sol2.f);

    solvetime_all2(iter) = toc(S2start);

    if strcmp(S2.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
        % extract solution
        sol2_x = sol2.x;
    else
        if iter <= 1
            fprintf('Infeasible in solver 2 in iteration %d. We try a third initial guess! \n',iter)
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
            useBisection = true;
            initial_guess_trials = 3;
        else
            fprintf('Infeasible in solver 2 in iteration %d; retaining previous iterate\n', iter);
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
        end

        % only break, in case we do not run the third initial guess
        if ~useBisection
        break
        end

    end


    if ~useBisection
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

stats.iterations            = iter;
stats.UNIFIED_RETURN_STATUS = UNIFIED_RETURN_STATUS;
stats.cost                  = cost_CD;
stats.initial_guess_trials  = initial_guess_trials;

stats.totalSolveTime = totalSolveTime;
stats.totalBuildTime = totalBuildTime;
stats.solver1_stats  = S1.stats;
stats.solver2_stats  = S2.stats;

% store final solution
sol.x.V     = sol2_x(1);
sol.x.kappa = sol1_x(end);
sol.x.b     = sol2_x(1);
sol.x.s1    = sol1_x(1);
sol.x.s2    = sol1_x(2);
sol.x.s3    = sol1_x(3);
sol.x.s4    = sol1_x(4);

solStruct.sol          = sol;
solStruct.stats        = stats;

% store information about subsolver e.g. for problem size
solStruct.inf0.solver1 =  S1.info;
solStruct.inf0.solver2 =  S2.info;


end
