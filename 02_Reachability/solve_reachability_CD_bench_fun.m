function solStruct = solve_reachability_CD_bench_fun(probStruct,commonOpts,x0)

% de-pack for better readability
x       = probStruct.x;
u       = probStruct.u;
t       = probStruct.t;
f       = probStruct.f;
g       = probStruct.g;
s1      = probStruct.s1;
s2      = probStruct.s2;
s3      = probStruct.s3;
s4      = probStruct.s4;
s5      = probStruct.s5;
s6      = probStruct.s6;
s7      = probStruct.s7;
s8      = probStruct.s8;
s9      = probStruct.s9;
gT      = probStruct.gT;
T       = probStruct.T;
V       = probStruct.V;
kappa   = probStruct.kappa;
b       = probStruct.b;

% time polynomial
hT = (t)*(T-t);

umin = probStruct.umin;
umax = probStruct.umax;


%% setup solver
startBuildTime = tic;

% solver for control law and multiplier
sos = struct('x',[kappa; s1; s2; s3; s4; s5; s6; s7; s8; s9],... % decision variables
             'p',[V;b]);                                     % parameter

% constraints
sos.('g') = [
             s1;
             s2;
             s3;
             s4;
             s5;
             s6;
             s7;
             s8;
             s9;
             s1*(V-b) - s2*hT - nabla(V,t) - nabla(V,x)*subs(f,u,kappa);
             s3*(V-b) - s4*hT + kappa - umin;
             s5*(V-b) - s6*hT + umax - kappa;
             s7*(V-b) - s8*hT - g;
             s9*(subs(V,t,T)-b) - gT
            ];

% states + constraint are linear/SOS cones
opts.Kx           = struct('lin', length(sos.x));
opts.Kc           = struct('sos', length(sos.g));
opts.newton_solver = [];
opts.error_on_fail = 0;  % we perform our check during the V-s-iterations

% options for the selected SDP solver
if  strcmpi(commonOpts.solver,'mosek')

    opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
    opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

elseif strcmpi(commonOpts.solver,'clarabel')

    opts.sdpsol_options.clarabel = commonOpts.sdp_solver_opts;

elseif strcmpi(commonOpts.solver,'scs')

    opts.sdpsol_options.scs = commonOpts.sdp_solver_opts;

end

% build solver 1
S1 = casos.sossol('S',commonOpts.solver,sos,opts);

% build solver for reachability storage function and to improve solution
% quality

% minimize distance between V(t = 0, x) and the constraint set
cost = dot( g-(subs(V,t,0)-b), g-(subs(V,t,0)-b) );

sos = struct('x',[V;b],...            % decision variables
    'f',cost,...                     % cost
    'p',[kappa; s1; s2; s3; s4; s5; s6; s7; s8; s9]);  % solution from solver 1

% constraints
sos.('g') = [
             b - 1e-6;
             s1*(V-b) - s2*hT - nabla(V,t) - nabla(V,x)*(subs(f,u,kappa));
             s3*(V-b) - s4*hT + kappa - umin;
             s5*(V-b) - s6*hT + umax - kappa;
             s7*(V-b) - s8*hT - g;
             s9*(subs(V,t,T)-b) - gT
            ];

% states + constraint are linear/SOS cones
opts      = struct;
opts.Kx   = struct('lin', length(sos.x));
opts.Kc   = struct('sos', length(sos.g));
opts.newton_solver = [];
opts.error_on_fail = 0;  % we perform our check during the V-s-iterations

if  strcmpi(commonOpts.solver,'mosek')

    opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
    opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

elseif strcmpi(commonOpts.solver,'clarabel')

    opts.sdpsol_options.clarabel = commonOpts.sdp_solver_opts;

elseif strcmpi(commonOpts.solver,'scs')

    opts.sdpsol_options.scs = commonOpts.sdp_solver_opts;

end

% solver 2
S2 = casos.sossol('S',commonOpts.solver,sos,opts);

% functions for KKT condition

% common constraints for Lagrangian (same order as sequential solver)
g_common  = [
             s1;
             s2;
             s3;
             s4;
             s5;
             s6;
             s7;
             s8;
             s9;
             b - 1e-6;
             s1*(V-b) - s2*hT - nabla(V,t) - nabla(V,x)*subs(f,u,kappa);
             s3*(V-b) - s4*hT + kappa - umin;
             s5*(V-b) - s6*hT + umax - kappa;
             s7*(V-b) - s8*hT - g;
             s9*(subs(V,t,T)-b) - gT
            ];

% same order as initial guess for decision variables
sos_x     = [V;kappa;s1;s2;s3;s4;s5;s6;s7;s8;s9;b];

% get lagrange multiplier
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

startSolve_mainLoop = tic;


start_init_time = tic;
% initial guess [V;b]
sol2_x = [x0(1);
         x0(end)];

fprintf('------------------------------------------------------------------------------------------\n')
fprintf('Custom Coordinate Descent for Reachability Analysis\n')
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

useBisection = false;

% number of initial-guess restarts needed until the problem became feasible
initial_guess_trials = 1;

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
meas_time_kkt = 0;
solved = 1;

init_time = toc(start_init_time);

%% Main-loop (V-s-iteration)
% real k-th iteration runs at iter = k+1 (cf. sequential)

maxIter = commonOpts.maxiter;


for iter = 2:maxIter+1


    %% s-step
    s1start = tic;
    % second fall-back
    if useBisection
        startBuildBisec = tic;
        useBisection = false;
        s1_ok = false;

        % solver 1: s1-step
        sos0 = struct('x',[kappa; s1; s2; s3; s4; s5; s6;s7;s8;s9],... % decision variables
                      'f',-b,...
                      'p',V);    % parameter

        % SOS constraints (same as for first solver)
        sos0.('g') = [
                     s1;
                     s2;
                     s3;
                     s4;
                     s5;
                     s6;
                     s7;
                     s8;
                     s9;
                     s1*(V-b) - s2*hT - nabla(V,t) - nabla(V,x)*subs(f,u,kappa);
                     s3*(V-b) - s4*hT + kappa - umin;
                     s5*(V-b) - s6*hT + umax - kappa;
                     s7*(V-b) - s8*hT - g;
                     s9*(subs(V,t,T)-b) - gT
                     ];


        % states + constraint are SOS cones
        opts = [];
        opts.Kx = struct('lin', length(sos0.x));
        opts.Kc = struct('sos', length(sos0.g));
        opts.sossol = commonOpts.solver;
        opts.error_on_fail = 0;
        opts.sossol_options.newton_solver =[];

        % bisection fallback uses mosek: pass its params only if configured
        if  strcmpi(commonOpts.solver,'mosek')

            opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
            opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

        elseif strcmpi(commonOpts.solver,'clarabel')

            opts.sdpsol_options.clarabel = commonOpts.sdp_solver_opts;

        elseif strcmpi(commonOpts.solver,'scs')

            opts.sdpsol_options.scs = commonOpts.sdp_solver_opts;

        end

        % build bisection solver
        S0 = casos.qcsossol('S','bisection',sos0,opts);
        buildBisecTime = toc( startBuildBisec );

        % evaluate
        solveBisecTime = tic;
            sol0 = S0('p',casos.PD(sol2_x(1)));
        solveBisecTime = toc(solveBisecTime);

        switch(S0.stats.UNIFIED_RETURN_STATUS)
            case {'SOLVER_RET_SUCCESS'}

                sol1_x = sol0.x;

            otherwise
                UNIFIED_RETURN_STATUS = S0.stats.UNIFIED_RETURN_STATUS;
                break
        end

    else
        sol1 = S1('p',casos.PD(sol2_x(1:2)));

        s1_ok = strcmp(S1.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS');
        if s1_ok
            sol1_x = sol1.x;
            UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
        else
            % if we are infeasible in the first step, we try the same
            % initial guess for sequential i.e. for the multiplier
            if iter <= 2

                % take the same initial guess as sequential and try
                sol1_x = x0(2:end-1); % all except V,b

                % second try for initial guess
                initial_guess_trials = 2;

                fprintf('Infeasible in solver 1 in iteration %d. We provide another initial guess. \n',iter-1)

            else
                fprintf('Infeasible in solver 1 in iteration %d; retaining previous iterate\n', iter-1);
                UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
                break
            end
        end
    end

    solvetime_all1(iter) =  toc(s1start);

    % V-step
    S2start = tic;

    % take the solution from the multiplier step
    sol2 = S2('p',sol1_x);

    cost_CD(iter) = full(sol2.f);

    solvetime_all2(iter) = toc(S2start);

    if strcmp(S2.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
        % extract solution
        sol2_x = sol2.x;
    else
        % if the first two attempts for initial guess fail, try a bisection
        if iter <= 2
            fprintf('Infeasible in solver 2 in iteration %d. We try a third initial guess! \n',iter-1)
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
            useBisection = true;
            initial_guess_trials = 3;
        else
            fprintf('Infeasible in solver 2 in iteration %d; retaining previous iterate\n', iter-1);
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
        end

        % only break, in case we do not run the third initial guess
        if ~useBisection
            break
        end

    end

    % only check if we did not run the bisection to get a feasible initial
    % guess
    if ~useBisection

        eval_kkt_con_start = tic;

        % lam of s1..s9 (s3..s6 are nu-vectors) from S1, rest from S2
        if s1_ok
            % in same order as g_common (pairs correctly in dot(lam_gs_it,g_common))
            lam_gs_it = [sol1.lam_g(1:length([s1;s2;s3;s4;s5;s6;s7;s8;s9]));  % s1..s9 = g_common(1:k) from S1
                         sol2.lam_g];                                          % s10 + product terms = g_common(k+1:end) from S2
        else
            % in case S1 was infeasible, set multiplier to zero
            lam_gs_it = [zeros(length([s1;s2;s3;s4;s5;s6;s7;s8;s9]),1);        % s1..s9 from S1 (infeasible -> zero)
                         sol2.lam_g];                                          % s10 + product terms from S2
        end

        % dec vars: V,s10,b from S2; kappa,s1..s9 from S1
        sol_x_it  = [sol2.x(1);     % V from S2
                     sol1_x;        % kappa,s1..s9 from S1
                     sol2.x(2)];    % b from S2

        % for scaled termination
        delta_x     = sol_x_it - x_prev;
        inf_delta_x = norm(full(poly2basis(delta_x)),inf);

        % evaluate first order KKT at current iterate
        nabla_L(iter) = full(eval_gradLag( sol_x_it,lam_gs_it));

        % inf. norm of lambda*constraint(s)
        inf_g   = full(eval_gradLag_inf_norm(sol_x_it,lam_gs_it));

        % just for plotting
        nabla_L_scaled(iter) = nabla_L(iter)*inf_delta_x/(max(1,abs(full(sol2.f)) + inf_g ));

        meas_time_kkt(iter) = toc(eval_kkt_con_start);

        % display output header every 10th iteration
        if ~mod(iter-1,10) && iter-1 > 0
            fprintf('%-8s%-15s%-15s\n','iter','obj','||dLdx||');
            fprintf('------------------------------------------------------------------------------------------\n');
        end

        % display output current iterate
        fprintf('%-8d%-15g%-15g\n',iter-1,cost_CD(iter),nabla_L(iter));

        % scaled termination criterion
        if nabla_L(iter) <= commonOpts.tolerance_opt*(max(1,abs(full(sol2.f)) + inf_g ))/inf_delta_x
            break
        else
            % set current iterate to previous iterate
            x_prev =  sol_x_it;
        end

    end % check if we need bisection

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
% infeasible after a minimum number of iterations (commonOpts.miniter)
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

% we are feasible, but are not converged until max. iter
if strcmp(UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS') && iter >= maxIter
    % check for progress
    if  abs(cost_CD(1) - cost_CD(end)) > commonOpts.tolerance_opt

        UNIFIED_RETURN_STATUS = 'SOLVER_RET_SUCCESS';
        solved = 3;

    else
        UNIFIED_RETURN_STATUS = 'RET_LIMITED';
    end

end

% summary output
fprintf('------------------------------------------------------------------------------------------ \n')
if strcmp(UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS') && solved == 1
    fprintf('Solution status: Optimal solution found\n')
elseif solved == 2
    fprintf('Solution status: Abort midway. Solved to feasible point, but did not reach KKT point\n')
elseif solved == 3
    fprintf('Solution status: Maximum number of iterations reached. Solved, but needs further checks. \n')
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
stats.init_time          = init_time;
stats.eval_kkt_time      = sum(meas_time_kkt);
stats.time_per_iteration = stats.totalSolveTime / stats.iterations;

% write to struct
stats.solver1_stats      = S1.stats;
stats.solver2_stats      = S2.stats;

% individual solver times
stats.totalSolveTime_S1     =  sum(solvetime_all1);
stats.totalSolveTime_S2     =  sum(solvetime_all2);
stats.totalSolveTime_bisect = solveBisecTime;


% store final solution (same order as the sequential solver's decision
% variables, cf. sol_x_it in the main loop)
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
