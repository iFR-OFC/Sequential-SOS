function solStruct = solve_reachability_CD_bench_fun(probStruct,commonOpts,x0)



x  = probStruct.x;
u  = probStruct.u;
t = probStruct.t;
f  = probStruct.f;
g  = probStruct.g;
s1 = probStruct.s1;
s2 = probStruct.s2;
s3 = probStruct.s3;
s4 = probStruct.s4;
s5 = probStruct.s5;
s6 = probStruct.s6;
s7 = probStruct.s7;
s8 = probStruct.s8;
s9 = probStruct.s9;
s10 = probStruct.s10;
gT  = probStruct.gT ;
T = probStruct.T;

V  = probStruct.V;
kappa = probStruct.kappa;
b     = probStruct.b;

% time polynomial
hT = (t)*(T-t);

umin = probStruct.umin;
umax = probStruct.umax;

% options for mosek
opts.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
opts.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

%% setup solver
startBuildTime = tic;

% solver for control law and multiplier
sos = struct('x',[kappa; s1; s2; s3; s4; s5; s6;s7;s8;s9],...
              'p',[V;b]);


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
             s9*(subs(V,t,T)-b)    - gT;
             ];

% states + constraint are linear/SOS cones
opts.Kx = struct('lin', length(sos.x));
opts.Kc = struct('sos', length(sos.g));
opts.error_on_fail = 0;  % we perform our check during the V-s-iterations
opts.newton_solver = [];

% build solver 1
S1  = casos.sossol('S',commonOpts.solver,sos,opts);


% build solver for reachability storage function and to improve solution
% quality

% minimize distance between V(t = 0, x) and the constraint set
cost = dot( g-(subs(V,t,0)-b), g-(subs(V,t,0)-b) );

opts = [];
sos = struct('x',[V;b;s10],...
              'f', cost , ...
              'p',[kappa; s1; s2; s3; s4; s5; s6;s7;s8;s9]); %solution from solver 1


% constraints
sos.('g') = [
             s10
             s1*(V-b) - s2*hT - nabla(V,t) - nabla(V,x)*(subs(f,u,kappa));
             s3*(V-b) - s4*hT + kappa - umin;
             s5*(V-b) - s6*hT + umax - kappa;
             s7*(V-b) - s8*hT - g; 
             s9*(subs(V,t,T)-b)    - gT;
             s10*gT  - subs(V,t,0) + b;
             ];

% states + constraint are linear/SOS cones
opts.Kx = struct('lin', length(sos.x));
opts.Kc = struct('sos', length(sos.g));

opts.error_on_fail = 0;  % we perform our check during the V-s-iterations
opts.newton_solver = [];

% solver 2
S2  = casos.sossol('S','mosek',sos,opts);

% common constraints for Lagrangian (same order as sequential solver)
g_common  = [s1;
             s2;
             s3;
             s4;
             s5;
             s6;
             s7;
             s8;
             s9;
             s10;
             s1*(V-b) - s2*hT - nabla(V,t) - nabla(V,x)*subs(f,u,kappa);
             s3*(V-b) - s4*hT + kappa - umin;
             s5*(V-b) - s6*hT + umax - kappa;
             s7*(V-b) - s8*hT - g;
             s9*(subs(V,t,T)-b) - gT;
             s10*gT - subs(V,t,0) + b];

% same order as initial guess for decision variables
sos_x     = [V;kappa;s1;s2;s3;s4;s5;s6;s7;s8;s9;s10;b];

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

% initial guess [V;b;s10]
sol2_x = [
          x0(1);
          x0(end)
          x0(end-1);
          ];

fprintf('-------------------------------------------------- \n')
fprintf('Custom Coordinate Descent \n')
fprintf('-------------------------------------------------- \n')

% x_prev for termination criterion; take initial guess as first "previous"
% iteration
x_prev = x0;


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
        s1_ok = false;
        % solver 1: s1-step
        sos0 = struct('x',[kappa; s1; s2; s3; s4; s5; s6;s7;s8;s9],... % decision variables
            'f',-b,...
            'p',V);    % parameter
        
        % SOS constraints ( same as for first solver)
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
                     s9*(subs(V,t,T)-b)    - gT;
                     ];


        % states + constraint are SOS cones
        opts = [];
        opts.Kx = struct('lin', length(sos0.x));
        opts.Kc = struct('sos', length(sos0.g));
        opts.sossol = 'mosek';
        opts.error_on_fail = 0;
        % opts.conf_interval 
        opts.tolerance_abs = 1e-6;
        opts.tolerance_rel = 1e-6;


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
  
                otherwise
                    UNIFIED_RETURN_STATUS = S0.stats.UNIFIED_RETURN_STATUS;
                    break
            end
 
    else
        sol1 = S1('p',casos.PD(sol2_x(1:2)));

        s1_ok = strcmp(S1.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS');
        if strcmp(S1.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
            sol1_x = sol1.x;
            UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
        else
            % if we are infeasible in the first step, we try the same
            % initial guess for sequential i.e. for th multiplier
            if iter <= 1
                
                % take the same initial guess as sequential and try
                sol1_x = x0(2:end-2); % all except V,b and s10
                
                % second try for initial guess
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
    
    % take the solution from the multiplier step 
    sol2 = S2('p',sol1_x);

    cost_CD(iter) = full(sol2.f);

    solvetime_all2(iter) = toc(S2start);

    if strcmp(S2.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
        % extract solution
        sol2_x = sol2.x;
    else
        % if the first two attempts for initial guess fail, try a bisection
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
                     sol2.x(3);     % s10 from S2
                     sol2.x(2)];    % b from S2

        % for scaled termination
        delta_x     = sol_x_it - x_prev;
        inf_delta_x = norm(full(poly2basis(delta_x)),inf);

        % evaluate first order KKT at current iterate
        nabla_L = full(eval_gradLag( sol_x_it,lam_gs_it));

        inf_g   = full(eval_gradLag_inf_norm(sol_x_it,lam_gs_it));

        % scaled termination criterion
        if nabla_L < commonOpts.tolerance_opt*(max(1,abs(full(sol2.f)) + inf_g ))/inf_delta_x
            break
        else
            % set current iterate to previous iterate
            x_prev =  sol_x_it;
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
% infeasible after a minimum number of iterations (commonOpts.miniter)
% and the cost function made some progress over these iterations;
% otherwise the run is marked as RET_LIMITED

if ~strcmpi(UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS')
    if iter < maxIter
        if iter >= commonOpts.miniter && ... % we made a least a few iterations
           abs(cost_CD(1) - cost_CD(end)) > commonOpts.tolerance_opt % some progress
            UNIFIED_RETURN_STATUS = 'SOLVER_RET_SUCCESS';
        else
            % means we got stuck
            UNIFIED_RETURN_STATUS = 'RET_LIMITED';
        end
    end
end

% summary output
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

sol.x = [sol2_x(1); sol1_x; sol2_x(2);sol2_x(3)];

solStruct.sol          = sol;
solStruct.stats        = stats;

% store information about the last subsolver e.g. for problem size
solStruct.inf0.solver1 =  S1.info;
solStruct.inf0.solver2 =  S2.info;


end
