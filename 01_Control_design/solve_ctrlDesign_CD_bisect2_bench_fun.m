function solStruct = solve_ctrlDesign_CD_bisect2_bench_fun(probStruct,commonOpts,x0)


x       = probStruct.x;
u       = probStruct.u;
f       = probStruct.f;
g       = probStruct.g;
s1      = probStruct.s1;
s2      = probStruct.s2;
s3      = probStruct.s3;
s4      = probStruct.s4;
V       = probStruct.V;
kappa   = probStruct.kappa;
b       = probStruct.b;


l = 1e-6*(x'*x);

umin = probStruct.umin;
umax = probStruct.umax;


%% setup solver
startBuildTime = tic;

% setup first solver
sos = struct('x',[s1;s2;s3;s4],... % decision variables
             'p',[V;kappa;b]);     % parameter

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

% solver for multiplier
S1 = casos.sossol('S1',commonOpts.solver,sos,opts);

% second solver


% setup second solver
cost = dot(g-(V-b),g-(V-b));

sos = struct('x',[V;kappa],...            % decision variables
    'f',1,...%cost,...                            % cost
    'p',[s1;s2;s3;s4;b]);  % parameter

% SOS constraints
sos.('g') = [
                V - l;
                s1*(V-b) - nabla(V,x)*subs(f,u,kappa) - l;
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

% setup solver: solveLyapunov function
S2 = casos.sossol('S2',commonOpts.solver,sos,opts);

totalBuildTime = toc(startBuildTime);


%% Initial guess

fval_old = [];

maxIter = commonOpts.maxiter;

% initial guess [V;kappa] (same ordering as in x0)
sol2_x = [x0(1);
          x0(6)];

fprintf('-------------------------------------------------- \n')
fprintf('Custom Coordinate Descent for Control Design \n')
fprintf('-------------------------------------------------- \n')


UNIFIED_RETURN_STATUS = 'SOLVER_RET_SUCCESS';


%% V-s-iteration
for iter = 1:maxIter


    %% s-step
    s1start = tic;


        % bisection tolerances
        relbistol = 1e-6;
        absbistol = 1e-6;
        lb = 0;
        ub = 100;
        while  (ub-lb>absbistol && ub-lb > relbistol*abs(lb))

            ptry = (ub+lb)/2;

            sol1 = S1('p',[sol2_x; ptry]);

            switch(S1.stats.UNIFIED_RETURN_STATUS)
                case {'SOLVER_RET_SUCCESS'}
                         lb = ptry;
                         sol1_x = sol1.x;
                case {'SOLVER_RET_INFEASIBLE','SOLVER_RET_NAN','SOLVER_RET_UNKNOWN','SOLVER_RET_LIMITED'}
                          ub = ptry;
                       
            end
        end

        % if strcmp(S1.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
        %     sol1_x = sol1.x;
        %     UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
        % else
        %     % if we just made a few iterations and become infeasible we
        %     % assume the initial guess is bad
        %     if iter <= 1
        % 
        %         sol1_x = x0(2:5);
        % 
        %         fprintf('Infeasible in solver 1 in iteration %d. We provide another initial guess. \n',iter)
        % 
        %     else
        %         fprintf('Infeasible in solver 1 in iteration %d; retaining previous iterate\n', iter);
        %         UNIFIED_RETURN_STATUS = S1.stats.UNIFIED_RETURN_STATUS;
        %         break
        %     end
        % end




    solvetime_all1(iter) =  toc(s1start);


    %% V-step
    S2start = tic;

    sol2 = S2('p',[sol1_x;lb]);

    cost_CD(iter) = full(sol2.f);

    solvetime_all2(iter) = toc(S2start);

    if strcmp(S2.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
        % extract solution
        sol2_x = sol2.x;
        UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
    else
        if iter <= 1
            fprintf('Infeasible in solver 2 in iteration %d. \n',iter)
            UNIFIED_RETURN_STATUS = S2.stats.UNIFIED_RETURN_STATUS;
            break
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


fprintf('----------------------------------- \n')
fprintf('Total solve time is %.2f s\n',totalSolveTime)
fprintf('Total build time is %.2f s\n',totalBuildTime)
fprintf('----------------------------------- \n')
fprintf('Total time is %.2f s\n',totalBuildTime+totalSolveTime)


if iter == maxIter
    UNIFIED_RETURN_STATUS = 'RET_LIMITED';
    stats.totalSolveTime = totalSolveTime;
end


% ************************************************************************
% check if coordinate-descent was infeasible after a few iterations
% Note: in this case we did not reach the termination criterion, but have a
% feasible solution

% Meaning: We made more than 1 iteration and one of the subsolver was
% infeasible we still return success
if iter > 1 && ~strcmp(UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS')
    UNIFIED_RETURN_STATUS= 'SOLVER_RET_SUCCESS';
end
% ************************************************************************

stats.iterations            = iter;
stats.UNIFIED_RETURN_STATUS = UNIFIED_RETURN_STATUS;
stats.cost                  = cost_CD;

stats.totalSolveTime = totalSolveTime;
stats.totalBuildTime = totalBuildTime;
stats.solver1_stats  = S1.stats;
stats.solver2_stats  = S2.stats;

% store final solution
sol.x.V     = sol2_x(1);
sol.x.kappa = sol2_x(2);
sol.x.b     = lb;
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
