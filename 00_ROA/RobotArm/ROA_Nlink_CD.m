%--------------------------------------------------------------------------
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       region-of-attraction for the N-link robot arm using
%                       coordinate descent.
%
%   Reference:  T. Cunis and J. Olucak, 
%   "CaΣoS: A nonlinear sum-of-squares optimization suite⋆," 
%    2025 American Control Conference (ACC), Denver, CO, 
%    USA, 2025, pp. 1659-1666, doi: 10.23919/ACC63710.2025.11107794.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%--------------------------------------------------------------------------

function [fval_array, solverTimes_total,buildTimes,iteration_array,solverStats] = ROA_Nlink_CD(noRuns,deg,Nmax)


buildTimes          = zeros(noRuns,Nmax-1);
solverTimes_total   = zeros(noRuns,Nmax-1);

% compute each N-link noRuns times
for jj = 1:noRuns

    % iterate of the Number of links
    for n = 2:Nmax

        disp(['Compute maximum ROA for the ' num2str(n) '-link nLinkRobot'])

        x = casos.PS('x',2*n,1);

        % get dynamics
        f = feval(['nLinkRobot_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);

        % load P matrix
        load(['data_n' num2str(n)],'P');

        % Lyapunov function candidate
        Vval = x'*P*x;


        % enforce positivity
        l = 1e-6*(x'*x);

        % start to measure build time of all parameterized solver
        buildTimes_start = tic;

        % Lyapunov function candidate
        V = casos.PS.sym('v',monomials(x,2));

        % SOS multiplier
        s1 = casos.PS.sym('s1',monomials(x,2));


        % build a safe set
        g = x'*P*0.01*x-2;

        % fixed gamma
        gval = 1;

        %% setup solver

        % solver 1: s1-step
        sos1 = struct('x',s1,'p',V);     % parameter

        % constraint
        sos1.('g') = [s1; s1*(V-gval)-nabla(V,x)*f-l];

        % states + constraint are SOS cones
        opts.Kx = struct('lin', 1);
        opts.Kc = struct('sos', 2);
        opts.error_on_fail = 0;
        opts.newton_solver = [];
        % build first solver
        S1 = casos.sossol('S1','mosek',sos1,opts);

        % solver 3: V-step
        sos2 = struct('x',V, ...        % dec.var
            'f',dot(g-V,g-V),...
            'p',s1); % parameter

        % constraints
        sos2.('g') = [V-l;
            s1*(V-gval)-nabla(V,x)*f-l];

        opts     = struct;
        opts.error_on_fail = 0;
        opts.newton_solver = [];
        opts.Kx = struct('sos', 0, 'lin', 1);
        opts.Kc = struct('sos', 2);

        % build third solver
        S2 = casos.sossol('S','mosek',sos2,opts);

        tempBuildTime = toc(buildTimes_start);

        % initialize arrays
        solvetime_all1 = zeros(100,1);
        solvetime_all2 = zeros(100,1);
        solvetime_all3 = zeros(100,1);


        fval_old = [];
        S_step_success = 1;
        V_step_success = 1;

        startTotalTime = tic;
        %% V-s-iteration
        for iter = 1:100
               
            %% s1-step
            s1start = tic;

            % call solver
            sol1 = S1('p',Vval);

            % check solution status
            if ~strcmp(S1.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS' )
                fprintf('Infeasible in s-step in iteration %d\n', iter)
                S_step_success = 0;
                
                % give CD gratitude in case first iteration is infeasible
                if iter > 1
                    % If we are infasible after first iteration break
                    solvetime_all1(iter) =  toc(s1start);
                    break
                end
            end

            solvetime_all1(iter) =  toc(s1start);

            % extract solution
            s1val = sol1.x;

            %% V-step
            S2start = tic;

            sol2 = S2('p',s1val);

            % check solution status; no gratitude for V-step
            if ~strcmp(S2.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS' )
                fprintf('Infeasible in V-step in iteration %d ',iter)
                V_step_success = 0;
                solvetime_all3(iter) = toc(S2start);
                break
            end
            solvetime_all3(iter) = toc(S2start);

            % extract solution
            Vval = sol2.x;

            % show progress
            fprintf('Iteration %d: f = %g, g = %g.\n',iter,full(sol2.f),full(1));


            if ~isempty(fval_old)
                if abs(full(sol2.f-fval_old)) <= 1e-4
                    break
                else
                    fval_old = sol2.f;
                end
            else
                fval_old = sol2.f;
            end

            % check if we are above totaltime; if we are above 10hrs break
            if toc(startTotalTime) >= 36000 
                fprintf('Maximum total time reached in iteration %d ',iter)
                break
            end


        end % end for-loop (inner)

        % store the last beta-value
        if S_step_success  && V_step_success
            fval_array(jj,n-1)        = full(sol2.f);
            iteration_array(jj,n-1)   = iter;
            % solStatus_array{jj,n-1}   = S2.stats.solutionStatus;
            solverStats{jj,n-1}       = {{S1.info.sdpsol.conic,S2.info.sdpsol.conic}};
            % total solver time over all iterations
            solverTimes_total(jj,n-1) = sum(solvetime_all1) + sum(solvetime_all2) + sum(solvetime_all3);
            buildTimes(jj,n-1)        = tempBuildTime;% + sum(callTime1) + sum(callTime2) + sum(callTime3);
        else
            if ~isempty(fval_old)
                fval_array(jj,n-1)        = full(fval_old);
            else
                fval_array(jj,n-1)        = -1;
            end
            iteration_array(jj,n-1)   = iter;
            % solStatus_array{jj,n-1}   = S2.stats.solutionStatus;
            if S_step_success == 0
                solverStats{jj,n-1}       = {{S1.info.sdpsol.conic}};
            else
                solverStats{jj,n-1}       = {{S1.info.sdpsol.conic,S2.info.sdpsol.conic}};
            end
            % total solver time over all iterations
            solverTimes_total(jj,n-1) = sum(solvetime_all1) + sum(solvetime_all2) + sum(solvetime_all3);
            buildTimes(jj,n-1)        = tempBuildTime;% + sum(callTime1) + sum(callTime2) + sum(callTime3);

        end

    end % end of for loop for N-link
end % end of for loop for noRuns
end % end-of-function