%% ------------------------------------------------------------------------
%
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       region-of-attraction for the longitudinal motion
%                       of the Nasa Generic Transport Model. To increase
%                       the size of the sublevel set we try to minimize the
%                       squared distance to a defined set. Additionally, we
%                       synthesis a linear control law at the same time.
%                       State and control constraints are considered.
%
%   Reference: Susceptibility of F/A-18 Flight Controllers to the 
%   Falling-Leaf Mode: Nonlinear Analysis Abhijit Chakraborty, 
%   Peter Seiler, and Gary J. Balas Journal of Guidance, Control, and 
%   Dynamics 2011 34:1, 73-85
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%--------------------------------------------------------------------------
close all
clear
clc


%% system states
x = casos.Indeterminates('x',12,1);

import casos.toolboxes.sosopt.cleanpoly


f = gen12StateDynSys(x);

g = x'*eye(12)*0.1*x-1;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,2));

% enforce positivity
l = 1e-6*(x'*x);

% start to measure build time of all parameterized solver
buildTimes_start = tic;


b= casos.PS.sym('b');

b0 = 0.1;

%% setup solver
% solver 1: s1-step
sos0 = struct('x',s1,'p',[V;b]);     % parameter

% constraint
sos0.('g') = [s1; s1*(V-b)-nabla(V,x)*f-l];

% states + constraint are SOS cones
opts.Kx = struct('lin', 1);
opts.Kc = struct('sos', 2);
opts.error_on_fail = 0;
opts.newton_solver =[];
% build first solver
S0 = casos.sossol('S0','mosek',sos0,opts);
% solver 1: s1-step
sos1 = struct('x',s1,'p',V);     % parameter

% constraint
sos1.('g') = [s1; s1*(V-.1)-nabla(V,x)*f-l];

% states + constraint are SOS cones
opts.Kx = struct('lin', 1);
opts.Kc = struct('sos', 2);
opts.newton_solver =[];
opts.error_on_fail = 1;

% build first solver
S1 = casos.sossol('S1','mosek',sos1,opts);

% solver 2: V-step
sos2 = struct('x',V, ...        % dec.var
    'f',dot(g-(V-b0),g-(V-1)),...
    'p',s1);                % parameter

% constraints
sos2.('g') = [V-l;
              s1*(V-b0)-nabla(V,x)*f-l];

opts     = struct;
opts.Kx = struct('sos', 0, 'lin', 1);
opts.Kc = struct('sos', 2);
opts.newton_solver =[];
% build third solver
S2 = casos.sossol('S','mosek',sos2,opts);

totalBuildTime = toc(buildTimes_start);

A_cl = subs(nabla(f,x),x,zeros(12,1));

P = lyap(full(A_cl)',eye(12));

Vval = x'*P*x;

fval_old = [];

% bisection tolerances (see default value of SOSOPT/GSOSOPT
relbistol = 1e-3;
absbistol = 1e-3;


%% V-s-iteration
for iter = 1:100

    %% s1-step
    s1start = tic;

    % if iter == 1
    %     %% initial guess for feasible level set; simple bisection
    %     lb = 0;
    %     ub = 100;
    %     while  (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
    % 
    %         ptry = (ub+lb)/2;
    % 
    %         sol0 = S0('p',[Vval;ptry]);
    % 
    %         switch(S0.stats.UNIFIED_RETURN_STATUS)
    %             case {'SOLVER_RET_SUCCESS'}
    %                       lb = ptry;
    %                      s1val  = sol0.x;
    %                      break; % take the first feasible level set
    %             case {'SOLVER_RET_INFEASIBLE','SOLVER_RET_NAN','SOLVER_RET_UNKNOWN'}
    %                       ub = ptry;
    %             otherwise 
    %                    error('Finding a feasible level set failed!')
    %         end
    %     end
    % 
    %     % s1val = x'*x;
    % 
    % else
            sol1 = S1('p',Vval);
            s1val  = sol1.x;
    % end

    solvetime_all1(iter) =  toc(s1start);


    %% V-step
    S2start = tic;

    sol2 = S2('p',s1val);

    solvetime_all2(iter) = toc(S2start);

    % extract solution
    Vval = sol2.x(1);

    % after initial guess, leave the level set constant
    bval =1;

    % show progress
    fprintf('Iteration %d: f = %g, g = %g.\n',iter,full(sol2.f),1);


    if ~isempty(fval_old)
        if abs(full(sol2.f-fval_old)) <= 1e-4
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

