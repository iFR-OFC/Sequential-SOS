
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

% bisection tolerances (see default value of SOSOPT/GSOSOPT
relbistol = 1e-3;
absbistol = 1e-3;
% lb = 0;
% ub = 100;
% while  (ub-lb>absbistol && ub-lb > relbistol*abs(lb))
% 
%     ptry = (ub+lb)/2;
% 
%     sol0 = S0('p',[Vval;ptry]);
% 
%     switch(S0.stats.UNIFIED_RETURN_STATUS)
%         case {'SOLVER_RET_SUCCESS'}
%                  lb = ptry;
%                  s1val  = sol0.x;
%                  break; % take the first feasible level set
%         case {'SOLVER_RET_INFEASIBLE','SOLVER_RET_NAN','SOLVER_RET_UNKNOWN'}
%                   ub = ptry;
%         otherwise 
%                error('Finding a feasible level set failed!')
%     end
% end
