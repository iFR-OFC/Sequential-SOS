function [solStruct] = solve_reachability_seq_bench_fun(probStruct,commonOpts,x0,hess_est_method)

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

% options
opts                              = struct('sossol',commonOpts.solver);
opts.max_iter                     = commonOpts.maxiter;

if  strcmpi(commonOpts.solver,'mosek')

    opts.sossol_options.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
    opts.sossol_options.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;

elseif strcmpi(commonOpts.solver,'clarabel')

    opts.sossol_options.sdpsol_options.clarabel = commonOpts.sdp_solver_opts;

elseif strcmpi(commonOpts.solver,'scs')

    opts.sossol_options.sdpsol_options.scs = commonOpts.sdp_solver_opts;

end

opts.hessian_approx               = hess_est_method;


% cost
cost = dot( g-(subs(V,t,0)-b), g-(subs(V,t,0)-b) ) ;

sos = struct('x',[V; kappa; s1; s2; s3; s4; s5; s6; s7; s8; s9; b],... % decision variables
    'f',cost, ...                                              % cost
    'p',[]);                                                   % parameter

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
             b - 1e-6;
             s1*(V-b) - s2*hT - nabla(V,t) - nabla(V,x)*subs(f,u,kappa);
             s3*(V-b) - s4*hT + kappa - umin;
             s5*(V-b) - s6*hT + umax - kappa;
             s7*(V-b) - s8*hT - g;
             s9*(subs(V,t,T)-b) - gT
            ];

% states + constraint cones
opts.Kx                           = struct('lin', length(sos.x));
opts.Kc                           = struct('sos', length(sos.g));
opts.sossol_options.newton_solver =[];
opts.tolerance_opt                = commonOpts.tolerance_opt;

% setup solver
S = casos.nlsossol('S','sequential',sos,opts);

% solve problem
sol = S('x0' ,x0);

% store output
solStruct.sol   = sol;
solStruct.stats = S.stats;
solStruct.info  = S.info;

end
