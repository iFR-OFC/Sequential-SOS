function [solStruct] = solve_ctrlDesign_seq_bench_fun(probStruct,commonOpts,x0,hess_est_method)



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

% options
opts                              = struct('sossol','mosek');
opts.max_iter                     = commonOpts.maxiter;
opts.sossol_options.newton_solver = [];
opts.hessian_approx               = hess_est_method;


% cost
cost = dot(g-(V-b),g-(V-b)) ;

sos = struct('x',[V; s1;s2;s3;s4;kappa;b],... % decision variables
    'f',cost, ...                   % cost
    'p',[]);                        % parameter

% SOS constraints
sos.('g') = [s1;
    s2;
    s3;
    s4;
    V-l;
    b - 1e-6;
    s1*(V-b)-nabla(V,x)*subs(f,u,kappa)-l;
    s2*(V-b) + kappa - umin;
    s3*(V-b) + umax - kappa;
    s4*(V-b) - g
    ];

% states + constraint cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));
opts.sossol_options.newton_solver =[];

% setup solver
S = casos.nlsossol('S1','sequential',sos,opts);

% solve problem
sol = S('x0' ,x0);

solStruct.sol   = sol;
solStruct.stats = S.stats;
solStruct.info  =  S.info;


end