
function [solStruct] = solve_ROA_seq_bench_fun(probStruct,commonOpts,x0,hess_est_method)



x  = probStruct.x;
f  = probStruct.f;
g  = probStruct.g;
s1 = probStruct.s1;
V  = probStruct.V;
b  = probStruct.b;

%% setup solver

% options
opts                              = struct('sossol',commonOpts.solver);
opts.max_iter                     = commonOpts.maxiter;
opts.sossol_options.newton_solver = [];
opts.hessian_approx               = hess_est_method;


% enforce positivity
l = 1e-6*(x'*x);

% cost
cost = dot(g-(V-b),g-(V-b)) ;

sos = struct('x',[V; s1],... % decision variables
    'f',cost, ...            % cost
    'p',[]);                 % parameter

% SOS constraints
sos.('g') = [s1;
    V-l;
    s1*(V-b)-nabla(V,x)*f-l;
    ];

% states + constraint cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));

% setup solver
S = casos.nlsossol('S1','sequential',sos,opts);

%% solve problem

% solve problem
sol = S('x0' ,x0);

solStruct.sol   = sol;
solStruct.stats = S.stats;
solStruct.info  =  S.info;


end


