function [solStruct] = solve_barrierCert_seq_bench_fun(probStruct,commonOpts,x0,hess_est_method)



x  = probStruct.x;
u  = probStruct.u;
f  = probStruct.f;
g  = probStruct.g;
s1 = probStruct.s1;
s2 = probStruct.s2;
s3 = probStruct.s3;
s4 = probStruct.s4;
B  = probStruct.B;
kappa = probStruct.kappa;
b     = probStruct.b;

l = 1e-6*(x'*x);

umin = probStruct.umin;
umax = probStruct.umax;

% options
opts                              = struct('sossol',commonOpts.solver);
opts.max_iter                     = commonOpts.maxiter;
opts.sossol_options.newton_solver = [];
opts.sossol_options.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
opts.sossol_options.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;
opts.hessian_approx               = hess_est_method;


% cost
cost = dot(g-(B-b),g-(B-b)) ;

sos = struct('x',[B; s1;s2;s3;s4;kappa;b],... % decision variables
    'f',cost, ...                   % cost
    'p',[]);                        % parameter

% SOS constraints
sos.('g') = [
            s2;
            s3;
            s4;
            s1*(B-b)-nabla(B,x)*subs(f,u,kappa);
            s2*(B-b) + kappa - umin;
            s3*(B-b) + umax - kappa;
            s4*(B-b) - g
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