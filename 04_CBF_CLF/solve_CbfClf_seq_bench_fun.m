function [solStruct] = solve_CbfClf_seq_bench_fun(probStruct,commonOpts,x0,hess_est_method)



x  = probStruct.x;
u  = probStruct.u;
f  = probStruct.f;
g  = probStruct.g;
V  = probStruct.V;
h  = probStruct.h;
kappa = probStruct.kappa;
s1 = probStruct.s1;
s2 = probStruct.s2;
s3 = probStruct.s3;
s4 = probStruct.s4;
s5 = probStruct.s5;
gammaV = probStruct.gammaV;
gammaW = probStruct.gammaW;
b  = probStruct.b;

umin = probStruct.umin;
umax = probStruct.umax;

% options
opts                                           = struct('sossol',commonOpts.solver);
opts.max_iter                                  = commonOpts.maxiter;
opts.con_vio_check                             = 'signed-distance';
opts.sossol_options.newton_solver              = [];
opts.sossol_options.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
opts.sossol_options.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;
opts.hessian_approx                            = hess_est_method;


% cost
cost = dot(g-(h-b),g-(h-b)) ;

% decision variables (sorted: V,h,s1,s2,s3,s4,s5,kappa,b,gammaW,gammaV)
sos = struct('x',[V;h;s1;s2;s3;s4;s5;kappa;gammaW;gammaV;b],...
    'f',cost, ...                   % cost
    'p',[]);                        % parameter

% SOS constraints
sos.('g') = [
            s1;
            s2;
            s3;
            s4;
            s5;
            gammaW - 1e-6;
            gammaV - 1e-6;
            V - 1e-6*(x'*x);                                      % CLF positivity
            b - 1e-6;                                             % keep the level set positive
            s1*(h-b) - g;                                         % state constraints
            s2*(h-b) - nabla(h,x)*subs(f,u,kappa) - gammaW*(h-b); % CBF dissipation
            s3*(h-b) + kappa - umin;                              % control constraints
            s4*(h-b) + umax - kappa;
            s5*(h-b) - nabla(V,x)*subs(f,u,kappa) - gammaV*V     % CLF dissipation
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
