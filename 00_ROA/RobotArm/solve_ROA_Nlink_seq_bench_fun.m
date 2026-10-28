%--------------------------------------------------------------------------
%
%   Short Description:  Calculate an inner-estimate of the
%                       region-of-attraction for the N-link robot arm
%                       using sequential SOS.
%
%                       Note: In contrast to the other ROA examples (see
%                       solve_ROA_seq_bench_fun.m), the squared distance of
%                       V itself to the safe set is minimized; the fixed
%                       level-set constant b enters the SOS constraint
%                       only. The structure, inputs, and outputs are
%                       identical to solve_ROA_seq_bench_fun.m.
%
%   Reference:  T. Cunis and J. Olucak,
%   "CaΣoS: A nonlinear sum-of-squares optimization suite⋆,"
%    2025 American Control Conference (ACC), Denver, CO,
%    USA, 2025, pp. 1659-1666, doi: 10.23919/ACC63710.2025.11107794.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%--------------------------------------------------------------------------

function [solStruct] = solve_ROA_Nlink_seq_bench_fun(probStruct,commonOpts,x0,hess_est_method)


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
opts.sossol_options.sdpsol_options.mosek_param = commonOpts.sdp_solver_opts.mosek_param;
opts.sossol_options.sdpsol_options.mosek_echo  = commonOpts.sdp_solver_opts.mosek_echo;
opts.hessian_approx               = hess_est_method;


% enforce positivity
l = 1e-6*(x'*x);

% cost: squared distance of V to the safe set (no level-set shift)
cost = dot(g-V,g-V) ;

sos = struct('x',[V; s1],... % decision variables
    'f',cost, ...            % cost
    'p',[]);                 % parameter

% SOS constraints
sos.('g') = [V-l;
    s1;
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
