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


%% setup solver
% options
opts = struct('sossol','mosek');
opts.verbose  = 1;
opts.almostOptCount = 10;
opts.sossol_options.newton_solver = [];

b = 0.1;

% cost
cost = dot(g-(V-b),g-(V-b)) ;

sos = struct('x',[V; s1],... % decision variables
    'f',cost, ...              % cost
    'p',[]);                   % parameter

% SOS constraints
sos.('g') = [s1;
    V-l;
    s1*(V-b)-nabla(V,x)*f-l;
    ];

% states + constraint cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));
opts.sossol_options.newton_solver =[];

% setup solver
S = casos.nlsossol('S1','sequential',sos,opts);

%% solve problem


A_cl = subs(nabla(f,x),x,zeros(12,1));

P = lyap(full(A_cl)',eye(12));

Vval = x'*P*x;

% initial guess
V0 = Vval;
s20 = (x'*x);
% b0  = 1;
x0 = casos.PD([V0;s20]);

% solve problem
sol = S('x0' ,x0);

% casos.postProcessSolver(S,true);

S.stats
