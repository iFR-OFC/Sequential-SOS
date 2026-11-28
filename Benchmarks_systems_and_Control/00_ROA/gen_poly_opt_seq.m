clear
clc

x = casos.Indeterminates('x',3);
u = casos.PS.sym('u',2);



%% setup solver

% options
opts                              = struct('sossol','mosek');



s1 = casos.PS.sym('s',monomials(x,0:2));


h = x'*x*u(1) - u'*u - 1;
g = x'*x+u'*u - 1;


c = [1;1];
P = eye(2);

% cost
cost = 1/2*u'*P*u + c'*u ;

sos = struct('x',[s1;u],... % decision variables
    'f',cost, ...            % cost
    'p',[]);                 % parameter

% SOS constraints
sos.('g') = [
    s1;
    s1*h - g
    ];

% states + constraint cones
opts.Kx                           = struct('lin', length(sos.x));
opts.Kc                           = struct('sos', length(sos.g));
opts.sossol_options.newton_solver =[];

% setup solver
S = casos.nlsossol('S1','sequential',sos,opts);

% solve problem

x0 = [x'*x; ones(2,1)];

sol = S('x0' ,x0);
sol.f
sol.x