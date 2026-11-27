
function [probStruct] = generateMIMO_probStruct(degV,degK)


x = casos.Indeterminates('x',3,1);
u = casos.Indeterminates('u',2,1);

% time
t  = casos.PS('t');

f0 = [-x(1) + x(2)^2; -2*x(2) ; 3*x(3) + x(2)^3];
g1 = [-2; 1 + x(3)^2 ; 1 + 4*x(1)^2];
g2 = [5*x(1); 1 - x(2)^2; 3];

f = f0 + g1*u(1) +  g2*u(2);


A0 = full(subs(nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

[~,P] = lqr(A0,B0,eye(3),eye(2));

umax =  1;
umin = -1;

% generate a desired safe set
g  = 3*x(2)^2 + x(1)^2 -1;  % g(x) \leq 0
gT = x'*P*x-1;


b  = casos.PS.sym('b');


V  = casos.PS.sym('v',monomials([x;t],degV));

kappa = casos.PS.sym('kappa',monomials([x;t],degK),[length(u) 1]);

s1 = casos.PS.sym('s1',monomials([x;t],0:2));
s2 = casos.PS.sym('s2',monomials([x;t],0:2));
s3 = casos.PS.sym('s3',monomials([x;t],0:2),[length(u) 1]);
s4 = casos.PS.sym('s4',monomials([x;t],0:2),[length(u) 1]);
s5 = casos.PS.sym('s5',monomials([x;t],0:2),[length(u) 1]);
s6 = casos.PS.sym('s6',monomials([x;t],0:2),[length(u) 1]);
s7 = casos.PS.sym('s7',monomials([x;t],0:2));
s8 = casos.PS.sym('s8',monomials([x;t],0:2));
s9 = casos.PS.sym('s9',monomials(x,0:2));

% Time horizon
T  = 1;

% fill structure
probStruct.x    = x;        % indeterminates
probStruct.u    = u;
probStruct.t    = t;

probStruct.f    = f;        % dynamics

probStruct.g    = g;        % constraint set
probStruct.gT   = gT;       % Terminal set

probStruct.T    = T;      % time horizon

probStruct.umin = umin; % control bounds
probStruct.umax = umax;

probStruct.V    = V;        % Reachability storage function

probStruct.s1   = s1;       % multiplier
probStruct.s2   = s2;
probStruct.s3   = s3;
probStruct.s4   = s4;
probStruct.s5   = s5;
probStruct.s6   = s6;
probStruct.s7   = s7;
probStruct.s8   = s8;
probStruct.s9   = s9;


probStruct.b    = b;  % level set

probStruct.kappa = kappa; % control law






probStruct.Name = 'VDP'; % Name of example

end

 