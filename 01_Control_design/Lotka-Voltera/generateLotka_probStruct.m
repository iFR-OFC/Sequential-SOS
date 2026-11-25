
function [probStruct] = generateLotka_probStruct(degV,degK)

x = casos.Indeterminates('x',2,1);
u = casos.Indeterminates('u',1,1);


x1 =x(1)*1;
x2 =x(2)*1;

% System parameters
alpha = 0.1;
beta  = 0.2;
c1    = 0.01;
gamma = 0.04;
delta = 0.03;
c2    = 0.01 ;


f = [x1*alpha - x1*x2*beta  - x1*c1*u
    -x2*gamma + x2*x1*delta - x2*c2*u ];


umax =  0.5;
umin = -0.5;

x0 = [gamma/delta;alpha/beta];
u0 = 0;

f = subs(f,[x;u],[x+x0;u+u0]);




% generate a desired safe set
g = 1*x(1)^2 + 1*x(2)^2 - 1;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,degV));


% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,2));
s2    = casos.PS.sym('s2',monomials(x,2),[length(u) 1]);
s3    = casos.PS.sym('s3',monomials(x,2),[length(u) 1]);
s4    = casos.PS.sym('s4',monomials(x,2));
b     = casos.PS.sym('b');
kappa = casos.PS.sym('kappa',monomials(x,degK),[length(u) 1]);


% fill structure
probStruct.x    = x;        % indeterminates
probStruct.u    = u;
probStruct.f    = f;        % dynamics
probStruct.g    = g;        % safe set


probStruct.V    = V;        % Lyapunov function
probStruct.s1   = s1;       % multiplier
probStruct.s2   = s2;
probStruct.s3   = s3;
probStruct.s4   = s4;
probStruct.b    = b;
probStruct.umin = umin;
probStruct.umax = umax;

probStruct.kappa = kappa;

probStruct.Name = 'Lotka-Voltera'; % Name of example

end

 