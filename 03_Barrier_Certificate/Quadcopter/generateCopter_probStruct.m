
function [probStruct] = generateCopter_probStruct(degB,degK)

x = casos.Indeterminates('x',6,1);
u = casos.Indeterminates('u',2,1);


f = copter_dynamics(x,u);

gn = 9.81;
K0 = 0.89/1.4;

umax =  [+1.5 + gn*K0; pi/12];
umin =  [-1.5 + gn*K0; -pi/12];

% generate a desired safe set
n = 2;
g = (x(1)^2/1.7^2)^(n/2) + (x(2)^2/0.85^2)^(n/2) + (x(3)^2/0.8^2)^(n/2) + ...
     (x(4)^2/1^2)^(n/2) + (x(5)^2/0.2618^2)^(n/2) + (x(6)^2/1.5708^2)^(n/2) - 1;

% Lyapunov function candidate
B = casos.PS.sym('v',monomials(x,degB));


% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,2:4));
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


probStruct.B    = B;        % Lyapunov function
probStruct.s1   = s1;       % multiplier
probStruct.s2   = s2;
probStruct.s3   = s3;
probStruct.s4   = s4;
probStruct.b    = b;
probStruct.umin = umin;
probStruct.umax = umax;

probStruct.kappa = kappa;

probStruct.Name = 'Copter'; % Name of example

end

 