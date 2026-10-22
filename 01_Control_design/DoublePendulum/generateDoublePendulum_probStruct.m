
function [probStruct] = generateDoublePendulum_probStruct(degV,degS)

x = casos.Indeterminates('x',4,1);
u = casos.Indeterminates('u',2,1);

% generate closed-loop dynamics
f = doublePendulumDynamics(x,u);


% generate a desired safe set as superellipsoid
n = 2;
g = (x(1)^2/1^2)^(n/2) + (x(2)^2/2^2)^(n/2) + (x(3)^2/1^2)^(n/2) + ...
     (x(4)^2/2^2)^(n/2) - 1;


% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,degV));

% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,degS));



% fill structure
probStruct.x    = x;        % indeterminates
probStruct.f    = f;        % dynamics
probStruct.g    = g;        % safe set
probStruct.s1   = s1;       % multiplier
probStruct.V    = V;        % Lyapunov function
probStruct.b    = 1;        % constant stable level set
probStruct.Name = 'Double-pendulum'; % Name of example

end

 