
function [probStruct] = generateVDP_probStruct(degV,degS)

x = casos.Indeterminates('x',2,1);

f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)];


% generate a desired safe set
Q = eye(2);
g = x'*Q*x-1;

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
probStruct.Name = 'VDP'; % Name of example

end

 