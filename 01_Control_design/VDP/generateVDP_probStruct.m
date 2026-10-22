
function [probStruct] = generateVDP_probStruct(degV,degK)

x = casos.Indeterminates('x',2,1);
u = casos.Indeterminates('u',1,1);


f = [-x(2); 
      u*1 - x(1) + (x(1)^2 - 1)*x(2)];


umax = 1;
umin = -1;
% generate a desired safe set
Q = eye(2)*0.1;
g = x'*Q*x-1;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,degV));


% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,2:4));
s2    = casos.PS.sym('s2',monomials(x,2),[length(u) 1]);
s3    = casos.PS.sym('s3',monomials(x,2),[length(u) 1]);
s4    = casos.PS.sym('s4',monomials(x,2));
b     = casos.PS.sym('b');
kappa = casos.PS.sym('kappa',monomials(x,degK));


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

probStruct.Name = 'VDP'; % Name of example

end

 