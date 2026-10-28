
function [probStruct] = generateMiMo_probStruct(degV,degS)

x = casos.Indeterminates('x',3,1);
u = casos.Indeterminates('u',2,1);

f = [-x(1) + x(2)^2; -2*x(2) ; 3*x(3) + x(2)^3] + ...
    [-2; 1 + x(3)^2 ; 1 + 4*x(1)^2]*u(1) + ...
    [5*x(1); 1 - x(2)^2; 3]*u(2);


A0 = full(subs(nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

[K,P] = lqr(A0,B0,eye(3),eye(2));

fcl  = subs(f,u,-K*x);

% generate a desired safe set
g = x'*P*x-2;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,degV));

% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,degS));


% fill structure
probStruct.x    = x;        % indeterminates
probStruct.f    = fcl;      % dynamics
probStruct.g    = g;        % safe set
probStruct.s1   = s1;       % multiplier
probStruct.V    = V;        % Lyapunov function
probStruct.b    = 1;        % constant stable level set
probStruct.Name = 'MIMO';   % Name of example

end

 