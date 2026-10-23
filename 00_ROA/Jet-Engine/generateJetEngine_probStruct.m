
function [probStruct] = generateJetEngine_probStruct(degV,degS)

x = casos.Indeterminates('x',2,1);
u = casos.Indeterminates('u',2,1);

f0  = dynamics_JetEngine(x,u);

A0 = full(subs(nabla(f0,x),[x;u],zeros(length(x)+length(u),1)));
B0 = full(subs(nabla(f0,u),[x;u],zeros(length(x)+length(u),1)));


K = lqr(A0,B0,eye(2),eye(2));

f = subs(f0,u,-K*x);

% generate a desired safe set
Q = eye(2);
g = x'*Q*x-1;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,degV));

% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,degS));


% fill structure
probStruct.x    = x;        % indeterminates
probStruct.u    = u;        % indeterminates
probStruct.f    = f;        % dynamics
probStruct.g    = g;        % safe set
probStruct.s1   = s1;       % multiplier
probStruct.V    = V;        % Lyapunov function
probStruct.b    = 0.1;        % constant stable level set

probStruct.Name = 'Jet-Engine';    % Name of example

end

 