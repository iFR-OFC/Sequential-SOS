
function [probStruct] = generateLotkaVol_probStruct(degV,degS)

x = casos.Indeterminates('x',5,1);

f = dynamics_LotkaVolterra5(x);


% generate a desired safe set

A_cl = subs(nabla(f,x),x,zeros(length(x),1));

P = lyap(full(A_cl)',eye(length(x)));

g = x'*P*x - 1;

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

probStruct.Name = 'Lotka-Voltera'; % Name of example

end

 