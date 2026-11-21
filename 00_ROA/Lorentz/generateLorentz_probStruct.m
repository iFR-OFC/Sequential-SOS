
function [probStruct] = generateLorentz_probStruct(degV,degS)

X = casos.Indeterminates('x',3,1);

x = 1*X(1);
y = 1*X(2);
z = 1*X(3);


sigma = 10;
rho   = 0.9;
beta  = 8/3;

f = [sigma*(y - x);
     x*(rho - z) - y; 
     x*y - beta*z];


% generate a desired safe set
A_cl = subs(nabla(f,X),X,zeros(length(X),1));

P = lyap(full(A_cl)',eye(length(X)));
g = x'*P*x-1;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,degV));

% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,degS));


% fill structure
probStruct.x    = X;        % indeterminates
probStruct.f    = f;        % dynamics
probStruct.g    = g;        % safe set
probStruct.s1   = s1;       % multiplier
probStruct.V    = V;        % Lyapunov function
probStruct.b    = 0.1;        % constant stable level set
probStruct.Name = 'Lorentz'; % Name of example

end
