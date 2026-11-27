
function [probStruct] = generateSpacecraft_orbit_probStruct(degB,degK)

% system states
x = casos.PS('x',6,1);
u = casos.PS('u',3,1);

% acceleration
umin = [-1 -1 -1]'*0.5;
umax = [ 1  1  1]'*0.5;

% dynamics
n = 0.0013;


A = [zeros(3,3),eye(3);
     3*n^2 0 0 0 2*n 0;
     0 0 0 -2*n 0 0;
     0 0 -n^2 0 0 0 ];


B = [zeros(3,3);eye(3)];
f =  A*x+B*u; % omega_dot


% state constraint
n = 2;
g0 = (x(1)^2/100^2)^(n/2) + (x(2)^2/100^2)^(n/2) + (x(3)^2/100^2)^(n/2) + ...
     (x(4)^2/1^2)^(n/2) + (x(5)^2/1^2)^(n/2) + (x(6)^2/1^2)^(n/2) - 1;

% re-scale input of state constraints
g = subs(g0,x,x); 


% Lyapunov function candidate
B = casos.PS.sym('v',monomials(x,degB));


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


probStruct.B    = B;        % Lyapunov function
probStruct.s1   = s1;       % multiplier
probStruct.s2   = s2;
probStruct.s3   = s3;
probStruct.s4   = s4;
probStruct.b    = b;

probStruct.umin = umin;
probStruct.umax = umax;

probStruct.kappa = kappa;

probStruct.Name = 'Satellite'; % Name of example

end

 