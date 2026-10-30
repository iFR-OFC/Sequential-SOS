
function [probStruct] = generatePendubot_probStruct(degB,degK)

x = casos.Indeterminates('x',4,1);
u = casos.Indeterminates('u',1,1);


f = pendubot_dynamic(x,u);

umax =  1;
umin = -1;


theta1max = 30*pi/180;
theta1dotmax = 1*pi/180;
theta2max = 30*pi/180;
theta2dotmax = 1*pi/180;


Dx = diag([1/(theta1max*2),1/(theta1dotmax*2),1/(theta2max*2),1/(theta2dotmax*2)]);
Dx = eye(4);

f = Dx*subs(f,x,Dx\x);

% generate a desired safe set
n = 2;
g = (x(1)^2/theta1max^2)^(n/2) + (x(2)^2/theta1dotmax^2)^(n/2) + (x(3)^2/theta2max^2)^(n/2) + ...
     (x(4)^2/theta2dotmax^2)^(n/2) - 1;

g = subs(g,x,Dx\x);



g =  x'*diag([1/0.1^2, 1/0.35^2, 1/0.1^2, 1/0.35^2])*x-1;
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

probStruct.Name = 'Pendubot'; % Name of example

end

 