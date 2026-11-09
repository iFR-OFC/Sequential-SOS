
function [probStruct] = generateCopter_probStruct(degV,degK)

x = casos.Indeterminates('x',6,1);
u = casos.Indeterminates('u',2,1);


f = copter_dynamics(x,u);

gn = 9.81;
K0 = 0.89/1.4;

u0 =  gn/K0;

umax =  [+1.5 + u0; pi/12];
umin =  [-1.5 + u0; -pi/12];


xmax = [1.7;0.85;0.8;1;pi/12;pi/2];
xmin = -xmax;


Dx = diag([1./(xmax(1)-xmin(1)),1./(xmax(2)-xmin(2)),1./(xmax(3)-xmin(3)),...
           1./(xmax(4)-xmin(4)),1./(xmax(5)-xmin(5)),1./(xmax(6)-xmin(6))]);

Du = diag([1./(umax(1)-umin(1));1./(umax(2)-umin(2))]);


umax = Du*(umax)+u0;
umin = Du*(umin)+u0;

f = Dx*subs(f,[x;u],[inv(Dx)*x;inv(Du)*(1*u-u0)]);

% generate a desired safe set
n = 2;
g0 = (x(1)^2/1.7^2)^(n/2) + (x(2)^2/0.85^2)^(n/2) + (x(3)^2/0.8^2)^(n/2) + ...
     (x(4)^2/1^2)^(n/2) + (x(5)^2/(pi/12)^2)^(n/2) + (x(6)^2/(pi/2)^2)^(n/2) - 1;


g = subs(g0,x,inv(Dx)*x);

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,degV));


% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,2:4));
s2    = casos.PS.sym('s2',monomials(x,2),[length(u) 1]);
s3    = casos.PS.sym('s3',monomials(x,2),[length(u) 1]);
s4    = casos.PS.sym('s4',monomials(x,2:4));
b     = casos.PS.sym('b');
kappa = casos.PS.sym('kappa',monomials(x,degK),[length(u) 1]);


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

probStruct.Name = 'Copter'; % Name of example

end

 