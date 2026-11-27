
function [probStruct] = generateSatellite_probStruct(degV,degK)

% system states
x = casos.Indeterminates('x',6,1);
u = casos.Indeterminates('u',3,1);
t = casos.Indeterminates('t');

% Hubble telescope parameter
J = diag([31046;77217;78754]);

% simple bounds on rates;
omegaMax1 = 0.5*pi/180;
omegaMax2 = 0.2*pi/180;
omegaMax3 = 0.2*pi/180;

x_low =  [-omegaMax1 -omegaMax2 -omegaMax3]';
x_up  =  [ omegaMax1  omegaMax2  omegaMax3]';


% control constraint; assumption is that the box is inside the full
% torque volume. This is roughly estimated visually.
umin = [-1 -1 -1]'*1.2;
umax = [ 1  1  1]'*1.2;

Dx   = diag([1/(x_up(1)-x_low(1)),1/(x_up(2)-x_low(2)),1/(x_up(3)-x_low(3)),1,1,1]); 


% dynamics
% cross-product matrix
cpm = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];

% dynamics
B = @(sigma) (1-sigma'*sigma)*eye(3)+ 2*cpm(sigma)+ 2*sigma*sigma';

f =  [-J\cpm(x(1:3))*J*x(1:3) + J\u;
      1/4*B(x(4:6))*x(1:3)]; % omega_dot

f = casos.PD(Dx*subs(f,x,Dx\x));

% state constraint
n = 2;
g0 = (x(1)^2/omegaMax1^2)^(n/2) + (x(2)^2/omegaMax2^2)^(n/2) + (x(3)^2/omegaMax3^2)^(n/2) + ...
     (x(4)^2/0.57^2)^(n/2) + (x(5)^2/0.57^2)^(n/2) + (x(6)^2/0.57^2)^(n/2) - 1;

% re-scale input of state constraints
g = subs(g0,x,Dx\x); 

A0 = full(subs( nabla(f,x),[x;u],zeros(9,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(9,1)));

[~,P] = lqr(A0,B0,eye(6),eye(3));

gT = x'*P*x - 1; 


V  = casos.PS.sym('v',monomials([x;t],degV));

kappa = casos.PS.sym('kappa',monomials([x;t],degK),[length(u) 1]);

s1 = casos.PS.sym('s1',monomials([x;t],0:2));
s2 = casos.PS.sym('s2',monomials([x;t],0:2));
s3 = casos.PS.sym('s3',monomials([x;t],0:2),[length(u) 1]);
s4 = casos.PS.sym('s4',monomials([x;t],0:2),[length(u) 1]);
s5 = casos.PS.sym('s5',monomials([x;t],0:2),[length(u) 1]);
s6 = casos.PS.sym('s6',monomials([x;t],0:2),[length(u) 1]);
s7 = casos.PS.sym('s7',monomials([x;t],0:2));
s8 = casos.PS.sym('s8',monomials([x;t],0:2));
s9 = casos.PS.sym('s9',monomials(x,0:2));
b  = casos.PS.sym('b');

% Time horizon
T  = 2;

% fill structure
probStruct.x    = x;        % indeterminates
probStruct.u    = u;
probStruct.t    = 1*t;
probStruct.Dx   = Dx;
probStruct.f    = f;        % dynamics

probStruct.g    = g;        % constraint set
probStruct.gT   = gT;       % Terminal set

probStruct.T    = T;      % time horizon

probStruct.umin = umin; % control bounds
probStruct.umax = umax;

probStruct.V    = V;        % Reachability storage function

probStruct.s1   = s1;       % multiplier
probStruct.s2   = s2;
probStruct.s3   = s3;
probStruct.s4   = s4;
probStruct.s5   = s5;
probStruct.s6   = s6;
probStruct.s7   = s7;
probStruct.s8   = s8;
probStruct.s9   = s9;

probStruct.b    = b;  % level set
probStruct.kappa = kappa; % control law

probStruct.Name = 'Satellite (Attitude)'; % Name of example

end

 