
function [probStruct] = generateSatellite_probStruct(degV,degS)

% system states
x = casos.PS('x',6,1);
u = casos.PS('u',3,1);

%% Hubble telescope parameter
J = diag([31046;77217;78754]);

% simple bounds on rates;
omegaMax1 = 0.5*pi/180;
omegaMax2 = 0.2*pi/180;
omegaMax3 = 0.2*pi/180;

x_low =  [-omegaMax1 -omegaMax2 -omegaMax3]';
x_up  =  [ omegaMax1  omegaMax2  omegaMax3]';

% scaling matrix for system states
Dx   = diag([1/(x_up(1)-x_low(1)),1/(x_up(2)-x_low(2)),1/(x_up(3)-x_low(3)),1,1,1]); 

%% system dynamics

% cross-product matrix
cpm = @(x) [   0  -x(3)  x(2); 
    x(3)   0   -x(1); 
    -x(2)  x(1)   0 ];

% MRP derivative
B = @(sigma) (1-sigma'*sigma)*eye(3)+ 2*cpm(sigma)+ 2*sigma*sigma';

% dynamics
f =  [-J\cpm(x(1:3))*J*x(1:3) + J\u; % omega_dot
    1/4*B(x(4:6))*x(1:3)];         % sigma_dot

% trim point/equilibrium point
x0    = [0 0 0 0 0 0]';
u0    = [0,0,0]';

A0 = full(casos.PD(subs(nabla(f,x),[x;u],[x0;u0])));
B0 = full(casos.PD(subs(nabla(f,u),[x;u],[x0;u0])));

% cost function weights
Q = diag([1, 1, 1, 20, 20 ,20]);
R = eye(3)*1;

% generate an initial guess for CLF
[K0,~] = lqr(full(A0),full(B0),Q,R);

f = subs(f,u,-K0*x);

% scale dynamics
f = casos.PD(Dx*subs(f,x,Dx\x));

% state constraint as super-ellipsoid
% Remark: The MRP parameter are selected such that \sigma'*\sigma <= 1
n = 2;
g0 = (x(1)^2/omegaMax1^2)^(n/2) + (x(2)^2/omegaMax2^2)^(n/2) + (x(3)^2/omegaMax3^2)^(n/2) + ...
    (x(4)^2/0.57^2)^(n/2) + (x(5)^2/0.57^2)^(n/2) + (x(6)^2/0.57^2)^(n/2) - 1;

% re-scale input of state constraints
g = subs(g0,x,Dx\x); 


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

probStruct.Name = 'Satellite-Attitude'; % Name of example

end

 