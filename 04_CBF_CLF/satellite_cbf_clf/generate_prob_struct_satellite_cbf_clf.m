% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for CBF-CLF problem
%                       using a Satellite dynamical system from [1].
%
%
%   Reference: [1]	J. L. Junkins und H. Schaub, Analytical Mechanics of
%                   Space Systems, Second Edition. Reston, VA: American
%                   Institute of Aeronautics and Astronautics, 2009. doi:
%                   10.2514/4.867231.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function [probStruct] = generate_prob_struct_satellite_cbf_clf(degV,degH,degK)

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


% control constraint; assumption is that the box is inside the full
% torque volume. This is roughly estimated visually.
umin = [-1 -1 -1]'*1.2;
umax = [ 1  1  1]'*1.2;

Dx   = diag([1/(x_up(1)-x_low(1)),1/(x_up(2)-x_low(2)),1/(x_up(3)-x_low(3)),1,1,1]);


%% dynamics
% cross-product matrix
cpm = @(x) [   0  -x(3)  x(2);
              x(3)   0   -x(1);
             -x(2)  x(1)   0 ];

% dynamics
B = @(sigma) (1-sigma'*sigma)*eye(3)+ 2*cpm(sigma)+ 2*sigma*sigma';

f_raw =  [-J\cpm(x(1:3))*J*x(1:3) + J\u;
      1/4*B(x(4:6))*x(1:3)]; % omega_dot

% scale dynamics
f = Dx*subs(f_raw,[x;u],[Dx\x;u]);

% state constraint
n = 2;
g0 = (x(1)^2/omegaMax1^2)^(n/2) + (x(2)^2/omegaMax2^2)^(n/2) + (x(3)^2/omegaMax3^2)^(n/2) + ...
     (x(4)^2/0.57^2)^(n/2) + (x(5)^2/0.57^2)^(n/2) + (x(6)^2/0.57^2)^(n/2) - 1;

% re-scale input of state constraints
g = subs(g0,x,Dx\x);

%% SOS decision variables

% CLF (we used adjusted monomial here i.e, we removed very small terms
% found by the initial guess)
V  = casos.PS.sym('v',monomials(x,degV));

% CBF
h  = casos.PS.sym('h',monomials(x,degH));

% control law (common polynomial vector)
kappa = casos.PS.sym('kappa',monomials(x,degK),[length(u) 1]);

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,2));
s2 = casos.PS.sym('s2',monomials(x,2));
s3 = casos.PS.sym('s3',monomials(x,2),[length(u) 1]);
s4 = casos.PS.sym('s4',monomials(x,2),[length(u) 1]);
s5 = casos.PS.sym('s5',monomials(x,2));

% rates of the comparison functions
gammaV = 0.001; %casos.PS.sym('gammaV');
gammaW = 0.001; %casos.PS.sym('gammaW');

% level set of terminal set
b = casos.PS.sym('b');


% fill structure
probStruct.x     = x;        % indeterminates
probStruct.u     = u;
probStruct.f     = f;        % dynamics (scaled)
probStruct.f_raw = f_raw;    % dynamics (unscaled, used for the initial guess)
probStruct.Dx    = Dx;       % state scaling
probStruct.g     = g;        % safe set

probStruct.V     = V;        % CLF
probStruct.h     = h;        % CBF
probStruct.kappa = kappa;    % control law

probStruct.s1    = s1;       % SOS multiplier
probStruct.s2    = s2;
probStruct.s3    = s3;
probStruct.s4    = s4;
probStruct.s5    = s5;

probStruct.gammaV = gammaV;
probStruct.gammaW = gammaW;

probStruct.b     = b;

probStruct.umin  = umin;
probStruct.umax  = umax;

probStruct.Name  = 'Satellite'; % Name of example

end
