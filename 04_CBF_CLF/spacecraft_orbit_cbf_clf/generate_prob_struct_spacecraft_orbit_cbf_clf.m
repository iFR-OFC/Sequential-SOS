% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for CBF-CLF problem
%                       using a Spacecraft (orbit) dynamical system from
%                       [1].
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

function [probStruct] = generate_prob_struct_spacecraft_orbit_cbf_clf(degV,degH,degK)

% system states
x = casos.PS('x',6,1);
u = casos.PS('u',3,1);

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
probStruct.f     = f;        % dynamics 
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

probStruct.Name  = 'Spacecraft (orbit)'; % Name of example

end
