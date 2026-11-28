% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for control design
%                       problem using a VDP dynamical system from [1].
%
%
%   Reference: [1]	T. Cunis und I. Kolmanovsky, „Viability, viscosity, and
%                   storage functions in model-predictive control with
%                   terminal constraints", Automatica, Bd. 131, S. 109748,
%                   Sep. 2021, doi: 10.1016/j.automatica.2021.109748.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function [probStruct] = generate_prob_struct_vdp_ctrl_design(degV,degK)

x = casos.Indeterminates('x',2,1);
u = casos.Indeterminates('u',1,1);


f = [x(2); 
      (1-x(1)^2)*x(2) - x(1) + u ];


umax =  1;
umin = -1;

% generate a desired safe set
g = x(1)^2 + 3*x(2)^2 - 1;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,degV));


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


probStruct.V    = V;        % Lyapunov function
probStruct.s1   = s1;       % multiplier
probStruct.s2   = s2;
probStruct.s3   = s3;
probStruct.s4   = s4;
probStruct.b    = b;
probStruct.umin = umin;
probStruct.umax = umax;

probStruct.kappa = kappa;

probStruct.Name = 'VDP'; % Name of example

end

 