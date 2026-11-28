% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for ROA problem using
%                       a VDP dynamical system from [1].
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

function [probStruct] = generate_prob_struct_vdp_roa(degV,degS)

x = casos.Indeterminates('x',2,1);

f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)];

% generate a desired safe set
Q = eye(2);
g = x'*Q*x-1;

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
probStruct.Name = 'VDP';    % Name of example

end
