% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for the 12 State
%                       Dynamical System.
%
%
%   Reference:
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function [probStruct] = generate_prob_struct_12state_roa()

import casos.toolboxes.sosopt.cleanpoly

x = casos.Indeterminates('x',12,1);

f = gen12StateDynamics(x);

g = x'*eye(12)*0.1*x-1;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,2));

% fill structure
probStruct.x    = x;        % indeterminates
probStruct.f    = f;        % dynamics
probStruct.g    = g;        % safe set
probStruct.s1   = s1;       % multiplier
probStruct.V    = V;        % Lyapunov function
probStruct.b    = 0.1;      % constant stable level set
probStruct.Name = '12 State Dynamical System'; % Name of example

end
