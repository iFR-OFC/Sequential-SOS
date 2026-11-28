% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for ROA problem using
%                       a VDP-coupled dynamical system from [1, p. 136].
%
%
%   Reference: [1]	X. Chen, „Reachability Analysis of Non-Linear Hybrid
%                   Systems Using Taylor Models", RWTH Aachen, Aachen,
%                   2015. [Online]. URL:
%                   https://publications.rwth-aachen.de/record/465295/files
%                   /465295.pdf
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function [probStruct] = generate_prob_struct_vdp_coupled_roa(degV,degS)

x = casos.Indeterminates('x',4,1);
u = casos.Indeterminates('u',2,1);

x1 = x(1);
y1 = x(2);
x2 = x(3);
y2 = x(4);

f = [-y1;
    y1*(1 - x1^2) - x1 + x2 - x1 + u(1);
    -y2;
    y2*(1 - x2^2) - x2 + x1 - x2+ u(2)];

A0 = full(subs(nabla(f,x),[x;u],zeros(length(x)+length(u),1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(length(x)+length(u),1)));

K = lqr(A0,B0,eye(4),eye(2));

f = subs(f,u,-K*x);

% generate a desired safe set
Q = eye(4);
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
probStruct.Name = 'VDP-coupled'; % Name of example

end
