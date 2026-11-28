% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for the VDP-coupled
%                       dynamical system.
%
%
%   Reference:
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function [probStruct] = generate_prob_struct_vdp_coupled_cbf_clf(degV,degH,degK)

x = casos.Indeterminates('x',4,1);
u = casos.Indeterminates('u',2,1);

x1 = x(1);
y1 = x(2);
x2 = x(3);
y2 = x(4);

f = [y1;
    y1*(1 - x1^2) - x1 + x2 - x1 + u(1);
    y2;
    y2*(1 - x2^2) - x2 + x1 - x2+ u(2)];


% A0 = full(subs(nabla(f,x),[x;u],zeros(length(x)+length(u),1)));
% B0 = full(subs(nabla(f,u),[x;u],zeros(length(x)+length(u),1)));
% 
% K = lqr(A0,B0,eye(4),eye(2));

umin = [-1 -1]';
umax = -umin;

% f = subs(f,u,-K*x);

% generate a desired safe set
Q = eye(4);
g = x'*Q*x-1;

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

probStruct.Name = 'VDP-coupled'; % Name of example

end


 