% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for CBF-CLF problem
%                       using an Aircraft (affine) dynamical system from
%                       [1].
%
%
%   Reference: [1]	A. Chakraborty, P. Seiler, und G. J. Balas, „Nonlinear
%                   region of attraction analysis for flight control
%                   verification and validation", Control Engineering
%                   Practice, Bd. 19, Nr. 4, S. 335–345, Apr. 2011, doi:
%                   10.1016/j.conengprac.2010.12.001.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_aircraft_affine_cbf_clf(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;
g = probStruct.g;

A0 = full(subs( nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

[K,P] = lqr(A0,B0,eye(4),1);

% initial guess
V0 = x'*P*x;
h0 = casos.PD(g);                 % start at the boundary of the safe set

s10 = x'*x;
s20 = x'*x;
s30 = ones(length(u),1);
s40 = ones(length(u),1);
s50 = x'*x;

kappa0 = -K*x;                    % LQR control law

b0 = 0.8;                         % level set of the terminal set

% gammaW0 = 0.01;
% gammaV0 = 0.01;

% initial guess for sequential / coordinate descent
x0 = casos.PD([V0; h0; s10; s20; s30; s40; s50; kappa0;b0]);



end