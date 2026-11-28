% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for CBF-CLF problem
%                       using a VDP dynamical system from [1].
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

function x0 = generate_init_guess_vdp_cbf_clf(probStruct)

f = probStruct.f;
x = probStruct.x;
u = probStruct.u;
g = probStruct.g;

A_0 = full(subs(nabla(f,x),[x;u],zeros(length(x)+length(u),1)));
B_0 = full(subs(nabla(f,u),[x;u],zeros(length(x)+length(u),1)));



[K,P] = lqr(A_0,B_0,eye(2),1);
% named initial guess for each decision variable
% (ordering: V,h,s1,s2,s3,s4,s5,kappa,b,gammaW,gammaV)

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