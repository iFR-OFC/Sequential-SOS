% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for the Power Converter
%                       dynamical system.
%
%
%   Reference:
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_power_conv_cbf_clf(probStruct)

f = probStruct.f;
x = probStruct.x;
u = probStruct.u;
g = probStruct.g;

A_0 = full(subs(nabla(f,x),[x;u],zeros(length(x)+length(u),1)));
B_0 = full(subs(nabla(f,u),[x;u],zeros(length(x)+length(u),1)));



[K,P] = lqr(A_0,B_0,eye(3),eye(2));


V0 = x'*P*x;  
h0 = casos.PD(g);                 % start at the boundary of the safe set

s10 = x'*x;
s20 = x'*x;
s30 = ones(length(u),1);
s40 = ones(length(u),1);
s50 = x'*x;

kappa0 = -K*x;                    % LQR control law

b0 = 0.1;                         % level set of the terminal set

% gammaW0 = 0.01;
% gammaV0 = 0.01;

% initial guess for sequential / coordinate descent
x0 = casos.PD([V0; h0; s10; s20; s30; s40; s50; kappa0;b0]);


end