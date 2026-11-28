% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for control design
%                       problem using a Lotka-Voltera dynamical system from
%                       [1].
%
%
%   Reference: [1]	A. Ibañez, „Optimal control of the Lotka–Volterra
%                   system: turnpike property and numerical simulations",
%                   Journal of Biological Dynamics, Bd. 11, Nr. 1, S.
%                   25–41, Jan. 2017, doi: 10.1080/17513758.2016.1226435.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_lotka_voltera_ctrl_design(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


A0 = full(subs( nabla(f,x),[x;u],zeros(3,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(3,1)));

[K,P] = lqr(A0,B0,eye(2),1);

% initial guess
V0 = x'*P*x;
s10 = (x'*x);
s20 = x'*x;
s30 = x'*x;
s40 = x'*x;
kappa0  = -K*x;
b0    = 1;


x0 = casos.PD([V0;s10;s20;s30;s40;kappa0;b0]);


end