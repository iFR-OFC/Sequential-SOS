% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for control design
%                       problem using a Pendubot dynamical system from [1].
%
%
%   Reference: [1]	H. Yin, M. Arcak, A. Packard, und P. Seiler, „Backward
%                   Reachability for Polynomial Systems on a Finite
%                   Horizon", IEEE Transactions on Automatic Control, Bd.
%                   66, Nr. 12, S. 6025–6032, Dez. 2021, doi:
%                   10.1109/TAC.2021.3056611.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_pendubot_ctrl_design(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


A0 = full(subs( nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

[K,P] = lqr(A0,B0,eye(4),1);

% initial guess
V0 = x'*P*x;
s10 = (x'*x);
s20 = x'*x;
s30 = x'*x;
s40 = x'*x;
kappa0  = -K*x;
b0  = 1;


x0 = casos.PD([V0;s10;s20;s30;s40;kappa0;b0]);


end