% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for barrier
%                       certificate problem using an Aircraft (affine)
%                       dynamical system from [1].
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

function x0 = generate_init_guess_aircraft_affine_barrier_cert(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


A0 = full(subs( nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

[K,P] = lqr(A0,B0,eye(4),1);

% initial guess
B0 = x'*P*x;
s10 = (x'*x)^2;
s20 = x'*x;
s30 = x'*x;
s40 = x'*x;
kappa0  = -K*x;
b0  = 1;


x0 = casos.PD([B0;s10;s20;s30;s40;kappa0;b0]);


end