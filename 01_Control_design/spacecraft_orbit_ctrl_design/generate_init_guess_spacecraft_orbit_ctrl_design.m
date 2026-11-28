% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for the
%                       Spacecraft_orbit dynamical system.
%
%
%   Reference:
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_spacecraft_orbit_ctrl_design(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


% trim point
x0    = [0 0 0 0 0 0]';
u0    = [0,0,0]';

A0 = full(casos.PD(subs(nabla(f,x),[x;u],[x0;u0])));
B0 = full(casos.PD(subs(nabla(f,u),[x;u],[x0;u0])));

% cost function weights
Q = diag([1, 1, 1, 1, 1 ,1]);
R = eye(3)*1;

[K,P] = lqr(A0,B0,Q,R);

% initial guess
V0 = x'*P*x;
s10 = (x'*x);
s20 = ones(length(u),1);
s30 = ones(length(u),1);
s40 = x'*x;
kappa0  = -K*x;
b0  = 0.1;


x0 = casos.PD([V0;s10;s20;s30;s40;kappa0;b0]);

end