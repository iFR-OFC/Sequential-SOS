% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for barrier
%                       certificate problem using a Power Converter
%                       dynamical system [1].
%
%
%   Reference: [1]	M. Schneeberger, F. Dörfler, und S. Mastellone, „SOS
%                   Construction of Compatible Control Lyapunov and Barrier
%                   Functions", IFAC-PapersOnLine, Bd. 56, Nr. 2, S.
%                   10428–10434, Jan. 2023, doi:
%                   10.1016/j.ifacol.2023.10.1058.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_power_conv_barrier_cert(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


A0 = full(subs( nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

K = lqr(A0,B0,eye(3),eye(2));

% initial guess
B0      = probStruct.g;
s10     = (x'*x);
s20     = ones(length(u),1)*(x'*x);
s30     = ones(length(u),1)*x'*x;
s40     = x'*x;
kappa0  = -K*x;
b0      = 0.1;


x0 = casos.PD([B0;s10;s20;s30;s40;kappa0;b0]);


end