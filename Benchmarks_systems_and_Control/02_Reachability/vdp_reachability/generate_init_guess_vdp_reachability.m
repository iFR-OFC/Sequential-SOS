% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for reachability
%                       problem using a VDP dynamical system from [1].
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

function x0 = generate_init_guess_vdp_reachability(probStruct)


x = probStruct.x;
u = probStruct.u;
t = probStruct.t;
f = probStruct.f;
g = probStruct.g;

A0 = full(subs( nabla(f,x),[x;u],zeros(3,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(3,1)));

[K,P] = lqr(A0,B0,eye(2),1);




V0      =    x'*P*x; 
kappa0  =    -K*x; 
s10     =    (x'*x) + t^2; 
s20     =    x'*x + t^2; 
s30     =    (x'*x) + t^2;
s40     =    x'*x + t^2;
s50     =    (x'*x) + t^2;
s60     =    x'*x + t^2;
s70     =    (x'*x) + t^2;
s80     =    x'*x + t^2;
s90     =    (x'*x);
b0      =    0.8;

% initial guess
x0 = casos.PD([ V0;     ...
                kappa0 ; ...
                s10    ;...
                s20    ;...
                s30    ;...
                s40    ;...
                s50    ;...
                s60    ;...
                s70    ;...
                s80    ;...
                s90    ;...
                b0]);


end