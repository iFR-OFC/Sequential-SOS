% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for the 12 State
%                       Dynamical System.
%
%
%   Reference:
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_12state_roa(probStruct)

f = probStruct.f;
x = probStruct.x;

A_cl = subs(nabla(f,x),x,zeros(12,1));

P = lyap(full(A_cl)',eye(12));

Vval = x'*P*x;

s10 = x'*x;

% initial guess
V0 = Vval;
x0 = casos.PD([V0;s10]);

end
