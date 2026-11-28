% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for the Aircraft
%                       dynamical system.
%
%
%   Reference:
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_aircraft_roa(probStruct)

f = probStruct.f;
x = probStruct.x;

A_cl = subs(nabla(f,x),x,zeros(4,1));

P = lyap(full(A_cl)',eye(4)*0.1);

Vval = x'*P*x;

% initial guess
V0 = Vval;
s10 = (x'*x);

x0 = casos.PD([V0;s10]);

end
