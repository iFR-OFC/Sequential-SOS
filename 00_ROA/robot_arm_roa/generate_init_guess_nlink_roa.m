% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for the N-link
%                       robot arm dynamical system.
%
%
%   Reference:
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_nlink_roa(probStruct)

x = probStruct.x;
P = probStruct.P;

% Lyapunov function candidate from the LQR solution
V0 = x'*P*x;

% unit polynomial as initial multiplier
s10 = (x'*x);

% initial guess
x0 = casos.PD([V0;s10]);

end
