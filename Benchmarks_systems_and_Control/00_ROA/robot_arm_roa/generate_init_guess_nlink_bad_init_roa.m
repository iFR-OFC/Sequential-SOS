% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for ROA problem using
%                       a N-link robot arm dynamical system.
%
%
%   Reference:
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_nlink_bad_init_roa(probStruct)

x = probStruct.x;

% clearly infeasible initial guess: negative unit polynomials for both the
% Lyapunov function and the SOS multiplier (cf. ROA_Nlink_sequential_badInit)
V0  = -(x'*x);
s10 = -(x'*x);

% initial guess
x0 = casos.PD([V0;s10]);

end
