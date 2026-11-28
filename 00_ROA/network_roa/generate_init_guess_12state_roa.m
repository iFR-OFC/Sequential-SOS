% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for ROA problem using
%                       a 12 State Dynamical System dynamical system [1].
%
%
%   Reference: [1]	J. Anderson und A. Papachristodoulou, „Dynamical 
%                   system decomposition for efficient, sparse analysis", 
%                   in 49th IEEE Conference on Decision and Control (CDC), 
%                   Dez. 2010, S. 6565–6570. doi: 10.1109/CDC.2010.5717269.
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
