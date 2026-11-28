% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for ROA problem using
%                       an aircraft dynamical system from [1].
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
