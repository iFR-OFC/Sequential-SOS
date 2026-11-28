% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for ROA problem using
%                       a F-A/18 dynamical system from [1].
%
%
%   Reference: [1] Susceptibility of F/A-18 Flight Controllers to the 
%                  Falling-Leaf Mode: Nonlinear Analysis Abhijit 
%                  Chakraborty, Peter Seiler, and Gary J. Balas Journal of 
%                  Guidance, Control, and Dynamics 2011 34:1, 73-85
%
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_f18_roa(probStruct_F18)

f = probStruct_F18.f;
x = probStruct_F18.x;
A_cl = subs(nabla(f,x),x,zeros(7,1));

P = lyap(full(A_cl)',eye(7));

Vval = x'*P*x;

% initial guess
V0 = Vval;
s10 = (x'*x);

x0 = casos.PD([V0;s10]);

end