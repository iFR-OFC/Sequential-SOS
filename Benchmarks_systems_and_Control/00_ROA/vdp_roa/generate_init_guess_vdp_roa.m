% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for ROA problem using
%                       a VDP dynamical system from [1].
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

function x0 = generate_init_guess_vdp_roa(probStruct)

f = probStruct.f;
x = probStruct.x;

A_cl = subs(nabla(f,x),x,zeros(length(x),1));

P = lyap(full(A_cl)',eye(length(x)));

Vval = x'*P*x;

% initial guess
V0 = Vval;
s10 = (x'*x);
s1_sparsity = sparsity(probStruct.s1);
s10 = casos.PD(s1_sparsity, ones(s1_sparsity.nnz,1));

x0 = casos.PD([V0;s10]);

end
