% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for ROA problem using
%                       a VDP-coupled dynamical system from [1].
%
%
%   Reference: [1]	X. Chen, „Reachability Analysis of Non-Linear Hybrid
%                   Systems Using Taylor Models", RWTH Aachen, Aachen,
%                   2015. [Online]. URL:
%                   https://publications.rwth-aachen.de/record/465295/files
%                   /465295.pdf
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_vdp_coupled_roa(probStruct)

f = probStruct.f;
x = probStruct.x;

A_cl = subs(nabla(f,x),x,zeros(length(x),1));

P = lyap(full(A_cl)',eye(length(x)));

Vval = x'*P*x;

% initial guess
V0 = Vval;
% s10 = (x'*x);
s1_sparsity = sparsity(probStruct.s1);
s10 = casos.PD(s1_sparsity, ones(s1_sparsity.nnz,1));

x0 = casos.PD([V0;s10]);

end
