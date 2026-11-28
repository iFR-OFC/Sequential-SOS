% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for ROA problem using
%                       a Double-pendulum dynamical system from [1].
%
%
%   Reference: [1]	S. W. Smith, H. Yin, und M. Arcak, „Continuous
%                   Abstraction of Nonlinear Systems using Sum-of-Squares
%                   Programming", in 2019 IEEE 58th Conference on Decision
%                   and Control (CDC), Dez. 2019, S. 8093–8098. doi:
%                   10.1109/CDC40024.2019.9029776.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_double_pendulum_roa(probStruct)

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
