% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for the
%                       Satellite-Attitude dynamical system.
%
%
%   Reference:
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_spacecraft_att_roa(probStruct)

x = probStruct.x;

% we want to start as close as possible at the constraint set
Vval = probStruct.g;

% initial guess
V0 = Vval;
s10 = (x'*x);
s1_sparsity = sparsity(probStruct.s1);
s10 = casos.PD(s1_sparsity, ones(s1_sparsity.nnz,1));

x0 = casos.PD([V0;s10]);

end
