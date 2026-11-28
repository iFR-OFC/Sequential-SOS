% *****+*******************************************************************
%
%   Short Description:  Implements the initial guess for control design
%                       problem using a Jet-Engine dynamical system from
%                       [1].
%
%
%   Reference: [1]	I. R. Manchester und J.-J. E. Slotine, „Output-Feedback
%                   Control of Nonlinear Systems using Control Contraction
%                   Metrics and Convex Optimization", 5. Juni 2014, arXiv:
%                   arXiv:1406.1256. doi: 10.48550/arXiv.1406.1256.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function x0 = generate_init_guess_jet_engine_ctrl_design(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


A0 = full(subs( nabla(f,x),[x;u],zeros(3,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(3,1)));

[K,P] = lqr(A0,B0,eye(2),1);

% initial guess
V0 = x'*P*x;
s10 = (x'*x)^2;
s20 = x'*x;
s30 = x'*x;
s40 = x'*x;
kappa0  = -K*x;
b0  = 0.1;


x0 = casos.PD([V0;s10;s20;s30;s40;kappa0;b0]);


end