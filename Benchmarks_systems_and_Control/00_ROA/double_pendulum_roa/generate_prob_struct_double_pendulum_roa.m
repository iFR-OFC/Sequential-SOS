% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for ROA problem using
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

function [probStruct] = generate_prob_struct_double_pendulum_roa(degV,degS)

x = casos.Indeterminates('x',4,1);
u = casos.Indeterminates('u',2,1);

% generate closed-loop dynamics
f = doublePendulumDynamics(x,u);

% generate a desired safe set as superellipsoid
n = 2;
g = (x(1)^2/1^2)^(n/2) + (x(2)^2/2^2)^(n/2) + (x(3)^2/1^2)^(n/2) + ...
     (x(4)^2/2^2)^(n/2) - 1;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,degV));

% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,degS));

% fill structure
probStruct.x    = x;        % indeterminates
probStruct.f    = f;        % dynamics
probStruct.g    = g;        % safe set
probStruct.s1   = s1;       % multiplier
probStruct.V    = V;        % Lyapunov function
probStruct.b    = 1;        % constant stable level set
probStruct.Name = 'Double-pendulum'; % Name of example

end

function [xdot,P] = doublePendulumDynamics(x,u)



f2 = -3.447*x(1)^3 + 2.350*x(1)^2*x(3) + 1.303*x(1)*x(2)^2 ...
    + 3.939*x(3)^3 + 21.520*x(1) - 5.000*x(3);
f4 =  4.023*x(1)^3 - 36.551*x(1)^2*x(3) - 4.131*x(2)^2*x(3) ...
    - 27.060*x(3)^3 - 25.115*x(1) + 77.700*x(3);


B = [  0     0  ;
    8   -31.2;
    0     0  ;
    -31.2 391.2];

xdot = [x(2); f2; x(4); f4] + B*u;

% linearize system

A0 = full(subs(nabla(xdot,x),[x;u],zeros(6,1)));
B0 = full(subs(nabla(xdot,u),[x;u],zeros(6,1)));

Q = eye(4);
R = eye(2);

[K,P] = lqr(A0,B0,Q,R);

xdot = subs(xdot,u,-K*x);

end
