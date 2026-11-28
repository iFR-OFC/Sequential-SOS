% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for barrier
%                       certificate problem using a Power Converter
%                       dynamical system.
%
%
%   Reference: [1]	M. Schneeberger, F. Dörfler, und S. Mastellone, „SOS
%                   Construction of Compatible Control Lyapunov and Barrier
%                   Functions", IFAC-PapersOnLine, Bd. 56, Nr. 2, S.
%                   10428–10434, Jan. 2023, doi:
%                   10.1016/j.ifacol.2023.10.1058.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function [probStruct] = generate_prob_struct_power_conv_barrier_cert(degV,degK)

x = casos.Indeterminates('x',3,1);
u = casos.Indeterminates('u',2,1);


f = power_converter_dynamics(x,u);
umax =  [1;1];
umin = -umax;

% generate a desired safe set
g = (x(1)/0.8)^2 + (x(2)/1.2)^2 + (x(3)/1.2)^2 - 1;

% Lyapunov function candidate
B = casos.PS.sym('v',monomials(x,degV));


% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,2));
s2    = casos.PS.sym('s2',monomials(x,2),[length(u) 1]);
s3    = casos.PS.sym('s3',monomials(x,2),[length(u) 1]);
s4    = casos.PS.sym('s4',monomials(x,2));
b     = casos.PS.sym('b');
kappa = casos.PS.sym('kappa',monomials(x,degK),[length(u) 1]);


% fill structure
probStruct.x    = x;        % indeterminates
probStruct.u    = u;
probStruct.f    = f;        % dynamics
probStruct.g    = g;        % safe set


probStruct.B    = B;        % Lyapunov function
probStruct.s1   = s1;       % multiplier
probStruct.s2   = s2;
probStruct.s3   = s3;
probStruct.s4   = s4;
probStruct.b    = b;
probStruct.umin = umin;
probStruct.umax = umax;

probStruct.kappa = kappa;

probStruct.Name = 'Power-Converter'; % Name of example

end


function f = power_converter_dynamics(x,u)

x1 = x(1);
x2 = x(2);
x3 = x(3);

f0 = [-0.05*x1 - 57.9*x2 + 0.00919*x3;
    1710*x1 + 314*x3;
    -0.271*x1 - 314*x2];

g0 = [0.05-57.9, -57.9*x3;
    1710+1710*x1, 0;
    0, 1710 + 1710*x1];


f = f0+g0*u*1;
end