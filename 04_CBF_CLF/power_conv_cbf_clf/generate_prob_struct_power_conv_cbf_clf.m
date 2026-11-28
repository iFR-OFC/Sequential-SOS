% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for CBF-CLF problem
%                       using a Power Converter dynamical system from [1].
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

function [probStruct] = generate_prob_struct_power_conv_cbf_clf(degV,degH,degK)


x = casos.Indeterminates('x',3,1);
u = casos.Indeterminates('u',2,1);


f = power_converter_dynamics(x,u);
umax =  [1;1];
umin = -umax;

% generate a desired safe set
g = (x(1)/0.8)^2 + (x(2)/1.2)^2 + (x(3)/1.2)^2 - 1;

%% SOS decision variables

% CLF (we used adjusted monomial here i.e, we removed very small terms
% found by the initial guess)
V  = casos.PS.sym('v',monomials(x,degV));

% CBF
h  = casos.PS.sym('h',monomials(x,degH));

% control law (common polynomial vector)
kappa = casos.PS.sym('kappa',monomials(x,degK),[length(u) 1]);

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,2));
s2 = casos.PS.sym('s2',monomials(x,2));
s3 = casos.PS.sym('s3',monomials(x,2),[length(u) 1]);
s4 = casos.PS.sym('s4',monomials(x,2),[length(u) 1]);
s5 = casos.PS.sym('s5',monomials(x,2));

% rates of the comparison functions
gammaV = 0.0001; %casos.PS.sym('gammaV');
gammaW = 0.0001; %casos.PS.sym('gammaW');

% level set of terminal set
b = casos.PS.sym('b');


% fill structure
probStruct.x     = x;        % indeterminates
probStruct.u     = u;
probStruct.f     = f;        % dynamics
probStruct.g     = g;        % safe set

probStruct.V     = V;        % CLF
probStruct.h     = h;        % CBF
probStruct.kappa = kappa;    % control law

probStruct.s1    = s1;       % SOS multiplier
probStruct.s2    = s2;
probStruct.s3    = s3;
probStruct.s4    = s4;
probStruct.s5    = s5;

probStruct.gammaV = gammaV;
probStruct.gammaW = gammaW;

probStruct.b     = b;

probStruct.umin  = umin;
probStruct.umax  = umax;

probStruct.Name = 'Power-Covnverter'; % Name of example

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