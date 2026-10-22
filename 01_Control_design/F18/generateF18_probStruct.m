%% ------------------------------------------------------------------------
%
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       region-of-attraction for the longitudinal motion
%                       of the Nasa Generic Transport Model. To increase
%                       the size of the sublevel set we try to minimize the
%                       squared distance to a defined set. Additionally, we
%                       synthesis a linear control law at the same time.
%                       State and control constraints are considered.
%
%   Reference: Susceptibility of F/A-18 Flight Controllers to the 
%   Falling-Leaf Mode: Nonlinear Analysis Abhijit Chakraborty, 
%   Peter Seiler, and Gary J. Balas Journal of Guidance, Control, and 
%   Dynamics 2011 34:1, 73-85
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%--------------------------------------------------------------------------

function [probStruct] = generateF18_probStruct(degV,degS)

x = casos.Indeterminates('x',7,1);

f = f18Dynamics(x);

d2r = pi/180;

% generate a desired safe set
Dmax = diag([10*d2r 25*d2r 35*d2r 30*d2r 15*d2r 25*d2r 20*d2r]);  
N = inv(Dmax^2);   %  Scale by inverse of max state values
g = x'*N*0.1*x-1;

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
probStruct.Name = 'F/A-18'; % Name of example

end

 