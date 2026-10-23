%--------------------------------------------------------------------------
%
%   Short Description:  Generate the problem structure (dynamics, safe
%                       set, polynomial decision variables) for the ROA
%                       estimation of the N-link robot arm with 2*N states.
%
%   Reference:  T. Cunis and J. Olucak,
%   "CaΣoS: A nonlinear sum-of-squares optimization suite⋆,"
%    2025 American Control Conference (ACC), Denver, CO,
%    USA, 2025, pp. 1659-1666, doi: 10.23919/ACC63710.2025.11107794.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%--------------------------------------------------------------------------

function [probStruct] = generateNlink_probStruct(n,deg)

% degree of the Taylor expansion of the dynamics and of the polynomial
% decision variables (cf. robot/nLinkRobot_dyn_poly_n*_d*.m)
if nargin < 2
    deg = 2;
end

x = casos.PS('x',2*n,1);

% folder of the robot-arm model (dynamics and LQR data)
robotDir = fullfile(fileparts(mfilename('fullpath')),'robot');
addpath(robotDir)

% get dynamics (Taylor expansion of degree deg)
f = feval(['nLinkRobot_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);

% load P matrix of the LQR-stabilized system
load(fullfile(robotDir,['data_n' num2str(n) '.mat']),'P');

% build a safe set
g = x'*P*0.01*x-2;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,deg));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,deg));

% fill structure
probStruct.x    = x;        % indeterminates
probStruct.f    = f;       % dynamics
probStruct.g    = g;       % safe set
probStruct.s1   = s1;      % multiplier
probStruct.V    = V;       % Lyapunov function
probStruct.b    = 1;       % fixed level-set constant (gamma) of the V-s-iteration
probStruct.P    = P;       % LQR matrix (used for the initial guess)
probStruct.Name = ['N-link robot arm (N=' num2str(n) ')']; % Name of example

end
