%% ------------------------------------------------------------------------
%   
%   Short Descirption:  Calculate a compatible pair of a CBF-CLF
%                       for a satellite using sequential SOS. This example
%                       is taken from [1] and [2].
%
%   Reference:  
%   [1] Olucak, Jan and Castello Branco de Oliveira, Arthur and 
%   Cunis, Torbjørn - Safe-by-Design: Approximate Nonlinear Model 
%   Predictive Control with Real Time Feasibility, submitted to IEEE
%   Transaction on Automatic Control, pre-print: 
%   https://arxiv.org/abs/2509.22422
%
%   [2] Olucak, Jan and Castello Branco de Oliveira, Arthur and 
%   Cunis, Torbjørn - Supplementary Material for: Safe-by-Design 
%   Approximate Nonlinear Model Predictive Control with Realtime 
%   Feasibility, doi = {10.18419/DARUS-5297}
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
% ------------------------------------------------------------------------

clear
close all
% clc

% system states
x = casos.PS('x',6,1);
u = casos.PS('u',3,1);

%% Hubble telescope parameter
J = diag([31046;77217;78754]);

% simple bounds on rates;
omegaMax1 = 0.5*pi/180;
omegaMax2 = 0.2*pi/180;
omegaMax3 = 0.2*pi/180;

x_low =  [-omegaMax1 -omegaMax2 -omegaMax3]';
x_up  =  [ omegaMax1  omegaMax2  omegaMax3]';


% control constraint; assumption is that the box is inside the full
% torque volume. This is roughly estimated visually.
umin = [-1 -1 -1]'*1.2;
umax = [ 1  1  1]'*1.2;

Dx   = diag([1/(x_up(1)-x_low(1)),1/(x_up(2)-x_low(2)),1/(x_up(3)-x_low(3)),1,1,1]); 

Dxin = inv(Dx);

%% dynamics
% cross-product matrix
cpm = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];

% dynamics
B = @(sigma) (1-sigma'*sigma)*eye(3)+ 2*cpm(sigma)+ 2*sigma*sigma';

f =  [-J\cpm(x(1:3))*J*x(1:3) + J\u;
      1/4*B(x(4:6))*x(1:3)]; % omega_dot

% trim point
x0    = [0 0 0 0 0 0]';
u0    = [0,0,0]';

A0 = full(casos.PD(subs(nabla(f,x),[x;u],[x0;u0])));
B0 = full(casos.PD(subs(nabla(f,u),[x;u],[x0;u0])));

% cost function weights
Q = diag([1, 1, 1, 1, 1 ,1]);
R = eye(3)*1;

% generate an initial guess 
[K0,P0] = lqr(full(A0),full(B0),Q,R);

% scaled initial guess for terminal penalty (Lyapunov linear system)
Wval = (inv(Dx)*x)'*P0*(inv(Dx)*x);

% scale dynamics
f = Dx*subs(f,[x;u],[Dx\x;u]);

% state constraint
n = 2;
g0 = (x(1)^2/omegaMax1^2)^(n/2) + (x(2)^2/omegaMax2^2)^(n/2) + (x(3)^2/omegaMax3^2)^(n/2) + ...
     (x(4)^2/0.57^2)^(n/2) + (x(5)^2/0.57^2)^(n/2) + (x(6)^2/0.57^2)^(n/2) - 1;

% re-scale input of state constraints
g = subs(g0,x,Dx\x); 

%% setup SOS problem


% CBF
W  = casos.PS.sym('w',monomials(x,2));

% CLF
V  = casos.PS.sym('v',monomials(x,2:4)); % we used adjusted monomial here i.e, we removed very small terms found by the intial guess

Kd = casos.PS.sym('kd',[3,1]);
Kp = casos.PS.sym('kp',[3,1]);

K = -Kd.*x(1:3)-Kp.*x(4:6);

% SOS mulitplier
s1 = casos.PS.sym('s1',monomials(x,2));
s2 = casos.PS.sym('s2',monomials(x,2));
s3 = casos.PS.sym('s3',monomials(x,0),[3 1]);
s4 = casos.PS.sym('s4',monomials(x,0),[3 1]);
s5 = casos.PS.sym('s5',monomials(x,2));


% fixed level set of terminal set
b = 0.9;

% options for sequential sos
opts = struct('sossol','mosek');
opts.con_vio_check = 'signed-distance'; 'projection';
opts.verbose       = 1;
opts.max_iter      = 100;
% opts.tolerance_con = 1e-1;
% opts.tolerance_opt = 1e-1;

% fixed constant for comparison functions
gammaV = 0.0025;
gammaW = 0.8;

cost = dot(g-(W-b),g-(W-b));

sos = struct('x', [W;V;Kd;Kp;s1;s2;s3; s4;s5],...  % decision variables
              'f', cost ,...                       % cost function
              'p',[]);                             % parameter

% constraints
sos.('g') = [
             s1
             s2;
             s3;
             s4;
             s5;
             V - 1e-6*(x'*x);                                      % CLF positivity              
             s1*(W-b) - g;                                         % State constraints
             s2*(W-b)  -  nabla(W,x)*subs(f,u,K) - gammaW*(W-b);   % CBF dissipation
             s3*(W-b)  + K-umin;                                   % control constraints
             s4*(W-b)  + umax-K;
             s5*(W-b) -  nabla(V,x)*subs(f,u,K) - gammaV*V;        % CLF dissipation
             ];

% states + constraint are linear/SOS cones
opts.Kx = struct('lin', length(sos.x));
opts.Kc = struct('sos', length(sos.g));
opts.sossol_options.newton_solver = []; 

% solver setup
S  = casos.nlsossol('S','sequential',sos,opts);

% initial guess for sequential
x0 = casos.PD([ g; ...
                Wval;
                diag(K0(1:3,1:3)); ...         
                diag(K0(1:3,4:6)); ...
                ones(3,1);
                ones(3,1);
                x'*x; ...
                x'*x; ...
                x'*x]);
 

%  solve
sol = S('x0',x0);

bsol = b;

% re-scale invariant set, terminal penalty and local control law
Wsol_re = subs(sol.x(1),x,Dx*x) - full(casos.PD(bsol)); % CBF as sublevel set
Vsol_re = subs(sol.x(2),x,Dx*x);

K       = -sol.x(3:5).*x(1:3)-sol.x(6:8).*x(4:6);
Ksol_re = subs(K,x,Dx*x);

%% Run a sos feasibility with the solution
% check if  solution is really SOS
con_fun       = casos.Function('f',{sos.x},{sos.g});

solg = remove_coeffs(con_fun(sol.x),1e-8);

for k = 1:length(solg)
sos_conv = struct('g',solg(k));

% constraint is scalar SOS cone
opts_conv = struct('Kc',struct('sos',length(solg(k))));
opts_conv.newton_solver =[];

% turn off error on fail, otherwise it throws an error if infeasible
% we want to manually check it below
opts_conv.error_on_fail = false;

% solve by relaxation to SDP
S_conv = casos.sossol('S','mosek',sos_conv,opts_conv);

% evaluate
sol_conv = S_conv();

if strcmp(S_conv.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS')
    disp('Given polynomial is SOS!')
else
    disp('No decomposition found!')
    k
end

end


% %% Check proj
% 
% gsol0 = remove_coeffs(gsol10,1e-6);
% 
% % Gram decision variable
% for k = 1:length(gsol0)
% 
% s = casos.PS.sym('q',grambasis(gsol0(k)));
% 
% 
%     e = s - gsol0(k);
% 
% 
% % projection error
% 
% 
% % define Q-SOS problem:
% %   min ||s-p||^2 s.t. s is SOS
% sos = struct('x',s,'f',dot(e,e));
% % states is scalar SOS cone
% opts = struct('Kx',struct('sos',length(s)));
% opts.error_on_fail = true;
% % solve by relaxation to SDP
% S = casos.sossol('S','mosek',sos,opts);
% % evaluate
% sol = S();
% 
% fprintf('Distance to SOS cone is %g.\n', full(sol.f))
% 
% end