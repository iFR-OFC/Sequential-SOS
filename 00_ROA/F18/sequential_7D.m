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
set(groot, ...
    'defaultTextInterpreter', 'latex', ...
    'defaultAxesTickLabelInterpreter', 'latex', ...
    'defaultLegendInterpreter', 'latex');

close all
clear
clc


%% system states
x = casos.Indeterminates('x',7,1);

import casos.toolboxes.sosopt.cleanpoly


f1  = @(beta,alpha,p,q,r,phi,xcR)...
  3.153345014754952e-06*alpha^3 - 0.2065171402380763*alpha^2*beta ...
  - 0.3010445810254201*alpha^2*p + 0.001958571170322954*alpha^2*phi ...
  + 0.3558189687845597*alpha^2*r + 0.02576468558743544*alpha^2*xcR ...
  - 0.001360100752311846*alpha*beta^2 + 0.05556920090648195*alpha*beta...
  *phi - 0.3634044556134873*alpha*beta*q - 0.0004813582585402263*alpha...
  *phi^2 + 0.0577180044959827*beta^3 + 0.02033521373538973*beta^2*p ...
  - 0.01963707488336307*beta^2*phi + 0.03254761688773859*beta^2*r ...
  - 0.02700350462763341*beta^2*xcR + 0.00356302754006645*beta*phi^2 ... 
  - 0.006644423779964928*phi^3 + 0.02403995325667869*alpha^2 ...
  + 0.1716908664450242*alpha*beta + 0.9247080521758517*alpha*p ...
  - 0.006328286888283467*alpha*phi + 0.4181243993346678*alpha*r ...
  - 0.05287123183982542*alpha*xcR - 0.0245428159495927*beta^2 ...
  + 0.02606202420536543*beta*phi + 0.2707928746853959*beta*q ...
  - 0.02771270046902687*phi^2 - 0.001433627264618562*alpha ...
  + 0.02750398858742245*beta + 0.3529242183394632*p + 0.07252762792843861...
  *phi - 0.9372303643204858*r + 0.01429020844937823*xcR;
 
f2 = @(beta,alpha,p,q,r,phi,xcR)...
  -0.2139842122459398*alpha^3 + 0.00755027103112707*alpha^2*beta ...
  + 0.03540521522954149*alpha^2*phi - 1.132109984228988*alpha^2*q ...
  - 0.01845498597078899*alpha*beta^2 + 0.6357034892076658*alpha*beta*p ...
  - 4.181742495065986e-05*alpha*beta*phi - 0.7499070186028067*alpha*beta...
  *r + 0.01029228436846721*alpha*phi^2 - 0.004364563088207579*beta^3 ...
  - 0.004154373502286115*beta^2*phi + 0.1988252203389781*beta^2*q ...
  - 6.882451189841364e-05*beta*phi^2 + 0.0125161493105963*phi^3 ...
  + 0.3637036100312517*alpha^2 - 0.05181028530013231*alpha*beta ...
  + 0.01364058233723217*alpha*phi + 0.6940980395563138*alpha*q... 
  - 0.02242776933969258*beta^2 - 0.9576178724308215*beta*p ...
  + 0.00010932571006722*beta*phi - 0.3619062637698365*beta*r ... 
  - 0.03648758223043212*phi^2 - 0.2299147285988287*alpha ...
  + 0.001870494870104267*beta - 0.04688445933824111*phi ...
  + 0.7259694904194895*q;  
 
f3 = @(beta,alpha,p,q,r,phi,xcR)...
  -0.4414814899115931*alpha^3 - 23.21988319558105*alpha^2*beta ...
  - 5.103450502189039*alpha^2*p - 0.747599021664727*alpha^2*phi ...
  + 7.452697515311997*alpha^2*r - 2.226986383574201*alpha^2*xcR ...
  - 0.2556161723216137*alpha*beta^2 + 20.20077533202648*beta^3 ...
  + 0.2031166069052403*alpha^2 + 20.65261864937506*alpha*beta ...
  + 7.495897993314138*alpha*p + 1.149052159373335*alpha*phi ...
  - 16.51903205850301*alpha*r + 0.2822677706670466*alpha*xcR... 
  + 0.06666998057975507*beta^2 - 0.031731027233146*p*q ...
  - 0.8150553426917272*q*r + 0.06122624680466708*alpha - 9.701122668749656...
  *beta - 3.923306096425065*p - 0.6102495853689766*phi ...
  - 0.03688825706037443*q + 9.365411189648313*r + 0.5311015383180122*xcR;
 
f4  = @(beta,alpha,p,q,r,phi,xcR) ...
  1.553816147021179*alpha^3 + 17.13290975887215*alpha^2*q ...
  - 2.174664690578164*alpha^2 + 4.404153897583718*alpha*q ...
  + 0.01963739144013944*p^2 + 0.971261062970545*p*r - 0.01963739726930736...
  *r^2 - 2.303192595016487*alpha + 0.04393012243878198*p ...
  - 14.55941295884036*q - 0.0202593409751493*r;
 
f5  = @(beta,alpha,p,q,r,phi,xcR) ... 
  -0.2469297732870805*alpha^3 - 2.324509437528913*alpha^2*beta ...
  - 0.9357096905063863*alpha^2*p + 1.715633429424422*alpha^2*r ...
  - 0.5427154642650244*alpha^2*xcR + 0.09538116895297807*alpha*beta^2 ... 
  - 0.04018209136994876*beta^3 + 0.1780879838058609*alpha^2  ...
  - 1.419260863469877*alpha*beta + 0.5264411344974868*alpha*p ...
  - 0.8898805459780623*alpha*r + 0.4694093264593382*alpha*xcR ...
  - 0.02518659933031223*beta^2 - 0.7543542152936564*p*q ...
  + 0.0317310272331471*q*r - 0.01344081050039049*alpha ...
  + 0.5455270560107076*beta + 0.04253549876779805*p + 0.007623736588361647...
  *phi + 0.0157913770978545*q - 0.6026906228309178*r - 0.311390752923031...
  *xcR;
 
f6  = @(beta,alpha,p,q,r,phi,xcR) ...
  -0.1481374591474449*phi^2*q - 0.0722557699944649*phi^2*r ...
  + 0.2921438187540301*phi*q - 0.218166704226111*phi*r ...
  + 0.9999999999999888*p + 0.1940618826520113*q + 0.2771490908942477*r;
 
f7  = @(beta,alpha,p,q,r,phi,xcR) 4.9*r - 1*xcR;



% get trim point 
f = @(x) [
                    f1(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),x(7,:))
                    f2(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),x(7,:))
                    f3(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),x(7,:))
                    f4(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),x(7,:))
                    f5(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),x(7,:))
                    f6(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),x(7,:))
                    f7(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:),x(7,:))
];


 f = f(x*1);
% [x0, u0] = findtrim(f,x0, u0);
d2r = pi/180; r2d = 1/ d2r;


Dmax = diag([10*d2r 25*d2r 35*d2r 30*d2r 15*d2r 25*d2r 20*d2r]);  
N = inv(Dmax^2);   %  Scale by inverse of max state values
g = x'*N*0.1*x-1;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2:4));

% SOS multiplier
s1    = casos.PS.sym('s1',monomials(x,2));


% enforce positivity
l = 1e-6*(x'*x);


%% setup solver
% options
opts = struct('sossol','mosek');
opts.verbose  = 1;
opts.max_iter = 100;
opts.almostOptCount = 10;
% opts.tolerance_con = 1e-4;
% opts.conVioCheck = 'sampling';

opts.sossol_options.newton_solver = []; 


% b     = casos.PS.sym('b');
b = 1;

% cost
cost = dot(g-(V-b),g-(V-b)) ;

sos = struct('x',[V; s1],... % decision variables
    'f',cost, ...              % cost
    'p',[]);                   % parameter

% SOS constraints
sos.('g') = [s1;
    V-l;
    s1*(V-b)-nabla(V,x)*f-l;
    ];

% states + constraint cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));
% opts.sossol_options.newton_solver =[];

% setup solver
S = casos.nlsossol('S1','sequential',sos,opts);

%% solve problem


A_cl = subs(nabla(f,x),x,zeros(7,1));

P = lyap(full(A_cl)',eye(7));

Vval = x'*P*x;

% initial guess
V0 = Vval;
s20 = (x'*x);
% b0  = 1;
x0 = casos.PD([V0;s20]);

% solve problem
sol = S('x0' ,x0);

% casos.postProcessSolver(S,true);

S.stats

S.info.sossol_convex.gram.Kc

%%
cost      = nan(S.stats.iterations,1);
maxConVio = nan(S.stats.iterations,1);
alpha     = nan(S.stats.iterations,1);
counterSoC = nan(S.stats.iterations,1);
counterLineserchSteps = nan(S.stats.iterations,1);
counterFeasRes = nan(S.stats.iterations,1);
for k = 1:S.stats.iterations
    cost(k) = S.stats.single_iterations{k}.cost;
    maxConVio(k) = S.stats.single_iterations{k}.maxConVio;
    alpha(k) = S.stats.single_iterations{k}.stepLength;
    counterFeasRes(k) = S.stats.single_iterations{k}.counterFeasRes;
    counterSoC(k) = S.stats.single_iterations{k}.counterSoC;
    counterLineserchSteps(k) = S.stats.single_iterations{k}.counterLineserchSteps;
end


load CD_iter_7states.mat
load CD_iter_7states2.mat
maxConVio(maxConVio == 0) = eps;
% Plot convergence diagnostics
figure
tiledlayout(1,2)
nexttile, semilogy(abs(cost(2:end)-cost(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,
semilogy(abs(cost_CD(2:end)-cost_CD(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,
semilogy(abs(cost_CD2(2:end)-cost_CD2(1:end-1)),'LineWidth',1.2), grid on, title('Cost'), hold on,
      plot([0 length(cost_CD)],[S.stats.options.tolerance_opt S.stats.options.tolerance_opt],'k--','LineWidth',1.2) 
      xlabel('Iterations')
      ylabel('$ |f(\xi_{k+1}) - f(\xi_{k}) |$')
      legend('Sequential','Coordinate-descent','$\epsilon_\mathrm{opt}$')
nexttile, semilogy(maxConVio,'LineWidth',1.2), grid on, title('Maximum constraint violation') , hold on 
plot([0 S.stats.iterations],[S.stats.options.tolerance_con S.stats.options.tolerance_con],'k--','LineWidth',1.2)
plot([0 S.stats.iterations],[eps eps],'k-.','LineWidth',1.2)
xlabel('Iterations')
ylabel('$\min(\mathrm{eps},\theta(\xi_{k}))$')
legend('$ \theta(\xi_{k})$','$\epsilon_\mathrm{con}$','eps')
