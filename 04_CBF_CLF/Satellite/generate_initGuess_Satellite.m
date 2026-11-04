function x0 = generate_initGuess_Satellite(probStruct)


x     = probStruct.x;
u     = probStruct.u;
f_raw = probStruct.f_raw;
Dx    = probStruct.Dx;
g     = probStruct.g;


% trim point
x0trim = [0 0 0 0 0 0]';
u0     = [0,0,0]';

A0 = full(casos.PD(subs(nabla(f_raw,x),[x;u],[x0trim;u0])));
B0 = full(casos.PD(subs(nabla(f_raw,u),[x;u],[x0trim;u0])));

% cost function weights
Q = diag([1, 1, 1, 1, 1 ,1]);
R = eye(3)*1;

% generate an initial guess
[K,P] = lqr(full(A0),full(B0),Q,R);


% named initial guess for each decision variable
% (ordering: V,h,s1,s2,s3,s4,s5,kappa,b,gammaW,gammaV)

V0 = (inv(Dx)*x)'*P*(inv(Dx)*x);  % scaled terminal penalty (Lyapunov of the linearized system)
h0 = casos.PD(g);                 % start at the boundary of the safe set

s10 = x'*x;
s20 = x'*x;
s30 = ones(length(u),1);
s40 = ones(length(u),1);
s50 = x'*x;

kappa0 = -K*x;                    % LQR control law

b0 = 0.8;                         % level set of the terminal set

gammaW0 = 0.01;
gammaV0 = 0.01;


% initial guess for sequential / coordinate descent
x0 = casos.PD([V0; h0; s10; s20; s30; s40; s50; kappa0; gammaW0; gammaV0;b0]);

end
