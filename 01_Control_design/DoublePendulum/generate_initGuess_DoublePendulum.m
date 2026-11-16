function x0 = generate_initGuess_DoublePendulum(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


A0 = full(subs(nabla(f,x),[x;u],zeros(6,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(6,1)));

Q = eye(4);
R = eye(2);

[K,P] = lqr(A0,B0,Q,R);

% initial guess
V0 = x'*P*x;
s10 = (x'*x);
s20 = ones(length(u),1)*x'*x;
s30 = ones(length(u),1)*x'*x;
s40 = x'*x;
kappa0  = -K*x;
b0  = 1;


x0 = casos.PD([V0;s10;s20;s30;s40;kappa0;b0]);


end