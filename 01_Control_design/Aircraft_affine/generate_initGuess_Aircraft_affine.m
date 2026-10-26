function x0 = generate_initGuess_Aircraft_affine(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


A0 = full(subs( nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

[K,P] = lqr(A0,B0,eye(4),1);

% initial guess
V0 = x'*P*x;
s20 = (x'*x)^2;
s30 = x'*x;
s40 = x'*x;
s50 = x'*x;
K0  = -K*x;
b0  = 1;


x0 = casos.PD([V0;s20;s30;s40;s50;K0;b0]);


end