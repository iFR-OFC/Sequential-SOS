function x0 = generate_initGuess_Pendubot(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


A0 = full(subs( nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

[K,P] = lqr(A0,B0,eye(4),1);

% initial guess
B0 = probStruct.g;
s10 = (x'*x);
s20 = x'*x;
s30 = x'*x;
s40 = x'*x;
kappa0  = -K*x;
b0  = 1;


x0 = casos.PD([B0;s10;s20;s30;s40;kappa0;b0]);


end