function x0 = generate_initGuess_MiMo(probStruct)

f = probStruct.f;
x = probStruct.x;
u = probStruct.u;


A0 = full(subs(nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

[K,P] = lqr(A0,B0,eye(3),eye(2));


% initial guess
B0  = x'*P*x;
s10 = (x'*x);
s20 = x'*x*ones(length(u),1);
s30 = x'*x*ones(length(u),1);
s40 = x'*x;
kappa0  = -K*x;
b0  = 1;


x0 = casos.PD([B0;s10;s20;s30;s40;kappa0;b0]);


end