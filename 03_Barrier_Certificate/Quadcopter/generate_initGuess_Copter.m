function x0 = generate_initGuess_Copter(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


% A0 = full(subs( nabla(f,x),[x;u],zeros(8,1)));
% B0 = full(subs(nabla(f,u),[x;u],zeros(8,1)));

% initial guess
B0 = probStruct.g;
s10 = (x'*x);
s20 = x'*x*ones(length(u),1);
s30 = x'*x*ones(length(u),1);
s40 = x'*x;
kappa0  = ones(2,1);
b0  = 1;


x0 = casos.PD([B0;s10;s20;s30;s40;kappa0;b0]);


end