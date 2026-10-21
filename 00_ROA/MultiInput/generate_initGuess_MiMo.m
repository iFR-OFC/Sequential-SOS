function x0 = generate_initGuess_MiMo(probStruct)

f = probStruct.f;
x = probStruct.x;

A_cl = subs(nabla(f,x),x,zeros(length(x),1));

P = lyap(full(A_cl)',eye(length(x)));

Vval = x'*P*x;

% initial guess
V0 = Vval;
s10 = (x'*x)^2;

x0 = casos.PD([V0;s10]);

end