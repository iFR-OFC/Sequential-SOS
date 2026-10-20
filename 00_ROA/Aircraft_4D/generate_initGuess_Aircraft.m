function x0 = generate_initGuess_Aircraft(probStruct)

f = probStruct.f;
x = probStruct.x;

A_cl = subs(nabla(f,x),x,zeros(4,1));

P = lyap(full(A_cl)',eye(4));

Vval = x'*P*x;

% initial guess
V0 = Vval;
s20 = (x'*x)^2;

x0 = casos.PD([V0;s20]);

end