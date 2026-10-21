function x0 = generate_initGuess_12State(probStruct)

f = probStruct.f;
x = probStruct.x;

A_cl = subs(nabla(f,x),x,zeros(12,1));

P = lyap(full(A_cl)',eye(12));

Vval = x'*P*x;

s20 = x'*x;

% initial guess
V0 = Vval;
x0 = casos.PD([V0;s20]);

end