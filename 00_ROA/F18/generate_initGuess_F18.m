function x0 = generate_initGuess_F18(probStruct_F18)

f = probStruct_F18.f;
x = probStruct_F18.x;
A_cl = subs(nabla(f,x),x,zeros(7,1));

P = lyap(full(A_cl)',eye(7));

Vval = x'*P*x;

% initial guess
V0 = Vval;
s20 = (x'*x);

x0 = casos.PD([V0;s20]);

end