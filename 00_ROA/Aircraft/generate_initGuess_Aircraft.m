function x0 = generate_initGuess_Aircraft(probStruct)

f = probStruct.f;
x = probStruct.x;

A_cl = subs(nabla(f,x),x,zeros(4,1));


P = lyap(full(A_cl)',eye(4)*0.1);

Vval = x'*P*x;

% initial guess
V0 = Vval;
% s10 = (x'*x);

s1_sparsity = sparsity(probStruct.s1);
s10 = casos.PD(s1_sparsity, ones(s1_sparsity.nnz,1));

x0 = casos.PD([V0;s10]);

end