function x0 = generate_initGuess_Aircraft(probStruct)

f = probStruct.f;
x = probStruct.x;
V = probStruct.V;
g = probStruct.g;

s1 = probStruct.s1;
A_cl = subs(nabla(f,x),x,zeros(4,1));


P = lyap(full(A_cl)',eye(4)*0.1);

Vval = x'*P*x;

%% initial guess

x0 = {};

% 1.) simple, but feasible
V00 = Vval;
s100 = (x'*x);

x0{end+1} = casos.PD([V00;s100]);

% 2.) simple, but feasible
V00 = x'*x;
s100 = (x'*x);

x0{end+1} = casos.PD([V00;s100]);


% 3.) squared feasible
x0{end+1} = casos.PD([V00^2;s100^2]);


% 4.) simple, i.e., unit polynomial of used polynomials
V01 = casos.PD(V.sparsity,ones(V.nnz,1));
s101 = casos.PD(s1.sparsity,ones(s1.nnz,1));

x0{end+1} = casos.PD([V01;s101]);

% 5.) random polynomial (same sparsity, random coefficients (uniform))
V02 = casos.PD(V.sparsity,rand(V.nnz,1));
s102 = casos.PD(s1.sparsity,rand(s1.nnz,1));

x0{end+1} = casos.PD([V02;s102]);

% 6.) random polynomial (same sparsity, random coefficients (normal))
V03 = casos.PD(V.sparsity,randn(V.nnz,1));
s103 = casos.PD(s1.sparsity,randn(s1.nnz,1));

x0{end+1} = casos.PD([V03;s103]);

% 7.) just positve scalar values
x0{end+1} = casos.PD([1;1]);

% 8.) take problem information (here: minimum dist. to safe set)
x0{end+1} = casos.PD([g;s100]);


% 9.) infeasibe
x0{end+1} = -casos.PD([V00;s100]);




end