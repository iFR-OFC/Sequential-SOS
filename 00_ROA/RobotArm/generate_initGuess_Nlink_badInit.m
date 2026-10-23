function x0 = generate_initGuess_Nlink_badInit(probStruct)

x = probStruct.x;

% clearly infeasible initial guess: negative unit polynomials for both the
% Lyapunov function and the SOS multiplier (cf. ROA_Nlink_sequential_badInit)
V0  = -(x'*x);
s10 = -(x'*x);

% initial guess
x0 = casos.PD([V0;s10]);

end
