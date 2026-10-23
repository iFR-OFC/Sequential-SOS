function x0 = generate_initGuess_Nlink(probStruct)

x = probStruct.x;
P = probStruct.P;

% Lyapunov function candidate from the LQR solution
V0 = x'*P*x;

% unit polynomial as initial multiplier
s10 = (x'*x);

% initial guess
x0 = casos.PD([V0;s10]);

end
