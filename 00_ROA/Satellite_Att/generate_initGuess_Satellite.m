function x0 = generate_initGuess_Satellite(probStruct)

x = probStruct.x;

% we want to start as close as possible at the constraint set
Vval = probStruct.g;

% initial guess
V0 = Vval;
s10 = (x'*x);

x0 = casos.PD([V0;s10]);

end