function x0 = generate_initGuess_Aircraft(probStruct)

f = probStruct.f;
x = probStruct.x;

A_cl = subs(nabla(f,x),x,zeros(4,1));

% platform flag: compute the Lyapunov solution on macOS and store it; on
% Windows and Linux load the stored values so that the initial guess is
% identical across platforms (lyap/LAPACK is platform dependent)
initfile = fullfile(fileparts(mfilename('fullpath')),'initGuess_Aircraft.mat');
isMac    = strncmp(computer,'MAC',3);

if isMac
    P = lyap(full(A_cl)',eye(4));
    save(initfile,'P');
else
    assert(isfile(initfile),'Initial guess file "%s" not found. Generate it on macOS first.',initfile)
    S = load(initfile);
    P = S.P;
end

Vval = x'*P*x;

% initial guess
V0 = Vval;
s20 = (x'*x)^2;

x0 = casos.PD([V0;s20]);

end