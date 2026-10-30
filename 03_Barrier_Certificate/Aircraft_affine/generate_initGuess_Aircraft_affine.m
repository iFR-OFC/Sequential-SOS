function x0 = generate_initGuess_Aircraft_affine(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


A0 = full(subs( nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

% platform flag: compute the LQR gain on macOS and store it; on Windows
% and Linux load the stored values so that the initial guess is identical
% across platforms (lqr/LAPACK is platform dependent)
initfile = fullfile(fileparts(mfilename('fullpath')),'initGuess_Aircraft_affine.mat');
isMac    = strncmp(computer,'MAC',3);

if isMac
    [K,P] = lqr(A0,B0,eye(4),1);
    save(initfile,'K','P');
else
    assert(isfile(initfile),'Initial guess file "%s" not found. Generate it on macOS first.',initfile)
    S = load(initfile);
    K = S.K;
    P = S.P;
end

% initial guess
B0 = x'*P*x;
s20 = (x'*x)^2;
s30 = x'*x;
s40 = x'*x;
s50 = x'*x;
K0  = -K*x;
b0  = 1;


x0 = casos.PD([B0;s20;s30;s40;s50;K0;b0]);


end