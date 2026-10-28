function x0 = generate_initGuess_MiMo(probStruct)

f = probStruct.f;
x = probStruct.x;
u = probStruct.u;


A0 = full(subs(nabla(f,x),[x;u],zeros(5,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(5,1)));

% platform flag: compute the LQR gain on macOS and store it; on Windows
% and Linux load the stored values so that the initial guess is identical
% across platforms (lqr/LAPACK is platform dependent)
initfile = fullfile(fileparts(mfilename('fullpath')),'initGuess_MiMo.mat');
isMac    = strncmp(computer,'MAC',3);

if isMac
    [K,P] = lqr(A0,B0,eye(3),eye(2));
    save(initfile,'K','P');
else
    assert(isfile(initfile),'Initial guess file "%s" not found. Generate it on macOS first.',initfile)
    S = load(initfile);
    K = S.K;
    P = S.P;
end


% initial guess
V0  = x'*P*x;
s10 = (x'*x);
s20 = x'*x*ones(length(u),1);
s30 = x'*x*ones(length(u),1);
s40 = x'*x;
K0  = -K*x;
b0  = 1;


x0 = casos.PD([V0;s10;s20;s30;s40;K0;b0]);


end