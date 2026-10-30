function x0 = generate_initGuess_Satellite(probStruct)


x = probStruct.x;
u = probStruct.u;
f = probStruct.f;


% trim point
x0    = [0 0 0 0 0 0]';
u0    = [0,0,0]';

A0 = full(casos.PD(subs(nabla(f,x),[x;u],[x0;u0])));
B0 = full(casos.PD(subs(nabla(f,u),[x;u],[x0;u0])));

% cost function weights
Q = diag([1, 1, 1, 1, 1 ,1]);
R = eye(3)*1;

% platform flag: compute the LQR gain on macOS and store it; on Windows
% and Linux load the stored values so that the initial guess is identical
% across platforms (lqr/LAPACK is platform dependent)
initfile = fullfile(fileparts(mfilename('fullpath')),'initGuess_Satellite.mat');
isMac    = strncmp(computer,'MAC',3);

if isMac
    [K,P] = lqr(A0,B0,Q,R);
    save(initfile,'K','P');
else
    assert(isfile(initfile),'Initial guess file "%s" not found. Generate it on macOS first.',initfile)
    S = load(initfile);
    K = S.K;
    P = S.P;
end

% initial guess
B0 = casos.PD(probStruct.g);
s10 = (x'*x);
s20 = ones(length(u),1);
s30 = ones(length(u),1);
s40 = x'*x;
K0  = -K*x;
b0  = 1;


x0 = casos.PD([B0;s10;s20;s30;s40;K0;b0]);

end