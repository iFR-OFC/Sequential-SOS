function x0 = generate_initGuess_VDP(probStruct)


x = probStruct.x;
u = probStruct.u;
t = probStruct.t;
f = probStruct.f;
g = probStruct.g;

A0 = full(subs( nabla(f,x),[x;u],zeros(3,1)));
B0 = full(subs(nabla(f,u),[x;u],zeros(3,1)));

% platform flag: compute the LQR gain on macOS and store it; on Windows
% and Linux load the stored values so that the initial guess is identical
% across platforms (lqr/LAPACK is platform dependent)
initfile = fullfile(fileparts(mfilename('fullpath')),'initGuess_VDP.mat');
isMac    = strncmp(computer,'MAC',3);

if isMac
    [K,P] = lqr(A0,B0,eye(2),1);
    save(initfile,'K','P');
else
    assert(isfile(initfile),'Initial guess file "%s" not found. Generate it on macOS first.',initfile)
    S = load(initfile);
    K = S.K;
    P = S.P;
end




V0      =    x'*P*x; 
kappa0  =    -K*x; 
s10     =    (x'*x) + t^2; 
s20     =    x'*x + t^2; 
s30     =    (x'*x) + t^2;
s40     =    x'*x + t^2;
s50     =    (x'*x) + t^2;
s60     =    x'*x + t^2;
s70     =    (x'*x) + t^2;
s80     =    x'*x + t^2;
s90     =    (x'*x);
b0      =    0.8;

% initial guess
x0 = casos.PD([ V0;     ...
                kappa0 ; ...
                s10    ;...
                s20    ;...
                s30    ;...
                s40    ;...
                s50    ;...
                s60    ;...
                s70    ;...
                s80    ;...
                s90    ;...
                b0]);


end