
function [probStruct] = generateAircraft_benchmark(degV,degS)

x = casos.Indeterminates('x',4);

% Polynomial Dynamics
f = GTM_dynamics(x(1),x(2),x(3),x(4));


G0 = [3.95382671347060	-0.230032507976836	 0.0316965275615692	 0.292992065909380
    -0.230032507976836	 0.904764915483640	-0.161191428789579	-1.32594376986430
    0.0316965275615692	-0.161191428789579	 0.0344002648214491	 0.242292666551058
    0.292992065909380	-1.32594376986430	 0.242292666551058	 2.02114797577027];

g = (x'*G0*x) - 1;


% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,degV));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,degS));


% fill structure
probStruct.x    = x;        % indeterminates
probStruct.f    = f;        % dynamics
probStruct.g    = g;        % safe set
probStruct.s1   = s1;       % multiplier
probStruct.V    = V;        % Lyapunov function
probStruct.b    = 1;        % constant stable level set
probStruct.Name = 'Aircraft'; % Name of example

end

 