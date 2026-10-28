% ------------------------------------------------------------------------
%
% The following polynomial dynamics are taken from [1, Example 5.2.3].
% The model is the Lotka-Volterra model of 5 variables. 
%
%   [1] Chen, Xin - Reachability Analysis of Non-Linear Hybrid Systems Using
%   Taylor Models, Dissertation 2015, URL:
%   https://publications.rwth-aachen.de/record/465295/files/465295.pdf
%
% ------------------------------------------------------------------------
function xdot = dynamics_LotkaVolterra5(x)

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);
x5 = x(5);

xdot0 = [   x1*(1 - (x1 + 0.85*x2 + 0.5*x5));
            x2*(1 - (x2 + 0.85*x3 + 0.5*x1));
            x3*(1 - (x3 + 0.85*x4 + 0.5*x2));
            x4*(1 - (x4 + 0.85*x5 + 0.5*x3));
            x5*(1 - (x5 + 0.85*x1 + 0.5*x4))];

xdot = xdot0;

end
