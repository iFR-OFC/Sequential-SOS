% ------------------------------------------------------------------------
%
% The following polynomial dynamics are taken from [1, Example 3.3.9] on
% page 71. It is an extension of the Moore-Greitzer jet engine model.
%
%
%   [1] Chen, Xin - Reachability Analysis of Non-Linear Hybrid Systems Using 
%   Taylor Models, Dissertation 2015, URL:  
%   https://publications.rwth-aachen.de/record/465295/files/465295.pdf
%
% ------------------------------------------------------------------------
function xdot = dynamics_JetEngine(x,u)

y = x(2);
x = x(1);

xdot0 = [- y - 1.5*x^2 - 0.5 *x^3 - 0.5;
          3*x - y];


q6 = [0.2*x^3*y^3;
      0.2*x^3*y^3];

q9 = [- 0.1*x^3*y^6 - 0.1*x^6*y^3;
      - 0.1*x^3*y^6 - 0.1*x^6*y^3];

xdot = xdot0 + q6 + q9 + u;


end