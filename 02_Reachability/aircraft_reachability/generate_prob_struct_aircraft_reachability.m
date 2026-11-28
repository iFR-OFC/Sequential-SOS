% *****+*******************************************************************
%
%   Short Description:  Implements the problem struct for reachability
%                       problem using an Aircraft dynamical system from
%                       [1].
%
%
%   Reference: [1]	H. Yin, M. Arcak, A. Packard, und P. Seiler, „Backward
%                   Reachability for Polynomial Systems on a Finite
%                   Horizon", IEEE Transactions on Automatic Control, Bd.
%                   66, Nr. 12, S. 6025–6032, Dez. 2021, doi:
%                   10.1109/TAC.2021.3056611.
%
%   License: GNU GENERAL PUBLIC LICENSE Version 3
%
% *****+*******************************************************************

function [probStruct] = generate_prob_struct_aircraft_reachability(degV,degK)

x = casos.Indeterminates('x',4,1);
u = casos.Indeterminates('u',1,1);

t  = casos.PS('t');

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

d2r = pi/180;
Dmax = diag([20 20*d2r 50*d2r 20*d2r]);

f0 = [ ...
      - 0.01955829859645207*x1^2*x2 + 0.0006116050110168922*x1^2*x3 - 0.4597074905072323*x1*x2^2 ...
      - 0.02143363124979007*x1*x2*x3 + 0.0913633506074555*x2^3 + 0.0104276150041391*x2^2*x4 ...
      - 0.0104276150041391*x2*x4^2 + 0.003475871668046367*x4^3 - 0.00975287614756193*x1^2 ...
      - 0.08801234368403432*x1*x2 - 0.001251382140034897*x1*x3 - 0.5201826057909725*x2^2 ...
      - 0.04982793763401992*x2*x3 + 0.00180695498489715*x3^2 - 0.04388794266402869*x1 ...
      + 0.07194181873717953*x2 - 0.00290915666614335*x3 - 0.1711592037553279*x4;...
      0.05879789171572614*x1^3 + 0.6755306062430397*x1^2*x2 + 0.07878176650294093*x1^2*x3 ...
      + 0.6715746030304725*x1*x2^2 - 0.03631259653379534*x1*x2*x4 - 0.006619604539373272*x1*x3^2 ...
      + 0.01815629826689767*x1*x4^2 - 0.1817374041405456*x2^3  + 0.1364168103441202*x1^2 ...
      - 1.371702883566865*x1*x2 + 0.005525130555690897*x1*x3 + 1.479555253341581*x2^2 ...
      + 0.002672981895904953*x2*x3 + 0.07915793925313672*x2*x4 + 0.0135837161533212*x3^2 ...
      - 0.03957896962656836*x4^2  - 0.5767816606949114*x1 - 3.236059303865032*x2 + 2.30669948822443*x3;...
      - 3.582346855286648*x1^2*x2 + 0.9194217971573437*x1^2*x3 + 2.279359758657476*x1*x2^2 ...
      - 0.3522980647489417*x2^3 - 0.07456233935706458*x1^2 - 16.12056084878992*x1*x2 ...
      - 1.881194554322757*x1*x3 + 2.56427972848966*x2^2  - 0.33553052710679*x1 ...
      - 18.13563095488866*x2 - 4.373316114187153*x3;...
      2.5*x3];
g = [ ...
     - 0.01200284890938563*x1^2- 0.0672488718061946*x1*x2- 0.05401282009223534*x1 ...
     - 0.07565498078196893*x2- 0.06076442260376476; ...
     0.190652302942119*x1^2 + 0.1033082169845372*x1*x2 - 0.3900865469718997*x1 ...
     + 0.240166275746567*x2 - 0.9068555816726751; ...
     -13.57875281932966*x1^2 + 14.75731515054161*x1*x2 - 61.10438768698349*x1 ...
     + 16.60197954435931*x2 - 68.74243614785642; ...
     0
    ];

f = f0 + g*u;

P = [4.54064415441163	0.139789072082587	-0.0107315729582360	-0.406387155887787
     0.139789072082587	0.165595957230149	-0.00115647767230559	-0.0841739454608939
    -0.0107315729582360	-0.00115647767230559	0.00732884319609208	0.00922486739143782
   -0.406387155887787	-0.0841739454608939	0.00922486739143782	0.514580693591350];

g  = x'*P*x-1;

gT = x'*Dmax'*blkdiag(1/(4)^2, 1/(pi/30)^2, 1/(pi/15)^2, 1/(pi/30)^2)*Dmax*x-1;

umax =   10*d2r - 0.0489;
umin = -(10*d2r + 0.0489);

T  = 3;

V  = casos.PS.sym('v',monomials([x;t],degV));

kappa = casos.PS.sym('kappa',monomials([x;t],degK),[length(u) 1]);

s1 = casos.PS.sym('s1',monomials([x;t],0:2));
s2 = casos.PS.sym('s2',monomials([x;t],0:2));
s3 = casos.PS.sym('s3',monomials([x;t],0:2),[length(u) 1]);
s4 = casos.PS.sym('s4',monomials([x;t],0:2),[length(u) 1]);
s5 = casos.PS.sym('s5',monomials([x;t],0:2),[length(u) 1]);
s6 = casos.PS.sym('s6',monomials([x;t],0:2),[length(u) 1]);
s7 = casos.PS.sym('s7',monomials([x;t],0:2));
s8 = casos.PS.sym('s8',monomials([x;t],0:2));
s9 = casos.PS.sym('s9',monomials(x,0:2));

b  = casos.PS.sym('b');

probStruct.x    = x;
probStruct.u    = u;
probStruct.t    = t;

probStruct.f    = f;

probStruct.g    = g;
probStruct.gT   = gT;

probStruct.T    = T;

probStruct.umin = umin;
probStruct.umax = umax;

probStruct.V    = V;

probStruct.s1   = s1;
probStruct.s2   = s2;
probStruct.s3   = s3;
probStruct.s4   = s4;
probStruct.s5   = s5;
probStruct.s6   = s6;
probStruct.s7   = s7;
probStruct.s8   = s8;
probStruct.s9   = s9;

probStruct.b    = b;

probStruct.kappa = kappa;

probStruct.Name = 'Aircraft';

end
