function xdot = copter_dynamics(x,u)


x = x*1;
u = u*1;

K0 = 0.89/1.4;
d0 = 70;
d1 = 17;
n0 = 55;
gn = 9.81;

cos_app = -0.498*x(5)^2+1;
sin_app = -0.166*x(5)^3+x(5);


xdot = [x(3);
         x(4);
        u(1)*K0*sin_app;
        u(1)*K0*cos_app - gn;
         x(6);
         -d0*x(5)-d1*x(6)+n0*u(2)];



% we need to shift the dynamics 
% u0 = 

end
