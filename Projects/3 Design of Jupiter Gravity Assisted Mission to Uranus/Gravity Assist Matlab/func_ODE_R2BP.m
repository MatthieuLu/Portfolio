% Equations of motion

function pdot = func_ODE_R2BP(~,p,mu)
x  = p(1);
y  = p(2);
Vx = p(3);
Vy = p(4);

r = sqrt(x^2+y^2);

x_dot  = Vx;
y_dot  = Vy;
Vx_dot = -mu/r^3*x;
Vy_dot = -mu/r^3*y;

pdot = [x_dot,y_dot,Vx_dot,Vy_dot]';