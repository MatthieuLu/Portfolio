clear
clc

syms ux1 uy1 ux2 uy2 ux3 uy3 ux4 uy4 fx1 fy1 fx2 fy2 fx3 fy3 fx4 fy4 E t real;
%reaction forces at each node
%fx1 = 22673.28;
fy1 = 0; 
fx2 = 0;
fy2 = -4723;
%fx3 = 0;
%fy3 = 18894.4;
%fx4 = -22673.28;
fy4 = 0;
t = 0.01;
E = 69*(10^9);

%truss parameters
L = 0.18;

B = 0.05;
H = 0.05;
A = (B + 2 * H) * t - 2 * (t^2);

forces = [fx1; fy1; fx2; fy2; fx3; fy3; fx4; fy4;];

%Displacements
ux1 = 0;
%uy1 unknown
%ux2 = 0; 
%uy2 unknown
ux3 = 0;
uy3 = 0;
ux4 = 0;
%uy4 unknown

u = [ux1; uy1; ux2; uy2; ux3; uy3; ux4; uy4;];

%beam angles
theta12 = 24.62;
theta13 = -24.62;
theta24 = -24.62;
theta34 = 24.62;

c12 = cosd(theta12);
s12 = sind(theta12);
c13 = cosd(theta13);
s13 = sind(theta13);
c24 = cosd(theta24);
s24 = sind(theta24);
c34 = cosd(theta34);
s34 = sind(theta34);

c = [c12  c13 c24 c34];
s = [s12 s13 s24 s34];

%stiffness matrices of each beam
K12 = [c12^2      c12 * s12    -(c12^2)        -(c12 * s12) 0 0 0 0;
       c12*s12    s12^2        -(c12 * s12)    -(s12^2) 0 0 0 0;
       -(c12^2)   -(c12 * s12) c12^2           c12*s12 0 0 0 0;
       -(c12*s12) -(s12^2)     c12*s12         s12^2 0 0 0 0;
       0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0;];

K13 = [c13^2      c13 * s13    0 0 -(c13^2)        -(c13 * s13) 0 0;
       c13*s13    s13^2        0 0 -(c13 * s13)     -(s13^2) 0 0;
       0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0;
       -(c13^2)   -(c13 * s13) 0 0 c13^2           c13*s13 0 0;
       -(c13*s13) -(s13^2)      0 0 c13*s13         s13^2 0 0;
       0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0;
       ];

K24 = [0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0;
       0 0 c24^2      c24 * s24  0 0   -(c24^2)        -(c24 * s24);
       0 0 c24*s24    s24^2     0 0   -(c24 * s24)    -(s24^2);
       0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0;
       0 0  -(c24^2)   -(c24 * s24) 0 0  c24^2           c24*s24;
       0 0 -(c24*s24) -(s24^2)     0 0 c24*s24         s24^2;];
 
K34 = [0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0;
       0 0 0 0 c34^2      c34 * s34    -(c34^2)        -(c34 * s34);
       0 0 0 0 c34*s34    s34^2        -(c34 * s34)    -(s34^2);
       0 0 0 0 -(c34^2)   -(c34 * s34) c34^2           c34*s34;
       0 0 0 0 -(c24*s34) -(s34^2)     c34*s34         s34^2;];

%combine the stiffness matrices
K = (K12 + K13 + K24 + K34);
vpa(K, 6);
K = (E * A / L) * K * u;
vpa(K, 6)

%set the matrix equal to boundary conditions and solve
for i = 1:8
    e(i, 1) = K(i) == forces(i);
end
vpa(e, 6)

sol = solve(e, [fx1 fx3 fy3 fx4 uy1 ux2 uy2 uy4]);
fleftx = vpa(sol.fx1, 5)
fbottomx = vpa(sol.fx3, 5)
fbottomy = vpa(sol.fy3, 5)
frightx = vpa(sol.fx4, 5)
dlefty = vpa(sol.uy1, 5)
dtopx = vpa(sol.ux2, 5)
dtopymeters = vpa(sol.uy2, 5)
drighty = vpa(sol.uy4, 5)

disp("--------------------------------");



   
 

   

