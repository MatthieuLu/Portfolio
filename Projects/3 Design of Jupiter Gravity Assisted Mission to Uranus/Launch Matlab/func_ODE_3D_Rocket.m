function Sdot = func_ODE_3D_Rocket(~,S,par)

Isp    = par(1);
Thrust = par(2);
d      = par(3);
CD     = par(4);
mf     = par(5);
RE     = par(6);
mu_E   = par(7);
z_po   = par(8);
gam_po = par(9);

z0     = 7.5;        %[km]       scale height of atmosphere
rho0   = 1.225E9;    %[kg/km^3]  atm. density at sea level
g0     = mu_E/RE^2;  %[km/s^2]   gravity acc. at sea level

z      = S(1,1);
v      = S(2,1);
gam    = S(3,1);
Az     = S(4,1);
phi    = S(5,1);
Lam    = S(6,1);
m      = S(7,1);

% Frontal Area
A = pi*(d/2)^2;

% Air Density
rho = rho0*exp(-z/z0);

% Drag
D = 0.5*CD*A*rho*v^2;

% m_dot
if m>mf
    T = Thrust;
else
    T = 0;
end

% Gravity acc. 
g = mu_E/(RE+z)^2;

if z<z_po && m>mf && gam>=gam_po
    % z_dot
    Sdot(1,1) = v;
    % v_dot
    Sdot(2,1) = T/m-D/m-g;
    % gam_dot
    Sdot(3,1) = 0;
    % Az_dot
    Sdot(4,1) = 0;
    % phi_dot
    Sdot(5,1) = 0;
    % Lambda_dot
    Sdot(6,1) = 0;
    % m_dot
    Sdot(7,1) = -T/(Isp*g0); 
else
    % z_dot
    Sdot(1,1) = v*sin(gam);
    % v_dot
    Sdot(2,1) = T/m-D/m -g*sin(gam);
    % gam_dot
    Sdot(3,1) = -(1/v)*(g - v^2/(RE+z))*cos(gam);
    % Az_dot
    Sdot(4,1) = (v*cos(gam)*sin(Az))/(cos(phi)*(RE+z))*sin(phi);
    % phi_dot
    Sdot(5,1) = (v*cos(gam)*cos(Az))/(RE+z);
    % Lambda_dot
    Sdot(6,1) = (v*cos(gam)*sin(Az))/(cos(phi)*(RE+z));
    % m_dot
    Sdot(7,1) = -T/(Isp*g0);
end
