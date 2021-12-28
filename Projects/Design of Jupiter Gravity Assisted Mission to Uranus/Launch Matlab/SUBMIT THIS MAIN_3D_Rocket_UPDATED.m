clc;clear 
close all

%% (EDIT)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  ROCKET parameters (ATLAS V - modified)
% gross weight of fully fueled Atlas V: 305,143 [kg]

mfuel = 150000; % Mass of fuel for Atlas V (variable)

mf = 21054; % Dry mass of AtlasV (fixed) [kg]
mpayload = 12762.16; % SEC-III+fuel for mission+scientific payload [kg] (fixed)
m_wet    = mf+mfuel+mpayload;  %[kg]
m_dry    = mpayload;   %[kg]     Dry mass
d        = 3.8E-3;    %[km]     Rocket diameter
T        = 3827;    %[kN] Avg between Sea-level+Vaccum thrust
Isp      = 311.3;   %[s]   "                          " Specific impulse

% Launch Parameters
Az_po    = 108.2;     %[deg]    Azimuth at pitchover (inertial)
gam_po   = 84.621 ;    %[deg]    FPA at pitchover
z_po     = 0.100;     %[km]    Altitude at pitchover


% Launch Window Data
alphLS = 10.789;
alphaV = -42.926;
D_Lam    = (alphLS-alphaV)/180*pi; 
d_vinf   = -1.19396/180*pi; % Declination of v_inf

% Launch Site
phi_LS   =  28.57;  %[deg]    Launch site latitude
Lam_LS   = -80.65;  %[deg]    Launch site longitude
z_LS     =      3;  %[m]      Launch site altitude Actually correct for KSC

t_Coast  = 55;     %[s]      Coasting time after b/o, SpaceX launch was 14s
z_park   = 380;     %[km]     Parking orbit altitude 

R_Planet  = 778.6E6; %[km]  Radius of Juipter's orbit from sun

% Show plots ('y'=yes;'n'=no)
plots='y';


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% (DO NOT EDIT)
% Other Parameters 
RE   = 6378;       %[km]       Earth's radius
mu_E = 398600;     %[km^3/s^2] Gravitational parameter for Earth
g0   = mu_E/RE^2;  %[km/s^2]   Gravitational acceleration at SL
TE   = 86164.1;    %[s]        Earth's rotational p eriod
CD   = 0.5;        %[--]       Drag coefficient

R_Earth   = 149.6E6;
mu_Sun = 1.32712E11;
e_H = abs(R_Planet-R_Earth)/(R_Planet+R_Earth);
h_H = sqrt(R_Earth*mu_Sun*(1+e_H));
v_inf  = abs(h_H/R_Earth - sqrt(mu_Sun/R_Earth));
rp = RE+z_park;
vp_hyp = sqrt(v_inf^2+2*mu_E/rp);
e_hyp = vp_hyp^2*rp/mu_E - 1;
beta = acosd(1/e_hyp);
xi =  atan2d(sqrt((-cos(phi_LS/180*pi)*sin(d_vinf)+sin(phi_LS/180*pi)*cos(d_vinf)*cos(D_Lam))^2+(cos(d_vinf)*sin(D_Lam))^2),-sin(phi_LS/180*pi)*sin(d_vinf)-cos(phi_LS/180*pi)*cos(d_vinf)*cos(D_Lam)); 

%% Thrusting Phase

%Initial Condition
z0    = z_LS/1000;               %[km]
v0    = 0;                       %[km/s]
gam0  = gam_po/180*pi;           %[rad]
Az0   = Az_po/180*pi;            %[rad]
phi0  = phi_LS/180*pi;           %[rad]
Lam0  = Lam_LS/180*pi;           %[rad]
m0    = m_wet;                   %[kg]

% Mass and Time @ burnout
m_bo = m_wet-m_dry;
t_bo = m_bo/(T/(Isp*g0));

S01 = [z0,v0,gam0,Az0,phi0,Lam0,m0];

% Initial and final times
t0 = 0;  tf = t_bo;  %[s]

% Options for integration
Options = odeset('RelTol',1e-10,'AbsTol',1e-12,'OutputFcn',@func_collision,'Refine',200);
par1 = [Isp,T,d,CD,m_dry,RE,mu_E,z_po,gam_po/180*pi];

% Integration of Equations of Motion
[t1,S1] = ode45(@(t1,S1)func_ODE_3D_Rocket(t1,S1,par1), [t0,tf] ,S01 ,Options);

% Elaboration of results
z1   = S1(:,1);   % Altitude
v1   = S1(:,2);   % Velocity
gam1 = S1(:,3);   % FPA
Az1  = S1(:,4);   % Azimuth
phi1 = S1(:,5);   % Latitude
Lam1 = S1(:,6);   % Longitude 
m1   = S1(:,7);   % mass

Arc0 = haversine([phi1(1) Lam1(1)]*180/pi, [phi1(end) Lam1(end)]*180/pi);

% Gravity and Drag Losses
g1 = g0*(1+z1./RE).^-2; % Gravity Acceleration
rho1 = 1.225E9*exp(-z1./7.5); % Air Density
D1 = CD*d^2/4*pi*0.5*rho1.*v1.^2; % Drag

for i=1:size(t1,1)
    if t1(i,1)<=t_bo
        i_bo=i;
    end
end
if i_bo<size(t1,1)
    Nb=i_bo;
else
    Nb=size(t1,1)-1;
end
vG_loss = zeros(Nb,1);
vD_loss = zeros(Nb,1);    
vG_loss(1,1) = 0;
vD_loss(1,1) = 0;
for i=1:Nb
    vG_loss(i+1) = (t1(i+1)-t1(i))*(g1(i+1)*sin(gam1(i+1,1))+g1(i)*sin(gam1(i)))/2 + vG_loss(i);
    vD_loss(i+1) = (t1(i+1)-t1(i))*(D1(i+1)/m1(i+1)+D1(i)/m1(i))/2 + vD_loss(i);
end
Delta_vG = vG_loss(end);
Delta_vD = vD_loss(end);
Delta_V1 = v1(i_bo)+vG_loss(end)+vD_loss(end);

%% Coasting Phase (1)

% Initial and final times
t0 = t1(end);  tf = t0+t_Coast;  %[s]

% Options for integration
Options = odeset('RelTol',1e-10,'AbsTol',1e-12,'OutputFcn',@func_FPA,'Refine',200);
parC1 = [Isp,0,d,CD,m_dry,RE,mu_E,z_po,gam_po/180*pi];

% Integration of Equations of Motion
[tC1,SC1] = ode45(@(tC1,SC1)func_ODE_3D_Rocket(tC1,SC1,parC1), [t0,tf], S1(end,:), Options);
    
for i=1:size(SC1,1)-1
    if i>1 && SC1(i,3)*SC1(i+1,3)<0
        I=i;
        break
    else
        I=size(SC1,1);
    end
end

% Elaboration of results
zC1   = SC1(2:I,1);   % Altitude
vC1   = SC1(2:I,2);   % Velocity
gamC1 = SC1(2:I,3);   % FPA
AzC1  = SC1(2:I,4);   % Azimuth
phiC1 = SC1(2:I,5);   % Latitude
LamC1 = SC1(2:I,6);   % Longitude 
mC1   = SC1(2:I,7);   % mass

t   = [t1',tC1(2:I)']';
z   = [z1',zC1']';
v   = [v1',vC1']'+2*pi/TE*RE*cosd(phi_LS);
gam = [gam1',gamC1']';
Az  = [Az1',AzC1']';
phi = [phi1',phiC1']';
Lam = [Lam1',LamC1']';
m   = [m1',mC1']';

Arc1 = haversine([phi1(end) Lam1(end)]*180/pi, [phiC1(end) LamC1(end)]*180/pi);
%% Coasting Phase (2)
% Initial and final times
t_Coast2 = (xi+beta-(Arc0+Arc1))/180*pi/sqrt(mu_E)*(RE+z(end))^1.5;
t0 = t(end);  tf = t0 + t_Coast2;  %[s]

% Options for integration
Options = odeset('RelTol',1e-10,'AbsTol',1e-12,'Refine',200);
parC2 = [Isp,0,d,0,m_dry,RE,mu_E,z_po,gam_po/180*pi];

Delta_V2 = sqrt(mu_E/(RE+z(end)))-v(end);

% Integration of Equations of Motion
IC = [z(end),sqrt(mu_E/(RE+z(end))),0,Az(end),phi(end),Lam(end),m(end)];
[tC2,SC2] = ode45(@(tC2,SC2)func_ODE_3D_Rocket(tC2,SC2,parC2), [t0,tf], IC, Options);

% Elaboration of results
zC2   = SC2(2:end,1);   % Altitude
vC2   = SC2(2:end,2);   % Velocity
gamC2 = SC2(2:end,3);   % FPA
AzC2  = SC2(2:end,4);   % Azimuth
phiC2 = SC2(2:end,5);   % Latitude
LamC2 = SC2(2:end,6);   % Longitude 
mC2   = SC2(2:end,7);   % mass


t   = [t',tC2(2:end)']';
z   = [z',zC2']';
v   = [v',vC2']';
gam = [gam',gamC2']';
Az  = [Az',AzC2']';
phi = [phi',phiC2']';
Lam = [Lam',LamC2']';
m   = [m',mC2']';

Arc2 = haversine([phiC1(end) LamC1(end)]*180/pi, [phiC2(end) LamC2(end)]*180/pi);
%%
h = (z+RE).*v.*cos(gam);           % Angular momentum
v_park = sqrt(mu_E/(z_park+RE));
h_park = sqrt((z_park+RE)*mu_E);
Dv = Isp*g0*log(m_wet/m_dry)+2*pi/TE*RE*cosd(phi_LS);

thetaG = 2*pi/TE*t;
for i=1:size(t,1)
    R_BI(:,:,i) = [   cos(phi(i))*cos(thetaG(i)+Lam(i))  cos(phi(i))*sin(thetaG(i)+Lam(i))  sin(phi(i)); ...
                           -sin(thetaG(i)+Lam(i))            cos(thetaG(i)+Lam(i))              0      ;  ...
                     -sin(phi(i))*cos(thetaG(i)+Lam(i)) -sin(phi(i))*sin(thetaG(i)+Lam(i))  cos(phi(i))];

    r_I(:,i) = R_BI(:,:,i)'*[RE+z(i) 0 0]';
end


fprintf('At Burnout Time (t_bo = %.2f)    s:\n\n',t_bo)

fprintf('          Thrust Gains = +%.3f    km/s\n',Isp*g0*log(m_wet/m_dry))
fprintf('Earth''s Rotation Gains = +%.3f    km/s\n',2*pi/TE*RE*cosd(phi_LS))
fprintf('           Drag Losses = -%.3f    km/s\n',Delta_vD)
fprintf('        Gravity Losses = -%.3f    km/s\n',Delta_vG)
fprintf('---------------------------------------\n')
fprintf('         Burnout Speed = +%.3f    km/s\n\n\n',Isp*g0*log(m_wet/m_dry)+2*pi/TE*RE*cosd(phi_LS)-Delta_vD-Delta_vG)

fprintf('Delta V #1 = %.3f  km/s  (non-instantaneous)\n',Delta_V1)
if Delta_V2>0
    fprintf('Delta V #2 = %.3f  km/s  (PROGRADE - instantaneous to circularize orbit)\n',Delta_V2)
else
    fprintf('Delta V #2 = %.3f  km/s  (RETROGRADE - instantaneous to circularize orbit)\n',-Delta_V2)
end
fprintf('Total DV   = %.3f  km/s\n\n',Delta_V1+abs(Delta_V2))

fprintf('At Hyperbolic Orbit Injection Point:\n\n')

fprintf('Target Speed        v = %.3f    km/s\n',v_park)
fprintf('Target Altitude     z = %.1f    km\n',z_park)
fprintf('Target Ang. Mom.    h = %.0f    km^2/s\n',h_park)
fprintf('Target FPA        gam = 0.00      deg\n')
fprintf('Target Arc            = %.1f     deg\n\n',xi+beta) 

fprintf(' Final Speed        v = %.3f    km/s\n',v(end))
fprintf(' Final Altitude     z = %.1f    km\n',z(end))
fprintf(' Final Ang. Mom.    h = %.0f    km^2/s\n',h(end))
fprintf(' Final FPA        gam = %.2f     deg\n',gam(end)*180/pi)
fprintf(' Final Arc            = %.1f     deg\n\n',Arc0+Arc1+Arc2) 

fprintf(' Final LONGITUDE  Lam = %.2f    deg\n'  ,Lam(end)*180/pi)
fprintf(' Final LATITUDE   phi = %.2f    deg\n',phi(end)*180/pi)


if strcmp(plots,'y')
    % Altitude, Velocity, FPA, Angular momentum
    figure(2),hold on,grid on
    set(gcf,'units','normalized','position',[0 0 0.5 0.9])
    subplot(4,1,1),hold on,grid on
    xlabel('t   [s]')
    ylabel('z   [km]')
    plot(t,z,'color','b','linewidth',2.0)
    line([t(1),t(end)],[z_park,z_park],'linestyle','--','color','k')
    line([t_bo,t_bo],[0,max(z)],'linestyle','-.','color','k')
    subplot(4,1,2),hold on,grid on
    xlabel('t   [s]')
    ylabel('v   [km/s]')
    plot(t,v,'color','r','linewidth',2.0)
    line([t(1),t(end)],[v_park,v_park],'linestyle','--','color','k')
    line([t_bo,t_bo],[0,max(v)],'linestyle','-.','color','k')
    subplot(4,1,3),hold on,grid on
    xlabel('t   [s]')
    ylabel('\gamma   [deg]')
    plot(t,gam*180/pi,'color','g','linewidth',2.0)
    line([t(1),t(end)],[0,0],'linestyle','--','color','k')
    line([t_bo,t_bo],[0,max(gam*180/pi)],'linestyle','-.','color','k')
    hold off
    subplot(4,1,4),hold on,grid on
    xlabel('t   [s]')
    ylabel('h   [km^2/s]')
    plot(t,h,'color','m','linewidth',2.0)
    line([t(1),t(end)],[h_park,h_park],'linestyle','--','color','k')
    line([t_bo,t_bo],[0,max(h)],'linestyle','-.','color','k')
    hold off
    
    % Groundtrack
    figure(3)
    set(gcf,'units','normalized','position',[0 0 0.5 0.9])
    title('Groundtrack')
    geoplot(phi*180/pi,Lam*180/pi,'yo','markersize',2.0,'markeredgecolor','y','markerfacecolor','y')
    geobasemap colorterrain

    earth_plot(2*pi/TE,RE,r_I)
end