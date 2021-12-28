clear
clc

%Variables: Sun
mu_S = 132712000000;

%Variables: Earth
mu_E = 398600;
r_E = 6371;
r_ES = 149.6E6;
T_ES = 365.256*24*60*60;
V_ES = sqrt(mu_S/r_ES);

%Variables: Jupiter
mu_J = 126686000;
r_J = 71490;
r_JS = 778.6E6;
T_JS = 11.86*365*24*60*60;
V_JS = sqrt(mu_S/r_JS)

%Variables: Uranus
mu_U = 5794000;
r_U = 25560;
r_US = 2.872E9;
T_US = 84.01*365*24*60*60;
r_U_p = 200420;
r_U_a = 3*r_U_p;
V_US = sqrt(mu_S/r_US);

%Gravity Assist
%Variables: Gravity Assist
DV1 = 5.242; %km/s
DV2 = 2.257; %km/s
DVtot = 7.499; %km/s
t = 10*365; %mission time in days

%Phasing Angle
Phi_0 = 180*(1-((r_ES + r_JS)/(2*r_JS))^1.5) %degrees

%Interplanetary Travel

%Hohmann Transfer (Earth to Jupiter)
r_Hohmann_p = r_ES; %periapses radius for Hohmann transfer
r_Hohmann_a = r_JS; %apoapses radius for Hohmann transfer
mu_Hohmann = mu_S;
a_Hohmann = (r_Hohmann_p + r_Hohmann_a)/2 %semimajor axis for Hohmann Transfer
En_Hohmann = -mu_Hohmann/(2*a_Hohmann) %Energy of Elliptical Orbit
V_SC_J = sqrt((En_Hohmann + (mu_Hohmann/r_Hohmann_a))*2) %V of SC at edge of Jupiter SOI (heliocentric RF)
V_SC_E = sqrt((En_Hohmann + (mu_Hohmann/r_Hohmann_p))*2) %V of SC at edge of Earth SOI (heliocentric RF)

%Departure
r_E_parking = 400+r_E; %parking orbit radius
V_Departure_circ = sqrt(mu_E/r_E_parking); %parking orbit velocity
V_Departure_hyp = V_Departure_circ + DV1; %hyperbolic orbit velocity @ perigee
V_Departure_inf = V_SC_E - V_ES %V_inf for hyperbolic orbit departure
a_Departure_hyp = mu_E / (V_Departure_inf^2)
e_Departure_hyp = (r_E_parking / a_Departure_hyp) + 1 %eccentricity of hyperbolic orbit departure
ar_Departure_hyp = a_Departure_hyp*sqrt(e_Departure_hyp^2 - 1) %aiming radius for departure orbit
ta_Departure_hyp = 2*asind(1/e_Departure_hyp) %turning angle for departure orbit

%Arrival
r_U_capture = 240000;
a_Capture = (r_U_p + r_U_a)/2
En_Capture = -mu_U/(2*a_Capture) %Energy of capture orbit
V_Capture_ell_p = sqrt((En_Capture + (mu_U/r_U_p))*2)
V_Capture_ell_a = sqrt((En_Capture + (mu_U/r_U_a))*2)
V_Arrival_hyp_p = V_Capture_ell_p + DV2 %V @ hyperbolic perigee




%Fuel
%Rocket:
syms I_sp thrust m_0 real
% I_sp = ;
% thrust = ;
% m_0 = ;
m_f = 1500;
g_0 = 9.8/1000; %km/s^2

