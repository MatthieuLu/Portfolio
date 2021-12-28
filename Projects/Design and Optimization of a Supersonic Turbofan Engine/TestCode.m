clear all
clc

%% Physical Variables
Qr = 45e06;
cp = 1004;
Ta = 217;
Pa = 7172;
Pe = Pa;
gamma = 1.4;
R = 287;
M = 1.7;
Pa_sl = 101.3e3;
Ta_sl = 288;

%% Fixed Engine Parameters
%Inlet/Diffuser
nd  = .95;
gd  = 1.4;
cpd = gd*R/(gd-1);
%Compressor (Polytropic)
npc = .9;
gpc = 1.37;
cppc = gpc*R/(gpc-1);
%Fan (Adiabatic)
nf  = .92;
gf  = 1.4;
cpf = gf*R/(gf-1);
%Burner 
nb  = .97;
bpr = .95; %Burner Pressure Recovery
gb  = 1.35;
cpb = gb*R/(gb-1);
%Turbine (Polytropic)
npt = .92;
gpt = 1.33;
cppt = gpt*R/(gpt-1);
%Jet Nozzle
nn1 = .98;
gn1 = 1.36;
cpn1 = gn1*R/(gn1-1);
%Fan Nozzle
nn2 = .99;
gn2 = 1.4;
cpn2 = gn2*R/(gn2-1);

%% Chosen Variables
D = 1.85;
beta = 0.5;
pi_f = 1.6;
T04 = 1650;
pi_c = 24;

%% Turbojet
u = M*sqrt(1.4*287*Ta);

%Inlet/Diffuser
T02 = Ta * (1+((gd-1)/2)*(M^2));
P02 = Pa * ((1 + nd*((T02/Ta)-1))^(gd/(gd-1)));
T0a = Ta * (1+((gamma-1)/2)*M^2);

%Mass Flow Rate
sigma_0 = P02/Pa_sl;
theta_0 = T02/Ta_sl;
A = pi*(D^2)/4;
mdot_inlet = 231.8*(sigma_0/sqrt(theta_0))*A;
mdot_a = mdot_inlet/(beta+1);

%Post-Compressor (03)
T03 = T02 * pi_c^((gpc-1)/(npc*gpc));
P03 = pi_c * P02; %polytropic

%Pre-Turbine (Burner, 04)
T04 = T04;
P04 = bpr * P03;

%Bypass (08)
T08 = T02 * (1 + (1/nf)*(pi_f^((gf-1)/gf) - 1));
P08 = P02 * pi_f;

%Post-Turbine (05)
T05 = T04 - (T03 - T02) - beta*(T08 - T0a);
pi_t = (T05/T04)^(gpt*npt/(gpt-1)); %polytropic
P05 = P04 * pi_t;

%Nozzle Inlet (06)
T06 = T05;
P06 = P05;

%Nozzle Exit (e)
ue = sqrt(2*nn1*(gn1/(gn1-1))*R*T06*(1-(Pa/P06)^((gn1-1)/gn1)));

%Fan Exit (ef)
uef = sqrt(2*nn2*(gn2/(gn2-1))*R*T08*(1-(Pa/P08)^((gn2-1)/gn2)));

%Fuel
f = (T04 - T03)/((nb*Qr/cpb)-T04);
mdot_f = f * mdot_a;

%ST
ST = ((1+f)*ue + beta*uef - (1+beta)*u)/1000

%TSFC
TSFC = f/ST

%Total Exhaust Mass Flow Rate
mdot_e = mdot_f + mdot_inlet;

%Thrust
T_bare = (((1+f)*ue + beta*uef - (1+beta)*u))*mdot_a;
T_eff = T_bare/(1.04+(0.01*beta^1.2))

%Thermal Efficiency
nth = ((0.5*ue.^2.*(mdot_f+mdot_a) + 0.5*uef.^2.*(mdot_inlet-mdot_a-mdot_f))-(0.5*u^2*mdot_inlet))./(mdot_f*Qr);

%Propulsion Efficiency
np = (T_eff*u)./((0.5*ue.^2.*(mdot_f+mdot_a) + 0.5*uef.^2.*(mdot_inlet-mdot_a-mdot_f))-(0.5*u^2*mdot_inlet));

%Overall Efficiency
no = nth.*np

