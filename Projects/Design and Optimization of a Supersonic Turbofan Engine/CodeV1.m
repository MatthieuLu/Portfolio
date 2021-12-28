clear all
close all
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
mdot_a = mdot_inlet ./(beta+1);

%Min ST and Max TSFC
thrust_min = 80
D_max = 1.8;
mdot_max = 231.8*(sigma_0/sqrt(theta_0))*(pi*(D_max/2)^2);
T04_max = 1800;
T03_max = 1000; %max from design specs
f_max = (T04_max - T03_max)/((nb*Qr/cpb)-T04_max);
mdot_a_max = mdot_max/(1+beta);
mdot_f_max = mdot_a_max*f_max;
ST_min = thrust_min/(mdot_a_max*(1+beta))
TSFC_max = mdot_f_max/thrust_min

figure
hold on
for i = 1:9
    T04(i) = 1400 + 50*(i-1);
    
    for j = 1:13
        pi_c(i,j) = 16 + 2*(j-1);
        
        %Post-Compressor (03)
        T03(i,j) = T02 * pi_c(i,j)^((gpc-1)/(npc*gpc));
        P03 = pi_c(i,j) * P02; %polytropic
        
        %Pre-Turbine (Burner, 04)
        T04(i) = T04(i);
        P04 = bpr * P03;
        
        %Bypass (08)
        T08 = T02 * (1 + (1/nf)*(pi_f^((gf-1)/gf) - 1));
        P08 = P02 * pi_f;
        
        %Post-Turbine (05)
        T05 = T04(i) - (T03(i,j) - T02) - beta*(T08 - T0a);
        pi_t = (T05/T04(i)).^(gpt*npt/(gpt-1)); %polytropic
        P05 = P04 * pi_t;
        
        %Nozzle Inlet (06)
        T06 = T05;
        P06 = P05;
        
        %Nozzle Exit (e)
        ue(i,j) = sqrt(2*nn1*(gn1/(gn1-1))*R*T06*(1-(Pa/P06)^((gn1-1)/gn1)));
        
        %Fan Exit (ef)
        uef(i,j) = sqrt(2*nn2*(gn2/(gn2-1))*R*T08*(1-(Pa/P08)^((gn2-1)/gn2)));
        
        %Fuel
        f = (T04(i) - T03(i,j))/((nb*Qr/cpb)-T04(i));
        mdot_f(i,j) = f * mdot_a;
        
        %ST
        ST(i,j) = ((1+f)*ue(i,j) + beta*uef(i,j) - (1+beta)*u)/1000;
        
        %TSFC
        TSFC(i,j) = f/ST(i,j);
        
        %Total Exhaust Mass Flow Rate
        mdot_e(i,j) = mdot_f(i,j) + mdot_inlet;
        
        %Thrust
        T_bare(i,j) = (((1+f)*ue(i,j) + beta*uef(i,j) - (1+beta)*u))*mdot_a;
        T_eff(i,j) = T_bare(i,j)/(1.04+(0.01*beta^1.2));
        
        plot(ST(:,j), TSFC(:,j), 'b')

    end
    text(ST(i,1),TSFC(i,1),['T04 = ' num2str(T04(i))])
    plot(ST(i, :), TSFC(i, :), 'r')
    
end

%Thermal Efficiency
nth = ((0.5*ue.^2.*(mdot_f+mdot_a) + 0.5*uef.^2.*(mdot_inlet-mdot_a-mdot_f))-(0.5*u^2*mdot_inlet))./(mdot_f*Qr);

%Propulsion Efficiency
np = (T_eff*u)./((0.5*ue.^2.*(mdot_f+mdot_a) + 0.5*uef.^2.*(mdot_inlet-mdot_a-mdot_f))-(0.5*u^2*mdot_inlet));

%Overall Efficiency
no = nth.*np;

for x = 1:j
    text(ST(i,x),TSFC(i,x),['\pi_c = ' num2str(pi_c(i,x))])
end
title(['Carpet Plot of ST vs TSFC with a Bypass Ratio of ' num2str(beta) ' and a Fan Compression Ratio of ' num2str(pi_f)])
xlabel('ST (kN s / kg)')
ylabel('TSFC (kg / kN s)')
xline(ST_min,'--','ST_m_i_n')
yline(TSFC_max,'--','TSFC_m_a_x')

hold off

%Thrust Graph
figure
contour(pi_c(1,:),T04,T_eff,'ShowText','on')
title('Contour Plot of Thrust wrt T_0_4 and \pi_c')
xlabel('\pi_c')
ylabel('T_0_4')

%TSFC and ST Contours
figure
subplot(1,2,1)
contour(pi_c(1,:),T04,TSFC,'ShowText','on')
title('TSFC wrt \pi_c and T_0_4')
xlabel('\pi_c')
ylabel('T_0_4')

subplot(1,2,2)
contour(pi_c(1,:),T04,ST,'ShowText','on')
title('ST wrt \pi_c and T_0_4')
xlabel('\pi_c')
ylabel('T_0_4')

%nth, np, and no Contours
figure
subplot(2,2,1)
contour(pi_c(1,:),T04,nth,'ShowText','on')
title('Thermal Efficiency wrt \pi_c and T_0_4')
xlabel('\pi_c')
ylabel('T_0_4')

subplot(2,2,2)
contour(pi_c(1,:),T04,np,'ShowText','on')
title('Propulsion Efficiency wrt \pi_c and T_0_4')
xlabel('\pi_c')
ylabel('T_0_4')

subplot(2,2,3)
contour(pi_c(1,:),T04,no,'ShowText','on')
title('Overall Efficiency wrt \pi_c and T_0_4')
xlabel('\pi_c')
ylabel('T_0_4')


%% Beta and Pi_f Contours
clear T03

[beta,pi_f] = meshgrid(0:.1:10,1:0.1:1.6);
T04 = 1750;
pi_c = 20;
D = 1.8;

%Post-Compressor (03)
T03 = T02 * pi_c^((gpc-1)/(npc*gpc));
P03 = pi_c * P02; %polytropic

%Pre-Turbine (Burner, 04)
T04 = T04;
P04 = bpr * P03;

%Bypass (08)
T08 = T02 * (1 + (1/nf).*(pi_f.^((gf-1)/gf) - 1));
P08 = P02 .* pi_f;

%Post-Turbine (05)
T05 = T04 - (T03 - T02) - beta.*(T08 - T0a);
pi_t = (T05./T04).^(gpt*npt/(gpt-1)); %polytropic
P05 = P04 .* pi_t;

%Nozzle Inlet (06)
T06 = T05;
P06 = P05;

%Nozzle Exit (e)
ue = sqrt(2*nn1*(gn1/(gn1-1))*R.*T06.*(1-(Pa./P06).^((gn1-1)/gn1)));

%Fan Exit (ef)
uef = sqrt(2*nn2*(gn2/(gn2-1))*R.*T08.*(1-(Pa./P08).^((gn2-1)/gn2)));

%Fuel
f = (T04 - T03)/((nb*Qr/cpb)-T04);
mdot_f = f * mdot_a;

%ST
ST = ((1+f).*ue + beta.*uef - (1+beta).*u)./(1000*(1+beta));

%TSFC
TSFC = f./real(((1+f).*ue + beta.*uef - (1+beta).*u)/1000);

figure
subplot(2,1,1);
contour(beta,pi_f,ST,'ShowText','on');
title(['ST @ ' num2str(T04) 'K and a \pi_c of ' num2str(pi_c)])
xlabel('\beta');
ylabel('\pi_f');
subplot(2,1,2);
contour(beta,pi_f,TSFC,'ShowText','on');
xlabel('\beta');
ylabel('\pi_f');
title(['TSFC @ ' num2str(T04) 'K and a \pi_c of ' num2str(pi_c)])



