clear
clc
sympref('FloatingPointOutput',true);

%Constants
gamma = 1.4;
R = 287;
Cp = (gamma*R)/(gamma-1);

%Point 1
syms v1 P1 T1 P01 T01 M1 A1 a1 rho1
%v1  = 0;
P1  = 101325;
T1  = 285;
P01 = 101325;
T01 = 285;
%M1  = 0;
A1  = .5;
a1  = sqrt(gamma*R*T1);
rho1= 1.225;
rho01 = 1.225;

%Point 1
syms v2 P2 T2 P02 T02 M2 A2 a2 rho2
%v2  = 0;
%P2  = 101000;
%T2  = 296;
%P02 = 101000;
%T02 = 296;
% M2  = 0.025;
A2  = A1;
%a2  = sqrt(gamma*R*T1);
%rho2= 1.225;

%Point Crit
syms vc Pc Tc P0c T0c Mc Ac ac rhoc
%vc  = 0;
%Pc  = 101000;
%Tc  = 296;
%P0c = 101000;
%T0c = 296;
Mc  = 1;
%Ac  = ;
%ac  = sqrt(gamma*R*T1);
%rhoc= 1.225;

%Point 3
syms v3 P3 T3 P03 T03 M3 A3 a3 rho3
%v3  = 0;
% P3  = 17962.7;
% T3  = 164.476;
% P03 = 140553;
% T03 = 296.033;
M3  = 2;
A3  = .25;
%a3  = sqrt(gamma*R*T1);
%rho3= 1.225;
Ac = A3/1.687;

%Point Post Shock
syms vs2 Ps2 Ts2 P0s2 T0s2 Ms2 As2 as2 rhos2 s2
%vs2  = 0;
%Ps2  = 101000;
%Ts2  = 296;
%P0s2 = 101000;
%T0s2 = 296;
Ms2  = 0.6279;
As2  = 0.2056;
%as2  = sqrt(gamma*R*T1);
%rhos2= 1.225;
Asc = 0.1777;

%Point 4
syms v4 P4 T4 P04 T04 M4 A4 a4 rho4
%v4  = 0;
P4  = 101325;
T4  = 285;
% P0e = 101325;
% T0e = 296;
% M4  = .013225;
A4  = .75;
%a4  = sqrt(gamma*R*T1);
rho4= 1.225;

%Section 4
AMR4 = A4/Asc == ((5 + M4^2)^3)/(216*M4);
M4 = min(double(vpasolve(AMR4,M4,[0 Inf])));
PMR4 = (1+(((gamma-1)/2)*(M4^2)))^(gamma/(gamma-1));
P04 = PMR4*P4;
TMR4 = 1+(((gamma-1)/2)*(M4^2));
T04 = TMR4*T4;
RMR4 = (1+(((gamma-1)/2)*(M4^2)))^(1/(gamma-1));
rho04 = rho4*RMR4;
a4 = sqrt(gamma*R*T4);
v4 = M4*a4;

%Section Post Shock
P0s2 = P04;
T0s2 = T04;
rho0s2 = rho04;
PMRs2 = (1+(((gamma-1)/2)*(Ms2^2)))^(gamma/(gamma-1));
Ps2 = P0s2/PMRs2;
TMRs2 = 1+(((gamma-1)/2)*(Ms2^2));
Ts2 = T0s2/TMRs2;
RMRs2 = (1+(((gamma-1)/2)*(Ms2^2)))^(1/(gamma-1));
rhos2 = rho0s2/RMRs2;
as2 = sqrt(gamma*R*Ts2);
vs2 = Ms2*as2;

%Section 3
P0shock = ((((gamma+1)*(M3^2))/(2+((gamma-1)*(M3^2))))^(gamma/(gamma-1)))*(((gamma+1)/((2*gamma*(M3^2))-(gamma-1)))^(1/(gamma-1)));
P03 = P0s2/P0shock;
PMR3 = (1+(((gamma-1)/2)*(M3^2)))^(gamma/(gamma-1));
P3 = P03/PMR3;
Tshock = (1+((2*gamma*((M3^2)-1))/(gamma+1)))*((2+((gamma-1)*(M3^2)))/((gamma+1)*(M3^2)));
T3 = Ts2/Tshock;
TMR3 = 1+(((gamma-1)/2)*(M3^2));
T03 = T3*TMR3;
Rshock = ((gamma+1)*(M3^2))/(2+((gamma-1)*(M3^2)));
rho3 = rhos2/Rshock;
RMR3 = (1+(((gamma-1)/2)*(M3^2)))^(1/(gamma-1)); 
rho03 = rho3*RMR3;
a3 = sqrt(gamma*R*T3);
v3 = M3*a3;

%Section 2
P02 = P03;
T02 = T03;
rho02 = rho03;
AMR2 = A2/Ac == ((5 + M2^2)^3)/(216*M2);
M2 = min(double(vpasolve(AMR2,M2,[0 Inf])));
PMR2 = (1+(((gamma-1)/2)*(M2^2)))^(gamma/(gamma-1));
P2 = P02/PMR2;
TMR2 = 1+(((gamma-1)/2)*(M2^2));
T2 = T02/TMR2;
RMR2 = (1+(((gamma-1)/2)*(M2^2)))^(1/(gamma-1)); 
rho2 = rho02/RMR2;
a2 = sqrt(gamma*R*T2);
v2 = M2*a2;

%Mass Flow Rate
GammaM = sqrt(gamma*(2/(gamma+1))^((gamma+1)/(gamma-1)));
qM3 = M3*((2/(gamma+1))*(1+(((gamma-1)/2)*(M3^2))))^(-(gamma+1)/(2*(gamma-1)));
mdotflow = GammaM*qM3*A3*(P03/sqrt(R*T03));

%Fan Properties
TotalPressureRatio = P02/P01;
Power = (0.5*mdotflow*((M4*a4)^2))/.92
Power = Power*.001341022; %horsepower
ForceofThrust = mdotflow*v4;
Re = (rho3*v3*.5)/(1.789*(10^(-5)))
FanProperties = table(TotalPressureRatio,Power,ForceofThrust,mdotflow,Re)

%Section 1
qM1 = (mdotflow*sqrt(R*T01))/(GammaM*P01*A1);
AMR1 = 1/qM1 == ((5 + M1^2)^3)/(216*M1);
M1 = min(double(vpasolve(AMR1,M1,[0 Inf])));
P1 = P01/((1+(((gamma-1)/2)*(M1^2)))^(gamma/(gamma-1)));
T1 = T01/(1+(((gamma-1)/2)*(M1^2)));
a1 = sqrt(gamma*R*T1);
v1 = M1*a1;


%Display
Point = [1,2,3,0,4]';
MachNumber = [M1,M2,M3,Ms2,M4]';
Area = [A1,A2,A3,As2,A4]';
StagPres = [P01,P02,P03,P0s2,P04]';
StagTemp = [T01,T02,T03,T0s2,T04]';
StagDens = [rho01,rho02,rho03,rho0s2,rho04]';
Velocity = [v1,v2,v3,vs2,v4]';
SoundSpeed = [a1,a2,a3,as2,a4]';
Pressure = [P1,P2,P3,Ps2,P4]';
Temperature = [T1,T2,T3,Ts2,T4]';
Results = table(Point,StagPres,StagTemp,StagDens,MachNumber,Area,Velocity,SoundSpeed,Pressure,Temperature)

%--------------Tunnel Calculations--------------%


%Length of Tunnel
x = [0:0.01:9];
inlength = 3;
outlength = 4;
Mt2 = 0.175066
Mt2s = Ms2
At1=0.1482;
A_ratio=1/ 0.72087386; %p02/p01
At2=At1*A_ratio
Mt2=1.750666667;
Mt2s= 0.62792868;

M2=0.17513924050632912;
M3=2; 

A_exit=.75;
AstarS=At2/1.157
M_exit=0.1388592375366569;

A3=.5*.5;
At=A3/1.687; %area of throat

%define M(x)
x=linspace(0,9,10000);
M=[];
inlength=3;
outlength=4;

%Ma Number Piecewise
for i=1:length(x);
   if x(i) < inlength;
       M(i)=(M3-M2)/2*cos(1/3*pi*x(i)+pi)+(M3-M2)/2+M2;
   elseif x(i) >=inlength & x(i)<(inlength+1);
       M(i)=M3;
   elseif x(i) >=(inlength+1) & x(i)<=(inlength+1+outlength);
       M(i)=(M3-Mt2)/2*cos(1/4*pi*x(i)-pi)-(M3-Mt2)/2+M3;
   else
       M(i)=(Mt2s-M_exit)/2*cos(pi*x(i))-(Mt2s-M_exit)/2+Mt2s;
   end   
end

for n = 1:length(x);
    if x(n) <= (8)
        %Area
        Aratio(n) = (((2/(gamma+1))*(1+((gamma-1)/2)*(M(n)^2)))^((gamma+1)/(2*(gamma-1))))/M(n);
        A(n) = Aratio(n)*Ac;
        halfheight(n) = A(n);
        %Pressure
        Pratio(n) = (1+(((gamma-1)/2)*(M(n)^2)))^(gamma/(gamma-1));
        P(n) = (P03/Pratio(n))/1000;
        %Temperature
        Tratio(n) = 1+(((gamma-1)/2)*(M(n)^2));
        T(n) = T03/Tratio(n);
    else
        %Area
        Aratio(n) = (((2/(gamma+1))*(1+((gamma-1)/2)*(M(n)^2)))^((gamma+1)/(2*(gamma-1))))/M(n);
        A(n) = Aratio(n)*Asc;
        halfheight(n) = A(n);
        %Pressure
        Pratio(n) = (1+(((gamma-1)/2)*(M(n)^2)))^(gamma/(gamma-1));
        P(n) = (P04/Pratio(n))/1000;
        %Temperature
        Tratio(n) = 1+(((gamma-1)/2)*(M(n)^2));
        T(n) = T04/Tratio(n);
    end
    %Velocity
    v(n) = sqrt(R*gamma*T(n))*M(n);
    %Density
    rho(n) = rho03/(1+(((gamma-1)/2)*(M(n)^2)))^(1/(gamma-1));
end

%--------------Validation--------------%
for n = 1:length(x)
    qM(n) = 1/((((2/(gamma+1))*(1+((gamma-1)/2)*(M(n)^2)))^((gamma+1)/(2*(gamma-1))))/M(n));
    if x(n) <= (inlength+1)
        mdot(n) = GammaM*qM(n)*A(n)*(P03/sqrt(R*T03));
    else
        mdot(n) = GammaM*qM(n)*A(n)*(P04/sqrt(R*T04));
    end
    h0(n) = (Cp*T(n))+(0.5*(v(n)^2));
end

figure(1)
tiledlayout(2,1);

nexttile;
plot(x,mdot);
axis([0 9 25 75])
title('Mass Flow Rate wrt length');
xlabel('Length (m)');
ylabel('Mass Flow Rate');

nexttile;
plot(x,h0);
axis([0 9 200000 400000])
title('Energy wrt length');
xlabel('Length (m)');
ylabel('Energy');

%--------------Results--------------%
figure(2)
tiledlayout(3,2);

nexttile;
plot(x,M);
% axis([0 8 0 inf])
title('Mach Number wrt length');
xlabel('Length (m)');
ylabel('Ma');

nexttile;
plot(x,halfheight);
% axis([0 8 0 inf])
title('Shape of Wind Tunnel');
xlabel('Length (m)');
ylabel('Half Height (m), Depth = 0.5m');

nexttile;
plot(x,v);
% axis([0 8 30 220])
title('Velocity wrt length');
xlabel('Length (m)');
ylabel('Velocity (m/s)');

nexttile;
plot(x,P);
% axis([0 8 75 105])
title('Pressure wrt length');
xlabel('Length (m)');
ylabel('Presure (kPa)');

nexttile;
plot(x,T);
% axis([0 8 275 300])
title('Temperature wrt length');
xlabel('Length (m)');
ylabel('Temperature (K)');

nexttile;
plot(x,rho);
% axis([0 8 1 1.25])
title('Density wrt length');
xlabel('Length (m)');
ylabel('Density (kg/m^3)');

x = x';
halfheight = halfheight';
