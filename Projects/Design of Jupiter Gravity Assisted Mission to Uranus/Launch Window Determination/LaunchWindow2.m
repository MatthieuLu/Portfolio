clear
clc

J0 = 2459334.5; %Value of J0 for 2021 - 04 - 30
e = 0.0167; % Eccentricity 
ThetaEQ = 75.92305; % Degrees
Epsilon = 23.44; % Degrees
PhiLS = 28.5728722; % KSC launch site latitude, degrees
LamdaLS = -80.6489808; % KSC launch site longitude, degrees
VEquator = 0.460; % km/s
MuE = 398600; % km^3/s^2
rP = 380; % parking orbit radius, km
Vrot = VEquator*cos(PhiLS); % km/s
V = sqrt(MuE/rP); %km/s
% Part B: Theta_EQ = 75.92305 degrees at the vernal equinox of 2021/03/20 at 9:37:00 UT
for n = 0:96
    i(n+1,:) = [n];
    % Part C:
    DeltaUT(n+1,:) = [-24+(n/2)];
    JD(n+1,:) = [J0+(DeltaUT(n+1,1)/24)];
    TD(n+1,:) = [(JD(n+1,1)-2451545)/36525];
    % Part D: 
    aG(n+1,:) = [(100.4606184+(36000.77004*(TD(n+1,1)))+(0.000387933*((TD(n+1,1))^2))-((2.583*(10^(-8)))*((TD(n+1,1))^3)))+(360.98564724*(DeltaUT(n+1,1))/24)]; %Greenwich sidereal time 
    aLS(n+1,:) = [(aG(n+1,1)+LamdaLS)]; 
    % Part E: 
    Me(n+1,:) = [((100.464572+(35999.3725*TD(n+1,1)))-(102.937682+(0.32327364*TD(n+1,1))))];
    E(n+1,:) = [Me(n+1,1)+(e*sind(Me(n+1,1)))+(((e^2)/2)*sind(2*Me(n+1,1)))+(((e^3)/8)*((3*sind(3*Me(n+1,1)))-sind(Me(n+1,1))))];
    Theta(n+1,:) = [2*atand((sqrt((1+e)/(1-e)))*(tan(E(n+1,1)/2)))];
    UDeltav(n+1,:) = [atand(2*(((-sind(Epsilon)*sind(ThetaEQ)*sind(Theta(n+1,1)))-((e+cosd(Theta(n+1,1)))*sind(Epsilon)*cosd(ThetaEQ)))/(sqrt(((((cosd(ThetaEQ)*sind(Theta(n+1,1)))-((e+cosd(Theta(n+1,1)))*sind(ThetaEQ))))^2)+(((-cosd(Epsilon)*sind(ThetaEQ)*sind(Theta(n+1,1)))-((e+cosd(Theta(n+1,1)))*cosd(Epsilon)*cosd(ThetaEQ)))^2)))))]; % Uppercase delta v
    alphav(n+1,:) = [atand(2*(((-cosd(Epsilon)*sind(ThetaEQ)*sind(Theta(n+1,1)))-((e+cosd(Theta(n+1,1)))*cosd(Theta(n+1,1))*cosd(ThetaEQ)))/((cosd(ThetaEQ)*sind(Theta(n+1,1)))-((e+cosd(Theta(n+1,1)))*sind(ThetaEQ)))))];
    % Part F:
    dLambda(n+1,:) = [aLS(n+1,1)-alphav(n+1,1)];
    Az(n+1,:) = [(cosd(UDeltav(n+1,1))*sind(dLambda(n+1,1)))/((-cosd(PhiLS)*sind(UDeltav(n+1,1)))+(sind(PhiLS)*cosd(UDeltav(n+1,1))*cosd(dLambda(n+1,1))))];
    I(n+1,:) = [acosd(cosd(PhiLS)*sind(Az(n+1,1)))];
    Vrel(n+1,:) = [sqrt((V^2)+(Vrot^2)+(2*V*Vrot*sind(Az(n+1,1))))];
    AzEff(n+1,:) = [(Az(n+1,1)-acosd(sqrt(1-(((Vrot/Vrel(n+1,1))^2)*((cosd(Az(n+1,1)))^2)))))];
    M = [i, DeltaUT, JD, TD, aG, aLS, Me, E, Theta, UDeltav, alphav, dLambda, I, Vrel, Az, AzEff]; % Matrix with all desired variables at each value of n
end
format longG

disp( '[n, UT, JD, TD, aG, aLS, Me, E, Theta, UDeltav, alphav, dLambda, I, Vrel, Az, AzEff]:'); disp (M(37:38,:)); 
% AzMin<AzEff<AzMax at n = 36, so UT = 18 hr, adjusted for UTC-4 during
% daylight savings => UT = 14 hr. Therefore, our launch date and time is
% April 29th, 2021 at 2 PM UTC-4. 
