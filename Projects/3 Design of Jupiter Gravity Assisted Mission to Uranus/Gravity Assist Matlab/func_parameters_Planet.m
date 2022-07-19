% Planetary parameters

function [r,R,mu,T,n,v,r_SOI] = func_parameters_Planet(Planet)

mu_Sun = 1.3271244e11; %[km^3/s^2]

if strcmp(Planet,'Mercury')
    R = 2440;
    mu = 22032;
    r = 57.91e6;
elseif strcmp(Planet,'Venus')
    R = 6052;
    mu = 324900;
    r = 108.2e6;
elseif strcmp(Planet,'Earth')    
    R = 6378;
    r = 149.9E6;
    mu = 398600;
elseif strcmp(Planet,'Mars')   
    R = 3396;
    mu = 42828;
    r = 227.9e6;
elseif strcmp(Planet,'Jupiter')
    R = 71490;
    mu = 126686000;
    r = 778.6e6;
elseif strcmp(Planet,'Saturn') 
    R = 60270;
    mu = 37931000;
    r = 1433e6;
elseif strcmp(Planet,'Uranus') 
    R = 25560;
    mu = 5794000;
    r = 2872e6;
elseif strcmp(Planet,'Neptune')
    R = 24760;
    mu = 6835100;
    r = 4495e6;
elseif strcmp(Planet,'Pluto')  
    R = 1195;
    mu = 830;
    r = 5870e6;
end

T = 2*pi/sqrt(mu_Sun)*r^1.5; %[s]
n = 2*pi/T;                  %[rad/s]
v = sqrt(mu_Sun/r);          %[km/s]
r_SOI = r*(mu/mu_Sun)^(2/5); %[km]