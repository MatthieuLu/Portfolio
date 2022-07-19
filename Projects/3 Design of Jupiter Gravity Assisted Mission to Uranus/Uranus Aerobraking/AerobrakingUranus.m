% Assumption: One burn that needs to be done at the intended apoapsis radius.

mu_U = 5794000;
rP = 200420; %% Intended Periapsis Radius
rA = 3*rP; %% Intended Apoapsis Radius
a = 0.5*(rP + rA);
e = (rA - rP)/(rP + rA);
Energy_ell = -1*mu_U/(2*a);
vP = sqrt(2*(Energy_ell+(mu_U/rP)));
vA = sqrt(2*(Energy_ell+(mu_U/rA)));

rP_aero = 6000+25559; %Periapsis Radius of Aerobraking Orbit
a_aero = 0.5*(rA+rP_aero);
e_aero = (rA - rP_aero)/(rA+rP_aero);
Energy_aero = -mu_U/(2*a_aero);
vA_aero = sqrt(2*(Energy_aero+(mu_U/rA)));
delta_v_aero = abs(vA - vA_aero)
