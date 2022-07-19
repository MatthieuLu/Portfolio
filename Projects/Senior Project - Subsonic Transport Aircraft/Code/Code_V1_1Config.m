%% Start
close all
clear
clc

%% Prompts

prompt_0 = "Is this unique (1), sample (2), or quiz (3) configuration? ";
type_inputs = input(prompt_0);

prompt_1 = "Would you like to draw graphs and fit curves? 1 for no, 2 for yes: ";
type_graph = input(prompt_1);

prompt_2 = "Please input 1 for conventional wing, 2 for supercritical: ";
type_airfoil = input(prompt_2);

prompt_3 = "Please enter what type of engine you are using, 1 for JT9D Advanced, 2 for JT9D, and 3 for JT8D: ";
type_engine = input(prompt_3);

prompt_4 = "Please enter how many engines the aircraft has (2, 3, or 4): ";
type_engineqty = input(prompt_4);

if type_engineqty == 1
	type_engineqty = 2;
end

prompt_5 = "Are the engines mounted on the wing (1) or fuselage (2)? ";
type_mount = input(prompt_5);

prompt_6 = "Is the structure conventional (1), hybrid (2), or composites (3)? ";
type_structure = input(prompt_6);

prompt_7 = "Is this flight domestic (1) or international (2)? ";
type_flight = input(prompt_7);

if type_mount == 1
    K_w = 1;
    K_ts = 0.17;
elseif type_mount ==2
    K_w = 1.03;
    K_ts = 0.25;
end

%% Fixed Variables

if type_inputs == 1 %Unique
	AR = 12
	Sweep = 36
	Taper = 0.35;

	height = 35000; %ft
	M = 0.82;
	Range_scheduled = 4000;
	FF_max = 1;
	V_approach = 135;
	
	TOFL = 6000;
	
	PAX = 200; %passengers
	Abreast = 6; %number of passengers abreast
	Aisles = 1;
	Weight_cargo = 4000; %weight of cargo
	Number_flightcrew = 2; %number of flight crew
	Number_stewards = round((PAX/50)+0.49); %number of stewards
elseif type_inputs == 2 %Sample
	AR = 8
	Sweep = 35
	Taper = 0.35;

	height = 35000; %ft
	M = 0.82;
	Range_scheduled = 6000;
	FF_max = .75;
	V_approach = 140;
	TOFL = 9000;
	
	PAX = 275; %passengers
	Abreast = 8; %number of passengers abreast
	Aisles = 2;
	Weight_cargo = 12000; %weight of cargo
	Number_flightcrew = 2; %number of flight crew
	Number_stewards = round((PAX/50)+0.49); %number of stewards
	
	type_engineqty = 3;
elseif type_inputs ==3 %Quiz
	AR = 12.5
	Sweep = 36
	Taper = 0.35;

	height = 35000; %ft
	M = 0.82;
	Range_scheduled = 4000;
	FF_max = 1;
	V_approach = 135;
	TOFL = 6000;
	
	PAX = 200; %passengers
	Abreast = 6; %number of passengers abreast
	Aisles = 2;
	Weight_cargo = 4000; %weight of cargo
	Number_flightcrew = 2; %number of flight crew
	Number_stewards = round((PAX/50)+0.49); %number of stewards
end

gamma = 1.4;
R = 1718;
T = 218.9; %K
rho = .023751; %lb/ft^3
sigma = .953; %Ratio between local and sea level densities
delta = .2360; %Ratio between local and sea level pressures

%a = sqrt(gamma*R*T);
a = 576.4;

%% Digitized Plots

Table = readtable('Tables.xlsx');
Table = Table{:,:};

if type_airfoil == 1
% Figure 1a: Effect of Wing Sweep and Thickness on Mdiv (Conventional)
	x1 = Table(1:165,1);
	y1 = Table(1:165,3);
	z1 = Table(1:165,2);
	[fit_TC,info_TC] = fit([x1,y1],z1,'poly23'); 

% Figure 1b: Effect of Wing Sweep and Thickness on Mdiv (Supercritical)
elseif type_airfoil == 2
	x2 = Table(1:211,4);
	y2 = Table(1:211,6);
	z2 = Table(1:211,5);
	[fit_TC,info_TC] = fit([x2,y2],z2,'poly35');
end
% Figure 2: Wing Drag Divergence Change Due to Lift
if type_airfoil == 1
	x3 = Table(1:25,7);
	y3 = Table(1:25,8);
	[fit_MDIVdelta,info_MDIVdelta] = fit(y3,x3,'power2');
elseif type_airfoil == 2
	x4 = Table(1:29,9);
	y4 = Table(1:29,10);
	[fit_MDIVdelta,info_MDIVdelta] = fit(y4,x4,'power2');
end

% Figure 3: CLmax Conceptual Design Estimates
x5 = Table(1:34,11);
y5 = Table(1:34,12);
[fit_CLmax_takeoff,info_CLmax_takeoff] = fit(x5,y5,'smoothingspline');
x6 = Table(1:31,13);
y6 = Table(1:31,14);
[fit_CLmax_landing,info_CLmax_landing] = fit(x6,y6,'smoothingspline');

% Figure 4: JT8D-9 Fuel Fraction Estimate for All Out Range
x7 = Table(1:41,15);
y7 = Table(1:41,16);
[fit_FF,info_FF] = fit(x7,y7,'cubicinterp'); % FF is Fuel Fraction or W_f/W_to

% Figure 5: Jet Aircraft TOFL with Engine Failure
x8 = Table(1:58,17);
y8 = Table(1:58,18);
[fit_k_2ENG,info_k_2ENG] = fit(y8,x8,'poly2');
x9 = Table(1:63,19);
y9 = Table(1:63,20);
[fit_k_3ENG,info_k_3ENG] = fit(y9,x9,'poly2');
x10 = Table(1:67,21);
y10 = Table(1:67,22);
[fit_k_4ENG,info_k_4ENG] = fit(y10,x10,'poly2');

% Cf
x11 = Table(1:50,23);
y11 = Table(1:50,24);
y11 = y11.*(1e-3);
[fit_Cf, info_Cf] = fit(x11,y11,'power2');

% Body Form Factor
x12 = Table(1:12,25);
y12 = Table(1:12,26);
[fit_k, info_k] = fit(x12,y12,'exp2');

% Max Climb Thrust @ 15000ft
x13 = Table(1:34,27);
y13 = Table(1:34,28);
[fit_T_avbl15, info_T_avbl15] = fit(x13,y13,'smoothingspline');
x14 = Table(1:9,29);
y14 = Table(1:9,30);
[fit_SFC_15, info_SFC_15] = fit(x14,y14,'smoothingspline');

% Max Climb Thrust @ 25000ft
x15 = Table(1:30,31);
y15 = Table(1:30,32);
[fit_T_avbl25, info_T_avbl25] = fit(x15,y15,'smoothingspline');
x = Table(1:6,33);
y = Table(1:6,34);
[fit_SFC_25, info_SFC_25] = fit(x,y,'smoothingspline');

% Max Climb and Thrust @ 35000ft
x16 = Table(1:9,37);
y16 = Table(1:9,38);
[fit_SFC_35, info_SFC_35] = fit(x16,y16,'lin');

% Profile Drag at Takeoff
x17 = Table(1:60,39);
y17 = Table(1:60,40);
[fit_CDp_takeoff, info_CDp_takeoff] = fit(y17,x17,'poly5');

% Profile Drag at Landing
x18 = Table(1:60,41);
y18 = Table(1:60,42);
[fit_CDp_landing, info_CDp_landing] = fit(y18,x18,'poly5');

% Clean Wing CLmax
x19 = Table(1:69,43);
y19 = Table(1:69,45);
z19 = Table(1:69,44);
[fit_CLmax_clean, info_CLmax_clean] = fit([x19,y19],z19,'poly32');

% JT9D Dry Takeoff
x20 = Table(1:18,46);
y20 = Table(1:18,47);
[fit_JT9D_takeoffthrust, info_JT9D_takeoffthrust] = fit(x20,y20,'lin');

% JT9D Climb 
x21 = Table(1:26,48);
y21 = Table(1:26,49);
[fit_JT9D_climbthrust, info_JT9D_climbthrust] = fit(x21,y21,'lin');

%% Figures

if type_graph == 2

    figure
    plot(fit_TC_conv,[x1,y1],z1);
	
	figure
    plot(fit_TC_sc,[x2,y2],z2);
    hold on
    ylim([0 40]);
    zlim([0.07 .2]);
    hold off
	
	figure
    hold on
    plot(x3,y3)
    plot(x4,y4)
    xlim([-.04 .02]) 
    ylim([.4 .7])
    % plot(fit_MDIVdelta_conv)
    % plot(fit_MDIVdelta_sc)
    % xlim([-1 2]) 
    hold off
	
	figure
    hold on
    plot(x5,y5);
    plot(x6,y6);
    hold off
	
	figure
    plot(fit_FF);
	
	figure
    hold on
    plot(x8,y8,x9,y9,x10,y10)
    hold off
	
	figure
    plot(fit_Cf,x11,y11)
    hold on
    set(gca,'xscale','log')
    title('Skin Friction Coefficient wrt Reynolds Number')
    xlabel('Reynolds Number (Re)')
    ylabel('Coefficient of Skin Friction (k)')
    hold off
	
	figure
	plot(fit_k,x12,y12);
    hold on
    title('Body Form Factor wrt Fineness Ratio')
    xlabel('Fineness Ratio (l/d)')
    ylabel('Body Form Factor (k)')
    hold off
	
	figure
	plot(x17,y17,x18,y18)
	
	figure
	hold on
%  	plot([z19,x19],y19)
	plot(fit_CLmax_clean)
	xlim([.04 .16]);
	ylim([0 35]);
	zlim([.8 1.3]);
	hold off
end

%% Placeholders

SFC_jt9d = 0.61;
if type_engine == 1
	SFC_jt9d = 0.9 * SFC_jt9d;
end

%% Configuration Loop

Weightincrement = 0;
Thrust_reqtop_JT9D = 10001;
while Thrust_reqtop_JT9D > 10000
%% Range Loop
	Range = 0;
	Range_ao = 500000;
	Fuelincrement = 0;
	while abs(Range - Range_ao) > 2
%% CL Loop
		CL_guess = 0.5;
		CL_ic = 0;
		while abs(CL_ic - CL_guess) > 0.001
			Mdiv_delta = fit_MDIVdelta(CL_guess);
			Mdiv = (M + .004) - Mdiv_delta;
			TC = fit_TC(Mdiv,Sweep);
			CLmax_landing = fit_CLmax_landing((((cosd(Sweep))^2)*((TC)^2)*AR));
			CLmax_takeoff = fit_CLmax_takeoff((((cosd(Sweep))^2)*((TC)^2)*AR)); %5
			WS_landing = ((V_approach/1.3)^2) * (sigma * CLmax_landing) / 296; %6
			V_cruise = M * a;
			Range_ao = Range_scheduled + 200 + (V_cruise*(3/4)); %7
			FF_jt8d = fit_FF(Range_ao); %8
			if type_inputs == 1
				FF_jt9d = (FF_jt8d * (SFC_jt9d/.78)); %9
			elseif type_inputs == 2
				FF_jt9d = (FF_jt8d * (SFC_jt9d/.78)); %9
			elseif type_inputs == 3
				FF_jt9d = (FF_jt8d * (SFC_jt9d/.78)); %9
			end
			if type_engine == 1 || type_engine == 2
				FF = FF_jt9d + Fuelincrement;
			elseif type_engine == 3
				FF = FF_jt8d + Fuelincrement;
			end
			WS_takeoff = WS_landing/(1-(FF_max*FF));
			WS_ic = .965 * WS_takeoff;
			CL_ic = WS_ic / (1481 * delta * M^2);
			if abs(CL_ic - CL_guess) > 0.001
				CL_guess = (CL_guess + CL_ic)/2;
			end
		end
%% TOFL
		if type_engineqty == 2
			k_1 = fit_k_2ENG((TOFL/1000));
		elseif type_engineqty == 3
			k_1 = fit_k_3ENG((TOFL/1000));
		elseif type_engineqty == 4
			k_1 = fit_k_4ENG((TOFL/1000));
		end
		WT_7V_lo = (k_1/WS_takeoff) * sigma * CLmax_takeoff;
		V_lo = 1.2 *sqrt((296*WS_takeoff)/(sigma*CLmax_takeoff));
		M_lo = V_lo/(661/sqrt(sigma));
		Thrust_sls = fit_JT9D_takeoffthrust(0);
		Thrust_M = fit_JT9D_takeoffthrust(M_lo*0.7);
		WT = (WT_7V_lo * Thrust_M/Thrust_sls) + Weightincrement;
%% Weight
		eta = 1.5*2.5;
		% Wing
		TC_bar = TC + 0.03;
		W_w = (.00945 * (AR^0.8) * ((1+Taper)^0.25) * K_w * (eta^0.5))/((TC_bar^0.4) * cosd(Sweep) * (WS_takeoff^0.695)); %times W_takeoff^1.195
		% Fuselage
		K_f = 11.5;
		L_f = (3.76 * (PAX/Abreast) + 33.2);
		D_f = ((1.75*Abreast) + (1.58*Aisles) + 1);
		if type_flight == 2
			L_f = 1.1 * L_f;
			D_f = 1.1 * D_f;
		end
		W_fus = 0.6727 * K_f * (L_f^0.6) * (D_f^0.72) * (eta^0.3); %times W_takeoff^0.235
		% Landing Gear
		W_lg = 0.04; %times W_takeoff
		% Nacelle & Pylon
		W_np = (0.0555/WT); %times W_takeoff
		% Tail Surface
		W_ts = K_ts * W_w; %times W_takeoff^1.195
		% Power Plant
		W_pp = 1/(3.58*WT); %times W_takeoff
		if type_engine == 1
			W_pp = W_pp * 1.1;
		end
		% Fuel
		W_f = 1.0275 * (FF);
		% Payload
		W_pl = 215*PAX + Weight_cargo;
		% Fixed Equipment
		W_fe = (132*PAX) + (300*type_engineqty) + (260*Number_flightcrew) + (170*Number_stewards);
		W_fe_wto = (0.035); %times W_takeoff
		if type_structure == 2
			W_w = 0.7 * W_w;
			W_ts = 0.7 * W_ts;
			W_np = 0.8 * W_np;
		elseif type_structure == 3
			W_w = 0.7 * W_w;
			W_ts = 0.7 * W_ts;
			W_fus = 0.85 * W_fus;
			W_fe = 0.9 * W_fe;
			W_fe_wto = 0.9 * W_fe_wto;
			W_np = 0.8 * W_np;
		end
		a_W = W_w + W_ts;
		b_W = W_fus;
		c_W = W_lg + W_np + W_pp + W_f + W_fe_wto - 1;
		d_W = W_pl + W_fe;
		W_takeoff = 300000;
		difference_1 = W_takeoff;
		while abs(difference_1) > 0.1
			if abs(difference_1) > 0.5
				difference_1 = a_W*(W_takeoff^1.195) + b_W*(W_takeoff^0.235) + c_W*W_takeoff + d_W;
				W_takeoff = W_takeoff + 4*(difference_1);
			elseif abs(difference_1) <= 0.5
				difference_1 = a_W*(W_takeoff^1.195) + b_W*(W_takeoff^0.235) + c_W*W_takeoff + d_W;
				W_takeoff = W_takeoff + (difference_1);
			end
		end
		W_fuel = W_f * W_takeoff;
%% Wing Sizing
		S = W_takeoff/WS_takeoff;
		b = sqrt(AR*S);
		c_bar = S/b;
		Thrust_total = W_takeoff/WT;
		Thrust_engine = Thrust_total/type_engineqty;
%% Drag
		Re_l = 2.852e6 * 0.5; % 0.5 * rho * v / mu
		% Wing
		Re_w = Re_l * c_bar;
		Cf_w = fit_Cf(Re_w);
		Swet_w = 2 * S * 1.02;
		z_w = ((2-M^2)*cosd(Sweep))/sqrt(1-((M^2).*(cosd(Sweep))));
		k_w = 1+(z_w*TC)+100*(TC^4);
		f_w = k_w*Cf_w*Swet_w;
		% Fuselage
		Swet_fus = 0.9*pi*D_f*L_f;
		Re_fus = Re_l * L_f;
		Cf_fus = fit_Cf(Re_fus);
		fineness_fus = L_f/D_f;
		k_fus = fit_k(fineness_fus);
		f_fus = k_fus*Cf_fus*Swet_fus;
		% Tail
		f_t = 0.35*f_w;
		% Nacelles
		Re_n = Re_l * c_bar;
		Swet_n = 2.1 * sqrt(Thrust_engine) * type_engineqty;
		Cf_n = fit_Cf(Re_n);
		k_n = 1.25;
		f_n = k_n*Cf_n*Swet_n;
		% Pylons
		f_p = 0.2*f_n;
		%Total
		f = f_w + f_fus + f_t + f_n + f_p;
		CD_overall = f/S;
		e = 1/(1.035 + (.38*CD_overall*pi*AR));
		CD_gear = CD_overall;
%% Climb
		W_avgcl = (1.965) * W_takeoff / 2;
		h_avgcl = (20/35) * height;
		sigma_cl = .5702;
		V_ldmax = (12.9/(f*e)^(1/4)) * sqrt(W_avgcl/(sigma_cl*b));
		V_cl = 1.3 * V_ldmax;
		M_cl = V_cl/(973.1/1.68781);
		Thrust_reqcl = (sigma_cl*f*(V_cl^2)/296) + ((94.1/(sigma_cl*e))*((W_avgcl/b)^2)*(1/(V_cl^2)));
		Thrust_cl15 = fit_T_avbl15(M_cl) * 1000;
		Thrust_cl25 = fit_T_avbl25(M_cl) * 1000;
		SFC_15 = fit_SFC_15(M_cl);
		SFC_25 = fit_SFC_25(M_cl);
		Thrust_cl = (Thrust_cl15 + Thrust_cl25)/2;
		SFC_cl = (SFC_15 + SFC_25)/2;
		if type_engine == 1
			SFC_cl = SFC_cl * 0.9;
		end
		Thrust_avbl = (Thrust_engine/Thrust_sls) * Thrust_cl;
		RC = (101*((type_engineqty*Thrust_avbl)-Thrust_reqcl)*V_cl)/W_avgcl;
		Time_cl = height/RC;
		Range_cl = V_cl * (Time_cl/60);
		W_fcl = (type_engineqty*Thrust_avbl) * SFC_cl * (Time_cl/60);
%% Range
		W_0 = W_takeoff-W_fcl;
		W_I = (1-FF) * W_takeoff;
		CL_avg = ((W_0+W_I)/(2*S)) / (1481*delta*(M^2));
		CD_i = (CL_avg^2)/(pi*AR*e);
		CD_c = 0.001;
		CD = CD_overall + CD_i + CD_c;
		LD = CL_avg/CD;
		Thrust_reqcr = ((W_0 + W_I)/(2*LD));
		Thrust_reqcrjt9d = Thrust_reqcr * (Thrust_sls/Thrust_engine);
		Thrust_reqcreng = Thrust_reqcrjt9d / type_engineqty;
		SFC_cr = fit_SFC_35(Thrust_reqcreng);
		if type_engine == 1
			SFC_cr = SFC_cr * 0.9;
		end
		Range_cr = (V_cruise/SFC_cr) * LD * log(W_0/W_I);
		Range = Range_cl + Range_cr;
		difference_2 = Range-Range_ao;
		if     difference_2 <  -200 %200 or less
			Fuelincrement = Fuelincrement + 0.01;
		elseif difference_2 >= 200 %200 or more
			Fuelincrement = Fuelincrement - 0.01;
		elseif difference_2 <  -18 && (Range-Range_ao) >= -200 %Between -20 and -200
			Fuelincrement = Fuelincrement + 0.001;
		elseif difference_2 >  18  && (Range-Range_ao) <= 200 %Between 20 and 200
			Fuelincrement = Fuelincrement - 0.001;
		elseif difference_2 <  -3    && (Range-Range_ao) >= -18 %Between 0 and -20
			Fuelincrement = Fuelincrement + 0.0001;
		elseif difference_2 >= 3    && (Range-Range_ao) <= 18 %Between 0 and 20
			Fuelincrement = Fuelincrement - 0.0001;
		elseif difference_2 <  0    && (Range-Range_ao) >= -3 %Between 0 and -20
			Fuelincrement = Fuelincrement + 0.00001;
		elseif difference_2 >= 0    && (Range-Range_ao) <= 3 %Between 0 and 20
			Fuelincrement = Fuelincrement - 0.00001;
		end
	end
%% Top of Climb
	CL_ic_top = (W_0/S)/(1481*delta*M^2);
	CDi_top = (CL_ic_top^2) / (pi*AR*e);
	CD_top = CD_overall + CDi_top + 0.001;
	LD_top = CL_ic_top / CD_top;
	if type_inputs == 2
		LD_top = 15.05; % SAMPLE
	end
	Thrust_reqtop = W_0 / LD_top;
	Thrust_reqtop_engine = Thrust_reqtop / type_engineqty;
	Thrust_reqtop_JT9D = (Thrust_sls/Thrust_engine) * Thrust_reqtop_engine;

%% Climb Gradients
	% 1st Segment
	CL_seg1 = CLmax_takeoff/(1.2^2);
	CL_to_CLmax_to = CL_seg1/CLmax_takeoff;
	CDp_seg1 = fit_CDp_takeoff(CL_to_CLmax_to);
	CD_seg1 = CD_overall + CDp_seg1 + CD_gear + ((CL_seg1^2)/(pi*AR*e));
	LD_seg1 = CL_seg1/CD_seg1;
	Thrust_reqseg1 = W_takeoff/LD_seg1;
	Thrust_avblseg1 = (Thrust_engine/Thrust_sls) * fit_JT9D_takeoffthrust(M_lo);
	Grad_1 = ((type_engineqty-1)*Thrust_avblseg1-Thrust_reqseg1)*100/W_takeoff;
	% 2nd Segment
	CD_seg2 = CD_overall + CDp_seg1 + ((CL_seg1^2)/(pi*AR*e));
	LD_seg2 = CL_seg1/CD_seg2;
	Thrust_reqseg2 = W_takeoff/LD_seg2;
	Grad_2 = ((type_engineqty-1)*Thrust_avblseg1-Thrust_reqseg2)*100/W_takeoff;
	% 3rd Segment
	CLmax_clean = fit_CLmax_clean(TC,Sweep);
	V_seg3 = 1.2*sqrt((296*WS_takeoff)/(.925*CLmax_clean));
	M_seg3 = V_seg3/659;
	CL_seg3 = CLmax_clean/(1.2^2);
	CD_seg3 = CD_overall + ((CL_seg3^2)/(pi*AR*e));
	LD_seg3 = CL_seg3/CD_seg3;
	Thrust_reqseg3 = W_takeoff/LD_seg3;
	Thrust_avblseg3 = (Thrust_engine/Thrust_sls) * fit_JT9D_climbthrust(M_seg3);
	Grad_3 = ((type_engineqty-1)*Thrust_avblseg3-Thrust_reqseg3)*100/W_takeoff;
	% Approach
	CL_approach = CLmax_takeoff/(1.3^2);
	CL_approach_CLmax_landing = CL_approach/CLmax_takeoff;
	CDp_approach = fit_CDp_takeoff(CL_approach_CLmax_landing);
	CD_approach = CD_overall + CDp_approach + (CL_approach^2)/(pi*AR*e);
	LD_approach = CL_approach/CD_approach;
	W_landing = WS_landing * S;
	Thrust_reqapproach = W_landing/LD_approach;
	V_approach_grad = sqrt((296*WS_landing)/(.953*CL_approach));
	M_approach = V_approach_grad/659;
	Thrust_avblapproach = (Thrust_engine/Thrust_sls) * fit_JT9D_climbthrust(M_approach);
	Grad_approach = ((type_engineqty-1)*Thrust_avblapproach-Thrust_reqapproach)*100/W_landing;
	%Landing
	CL_landing = CLmax_landing/1.3^2;
	CDp_landing = fit_CDp_landing(1/1.3^2);
	CD_landing = CD_overall + CDp_landing + CD_gear + (CL_landing^2)/(pi*AR*e);
	LD_landing = CL_landing/CD_landing;
	Thrust_reqlanding = W_landing/LD_landing;
	V_landing = sqrt((296*WS_landing)/(0.953*CL_landing));
	M_landing = V_landing/659;
	Thrust_avbllanding = (Thrust_engine/Thrust_sls) * fit_JT9D_takeoffthrust(M_landing);
	Grad_landing = ((type_engineqty*Thrust_avbllanding) - Thrust_reqlanding) * 100 / W_landing;
	% Climb Requirements (pg. 39)
	if type_engineqty == 2
		Grad_1_min = 0;
		Grad_2_min = 2.4;
		Grad_3_min = 1.2;
		Grad_approach_min = 2.1;
		Grad_landing_min = 3.2;
	elseif type_engineqty == 3
		Grad_1_min = 0.3;
		Grad_2_min = 2.7;
		Grad_3_min = 1.5;
		Grad_approach_min = 2.4;
		Grad_landing_min = 3.2;
	elseif type_engineqty == 4
		Grad_1_min = 0.5;
		Grad_2_min = 3.0;
		Grad_3_min = 1.7;
		Grad_approach_min = 2.7;
		Grad_landing_min = 3.2;
	end

	if Thrust_reqtop_JT9D > 10000 || Grad_1 < Grad_1_min || Grad_2 < Grad_2_min || Grad_3 < Grad_3_min || Grad_approach < Grad_approach_min || Grad_landing < Grad_landing_min
		Weightincrement = Weightincrement - 0.1;
		disp('Gradient Failed')
	end
end

%% DOC
D = Range_scheduled * 1.15;
Time_gm = (D / (11.866+0.040669*D))/60;
Time_cl_hr = Time_cl/60;
Time_d = 0;
Time_am = 0.10;
if D <= 1400
	k_a = (7+0.015*D);
elseif D > 1400
	k_a = 0.02*D;
end
D_d = 0;
Time_cr = ((D+k_a+20)-(Range_cl*1.15+D_d))/(V_cruise*1.15);
V_B = D / (Time_gm + Time_cl_hr + Time_d + Time_cr + Time_am);
% Block Time
Time_b = Time_gm + Time_cl_hr + Time_d + Time_cr + Time_am;
% Block Fuel
F_cr_am = Thrust_reqcr * SFC_cr * (Time_cr+Time_am);
F_cl = W_fcl;
F_b = F_cr_am + F_cl;
% Flying Operations Cost
% a. Flight Crew
P = W_pl / 2000;
if Number_flightcrew == 2
	CostPerBlockHr = 17.849 * (V_cruise*W_takeoff/1e5)^0.3 + 40.83;
else
	CostPerBlockHr = 24.261 * (V_cruise*W_takeoff/1e5)^0.3 + 57.62;
end
C_TM_FlightCrew = CostPerBlockHr / (V_B * P);
% b. Fuel and Oil
if type_flight == 1
	C_ft = 0.0438;
elseif type_flight == 2
	C_ft = 0.0625;
end
C_ot = 2.15;
C_TM_FuelOil = ( 1.02*F_b*C_ft + type_engineqty*C_ot*Time_b*0.135 ) / (D*P);
% c. Hull Insurance Costs
W_a = W_takeoff * (1-FF) - W_pl - (W_pp*W_takeoff);
C_a = 2.4e6 + (87.5*W_a);
C_e = 590000 + (16*Thrust_engine);
if type_inputs == 2
	C_e = 1314942;
end
C_T = (type_engineqty*C_e) + C_a;
IR_A = 0.01;
U = 630 + 4000/( 1 + 1/(Time_b + 0.5)); % block hrs/yr
C_TM_HullInsurance = IR_A * C_T / (U * V_B * P);
% Direct Maintenance
% a. Airframe-Labor
K_FHa = 4.9169 * log10(W_a/1000) - 6.425;
K_FCa = 0.21256 * log10(W_a/1000)^3.7375; 
T_F = Time_b - Time_gm;
R_L = 8.6;
C_TM_a_AirframeLabor = (K_FHa * T_F + K_FCa) * R_L * (1 + 0.29*(M-1))^1.5 / (V_B * Time_b * P);
% b. Airframe Material
C_FHa = ((1.5994*C_a)/1e6) + 3.4263;
C_FCa = ((1.9229*C_a)/1e6) + 2.2504; 
C_TM_b_AirframeMaterial = (C_FHa*T_F+C_FCa) / (V_B*Time_b*P);
% c. Engine Labor
K_FHe = (type_engineqty*(Thrust_engine/1e3))/(0.82715*(Thrust_engine/1e3)+ 13.639);
K_FCe = 0.2 * type_engineqty;
C_TM_c_EngineLabor = ((K_FHe*T_F+K_FCe)*R_L) / (V_B*Time_b*P);
% d. Engine Material
C_FHe = (((28.2353*C_e)/1e6) - 6.5176) * type_engineqty;
C_FCe = (((3.6698*C_e)/1e6) + 1.3685) * type_engineqty;
C_TM_d_EngineMaterial = (C_FHe*T_F+C_FCe) / (V_B*Time_b*P);
if type_engine == 1
	C_TM_d_EngineMaterial = C_TM_d_EngineMaterial * 1.1;
end
% e. Total Maintenance - Burdened
C_TM_e_Maintenance = (C_TM_a_AirframeLabor + C_TM_b_AirframeMaterial + C_TM_c_EngineLabor + C_TM_d_EngineMaterial)*2;
% Depreciation
D_a = 14;
C_TM_Depreciation = (1/(V_B*P)) * (C_T + 0.06*(C_T-type_engineqty*C_e)+0.3*type_engineqty*C_e) / (D_a*U);
% Total DOC
DOC_perTonMile = C_TM_FlightCrew + C_TM_FuelOil + C_TM_HullInsurance + C_TM_e_Maintenance + C_TM_Depreciation;
DOC_perPAXMile = DOC_perTonMile * P/PAX;

