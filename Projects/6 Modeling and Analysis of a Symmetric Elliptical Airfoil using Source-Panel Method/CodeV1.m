clear
clc
close all

U = 1;
a = 4.2;
c = .8;

xst = (a+((c^2)/a));
yst = (a-((c^2)/a));

ellipse = @(x,y) (x./xst).^2 + (y./yst).^2 - 1;

%Streamline

[X1,Y1] = meshgrid(-10:.1:-.01,-10:.1:10);
gamma1 = -atan2((2*X1.*Y1) , ((X1.^2)-(Y1.^2)-(4*(c^2))));
mu1 = (4*(X1.^2).*(Y1.^2)) + ((4*(c^2))-(X1.^2)+(Y1.^2)).^2;
psi1 = (U/(2*(c^2))) * ((((a^2)+(c^2))*Y1) - (((a^2)-(c^2))*(mu1.^(1/4)).*sin(gamma1./2)));

[X2,Y2] = meshgrid(0:.1:10,-10:.1:10);
gamma2 = atan2((2*X2.*Y2) , ((X2.^2)-(Y2.^2)-(4*(c^2))));
mu2 = (4*(X2.^2).*(Y2.^2)) + ((4*(c^2))-(X2.^2)+(Y2.^2)).^2;
psi2 = (U/(2*(c^2)))*((((a^2)+(c^2))*Y2)-(((a^2)-(c^2))*(mu2.^(1/4)).*sin(gamma2./2)));

[X,Y] = meshgrid(-10:.1:10,-10:.1:10);
psi = [psi1 psi2];
gamma = [gamma1 gamma2];
mu = [mu1 mu2];

tiledlayout(1,2)

nexttile
hold on
axis equal
fimplicit(ellipse, [-10 10 -10 10]);
contour(X,Y,psi,[100])
hold off

%Coefficient of Pressure

%Body

xtop1 = linspace(-xst,-.01,xst*100);
xtop2 = linspace(0,xst,(xst*100)+1);
xtop  = linspace(-xst,xst,2*(xst*100)+1);
ytop1 = sqrt(1-(xtop1/xst).^2)*yst;
ytop2 = sqrt(1-(xtop2/xst).^2)*yst;
ytop  = sqrt(1-(xtop/xst).^2)*yst;

gammatop = [-atan2(2*xtop1.*ytop1,((xtop1.^2)-(ytop1.^2)-4*c^2)) atan2(2*xtop2.*ytop2,((xtop2.^2)-(ytop2.^2)-4*c^2))];
mutop = (4*(xtop.^2).*(ytop.^2))+(4*c^2-xtop.^2+ytop.^2).^2;

u = (U./(2*c^2*mutop.^(1/4))).*((a^2+c^2)*mutop.^(1/4)-(a^2-c^2)*(abs(xtop).*cos(gammatop/2)+ytop.*sin(gammatop/2)));
v = U.*(a^2-c^2)*(ytop.*cos(gammatop/2)-abs(xtop).*sin(gammatop/2))./(2*c^2*mutop.^(1/4));
Cp = 1 - (u.^2 + v.^2)/(U^2);

nexttile

plot(xtop,Cp)

n = 8

for i = [1:n]
    %Angle of ith panel vertices relative to center
    angle_i_1 = (360/(2*n)) + (360/n)*(i-1);
    angle_i_2 = (360/(2*n)) + (360/n)*(i);
    angle_i_m = (angle_i_1 + angle_i_2)/2;
    
    %ith panel vertices
    X_i_1 = xst * cosd(angle_i_1);
    Y_i_1 = yst * sind(angle_i_1);
    X_i_2 = xst * cosd(angle_i_2);
    Y_i_2 = yst * sind(angle_i_2);
    x_v_i(i) = X_i_1;   
    x_v_i(i+1) = X_i_2; %x-coordinate vertices of ith panel (Xi is i, Xi+1 is i+1)
    y_v_i(i) = Y_i_1;
    y_v_i(i+1) = Y_i_2; %y-coordinate vertices of ith panel 
    
    
    %ith panel midpoints
    x_i(i) = (X_i_1 + X_i_2)/2;
    y_i(i) = (Y_i_1 + Y_i_2)/2;
    
    %ith panel angle and normal angle
    run_i(i) = X_i_2 - X_i_1;
    rise_i(i) = Y_i_2 - Y_i_1;
    phi_i(i) = wrapTo360(atan2d(rise_i(i),run_i(i)));
    beta_i(i) = wrapTo360(phi_i(i) + 90);
    
    
    for j = [1:n]
        %Angle of jth panel vertices relative to center
        angle_j = (360/(2*n)) + (360/n)*(j-1);
        angle_j_1 = (360/(2*n)) + (360/n)*(j);
    
        %jth panel vertices
        X_j(j) = xst * cosd(angle_j);
        Y_j(j) = yst * sind(angle_j);
        X_j_1(j) = xst * cosd(angle_j_1);
        Y_j_1(j) = yst * sind(angle_j_1);

        %jth panel angle and normal angle
        run_j(j) = X_j_1(j) - X_j(j);
        rise_j(j) = Y_j_1(j) - Y_j(j);
        phi_j(j) = wrapTo360(atan2d(rise_j(j),run_j(j)));
        beta_j(j) = wrapTo360(phi_j(j) + 90);
        
        A(i,j) = -(x_i(i)-X_j(j))*cosd(phi_j(j)) - (y_i(i)-Y_j(j))*sind(phi_j(j));
        B(i,j) = (x_i(i)-X_j(j))^2 + (y_i(i)-Y_j(j))^2;
        C(i,j) = sind(phi_i(i)-phi_j(j));
        D(i,j) = (y_i(i)-Y_j(j))*cosd(phi_i(i)) - (x_i(i)-X_j(j))*sind(phi_i(i));
        E(i,j) = sqrt(B(i,j)-A(i,j)^2);
        Sj(i,j) = sqrt((X_j_1(j)-X_j(j))^2 + (Y_j_1(j)-Y_j(j))^2);
        
        I(i,j) = (C(i,j)/2)*log((Sj(i,j)^2+(2*A(i,j)*Sj(i,j))+B(i,j))/B(i,j)) + ((D(i,j)-(A(i,j)*C(i,j)))/E(i,j))*(atan((Sj(i,j)+A(i,j))/E(i,j))-atan(A(i,j)/E(i,j)));
    end
end

figure
hold on
fimplicit(ellipse, [-10 10 -10 10]);
plot(x_i,y_i)
hold off

    
    






