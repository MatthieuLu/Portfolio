clear all
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

subplot(2,2,3)
hold on
plot(xtop,Cp)

x = [8 16 64]
for y = [1:3]

    n = x(y)

    for i = [1:n]
        %ith Panel Vertex Angles (ccw)
        angle_i = (360/(2*n)) - (360/n)*(i-1);
        angle_i_1 = (360/(2*n)) - (360/n)*(i);

        %ith Panel Vertex Coordinates
        X_i(i) = xst * cosd(angle_i);
        X_i(i+1) = xst * cosd(angle_i_1);
        Y_i(i) = yst * sind(angle_i);
        Y_i(i+1) = yst * sind(angle_i_1);

        %ith Panel Midpoint
        x_i(i) = (X_i(i+1) + X_i(i))/2;
        y_i(i) = (Y_i(i+1) + Y_i(i))/2;

        %ith Panel Angle(phi) and Normal Angle(beta)
        rise_i = Y_i(i+1) - Y_i(i);
        run_i = X_i(i+1) - X_i(i);
        phi_i(i) = wrapTo360(atan2d(rise_i,run_i));
        beta_i(i) = wrapTo360(phi_i(i) + 90);

        %Placeholder
        p(i,1) =  -U*2*pi*cosd(beta_i(i));

        for j = 1:n
            %jth Panel Vertex Angles (ccw)
            angle_j = (360/(2*n)) - (360/n)*(j-1);
            angle_j_1 = (360/(2*n)) - (360/n)*(j);

            %jth Panel Vertex Coordinates
            X_j(j) = xst * cosd(angle_j);
            X_j(j+1) = xst * cosd(angle_j_1);
            Y_j(j) = yst * sind(angle_j);
            Y_j(j+1) = yst * sind(angle_j_1);

            %jth Panel Angle(phi) and Normal Angle(beta)
            rise_j = Y_j(j+1) - Y_j(j);
            run_j = X_j(j+1) - X_j(j);
            phi_j(j) = wrapTo360(atan2d(rise_j,run_j));
            beta_j(j) = wrapTo360(phi_j(i) + 90);

            %Integration Constants
            A(i,j)  = -(x_i(i) - X_j(j))*cosd(phi_j(j)) - (y_i(i) - Y_j(j))*sind(phi_j(j));
            B(i,j)  = (x_i(i)-X_j(j))^2 + (y_i(i)-Y_j(j))^2;
            C(i,j)  = sind(phi_i(i) - phi_j(j));
            D(i,j)  = (y_i(i)-Y_j(j))*cosd(phi_i(i)) - (x_i(i)-X_j(j))*sind(phi_i(i));
            Sj(i,j) = sqrt((X_j(j+1)-X_j(j))^2 + (Y_j(j+1)-Y_j(j))^2);
            E(i,j)  = sqrt(B(i,j) - A(i,j)^2);

            I(i,j)  = (C(i,j)/2)*log((Sj(i,j)^2+(2*A(i,j)*Sj(i,j))+B(i,j))/B(i,j)) + ((D(i,j)-(A(i,j)*C(i,j)))/E(i,j))*(atan((Sj(i,j)+A(i,j))/E(i,j))-atan(A(i,j)/E(i,j)));        
            if i == j
                I(i,j) = pi;
            end
        end
    end

    %Lambda
    lambda = linsolve(I,p);
    lambdaRatio = lambda./(2*pi*U);



    for i = 1:n
        V(i) = U*sind(beta_i(i));
        for j = 1:n
            if i ~= j
                V(i) = V(i) + (lambda(j)/(2*pi))*((((D(i,j)-A(i,j)*C(i,j))/(2*E(i,j)))*log((Sj(i,j)^2+(2*A(i,j)*Sj(i,j))+B(i,j))/B(i,j))) - C(i,j)*(atan((Sj(i,j)+A(i,j))/E(i,j)) - atan(A(i,j)/E(i,j))));
            end
        end
    end

    Cp_i = 1 - (V/U).^2;
    
    subplot(2,2,3)
    hold on
    plot(x_i,Cp_i,'.')
    title('Cp Distribution')


end

legend('Exact','8','16','64')
hold off

subplot(2,2,1)
hold on
axis equal
fimplicit(ellipse, [-10 10 -10 10],'b');
contour(X,Y,psi,[100])
title('Ellipse and Streamlines for a = 4.2, c = 0.8, and U = 1')
xlabel('x')
ylabel('y')
hold off

subplot(2,2,2)
hold on
axis equal
fimplicit(ellipse, [-10 10 -10 10]);
plot(X_i,Y_i)
plot(x_i,y_i,'k.','MarkerSize',5)
title('Ellipse and Discretization Representation (Panels and Control Points for n = 64)')
xlabel('x')
ylabel('y')
legend('Exact Ellipse','Panels','Control Points')
hold off






