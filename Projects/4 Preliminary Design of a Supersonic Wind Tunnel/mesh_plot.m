close all; clear all; clc;
format long e;
option = 2;                   % different options to generate the geometry
dx     = 0.01;               % size of the grid
x_span = 9;                  % total length of the nozzle
nx     = ceil(x_span/dx)+1;
x      = linspace(0,x_span,nx);
y      = geometry(x,option);
ny     = ceil(max(y)/dx)+1;
for i = 1:nx
    ypt = linspace(0,y(i),ny);
%     % laminar grid
%     Yn = 30;
%     for j = (ny-Yn):(ny+2*Yn)
%         Y = ((ny+2*Yn)-j)/3/Yn;
%         ypt(j) = ypt(ny-Yn) + (y(i) - ypt(ny-Yn))*(1 - Y^3);
%     end
    for j = 1:length(ypt)
        pt(i,j,1) = x(i);
        pt(i,j,2) = ypt(j);
    end
end
mesh = fopen('mesh.dat','w');
fprintf(mesh,'TITLE = "Nozzle Design" \n VARIABLES = x, y \n');
fprintf(mesh,'ZONE I = %d, J = %d, DATAPACKING = POINT \n',nx,length(ypt));
fprintf(mesh,'DT = (DOUBLE, DOUBLE) \n');
for j = 1:length(ypt)
    for i = 1:nx
        fprintf(mesh,'%.16e, \t %.16e   \n',[pt(i,j,1),pt(i,j,2)]);
    end
end
fclose(mesh);

function y = geometry(x,option)
switch option
    case 1
        % nozzle
        y = (3/4 - atan(x-tan(1))/4) .*(x>=0 & x<2*tan(1))...
            + 1/2 .*(x>=2*tan(1) & x<6)...
            + (3/4 + atan(x-tan(1)-6)/4) .*(x>=6 & x<=10);
        % compression wedge
        %         y = 8.*(x>=0 & x<4)+(8-(x-4).*tan(pi/9)).*(x>=4 & x<=10);
        % expansion wedge
        %         y = 8.*(x>=0 & x<4)+(8+(x-4).*tan(pi/9)).*(x>=4 & x<=10);
        figure; axis([0 max(x) 0 max(y)*2]); hold on; grid on;
        plot(x,y);
    case 2
        pt = xlsread('geometry');
        y  = interp1(pt(:,1),pt(:,2),x,'spline');
        figure; axis([0 max(x) 0 max(pt(:,2))*2]); hold on; grid on;
        plot(pt(:,1), pt(:,2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'r');
        plot(x,y);
    case 3
        figure; axis([0 max(x) 0 max(x)/3]); hold on; grid on;
        i = 0;
        while i < 100
            i = i + 1;
            k = waitforbuttonpress;
            pt_temp = get(gca, 'currentpoint');
            pt(i,1) = pt_temp(1,1);
            pt(i,2) = pt_temp(1,2);
            plot(pt(i,1), pt(i,2), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'r');
            if ~isempty(get(gcf,'CurrentCharacter'))
                break;
            end
        end
        p     = polyfit(pt(:,1),pt(:,2),7);
        y     = polyval(p,x);
        plot(x,y);
end
end