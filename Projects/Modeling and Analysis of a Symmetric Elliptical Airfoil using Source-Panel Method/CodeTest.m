clear
clc
syms A B C D E Sj j n Xb Yb a c U mu gamma x y psi u v Cp x_i X_j phi_j y_i Y_j phi_i X_j1 Y_j1 I_ij V_i lambda_j V_inf beta_i

eq1 = latex((Xb/(a+(c^2/a)))^2+(Yb/(a-(c^2/a)))^2 == 1)
eq2 = latex(psi == (U/(2*c^2))*((a^2+c^2)*y - (a^2-c^2)*(mu^(1/4))*sin(gamma/2)))
eq3 = latex(mu == (4*x^2*y^2)+(4*c^2-x^2+y^2)^2)
eq4 = latex(gamma == atan((2*x*y)/(x^2-y^2-(4*c^2))))

eq5 = latex(u == (U/(2*c^2*mu^(1/4)))*((a^2+c^2)*mu^(1/4)-(a^2-c^2)*(abs(x).*cos(gamma/2)+y*sin(gamma/2))))
eq6 = latex(v == (U*(a^2-c^2)*(y*cos(gamma/2)-abs(x)*sin(gamma/2)))/(2*c^2*mu^(1/4)))

eq7 = latex(Cp == 1 - (u.^2 + v.^2)/(U^2))

eq8 = latex(A == -(x_i-X_j)*cos(phi_j)-(y_i-Y_j)*sin(phi_j))
eq9 = latex(B == (x_i-X_j)^2+(y_i-Y_j)^2)
eq10= latex(C == sin(phi_i - phi_j))
eq11= latex(D == (y_i-Y_j)*cos(phi_i)-(x_i-X_j)*sin(phi_i))
eq12= latex(E == sqrt(B-A^2))
eq13= latex(Sj == sqrt((X_j1-X_j)^2+(Y_j1-Y_j)^2))
eq14= latex(I_ij == (((D-A*C)/(2*E))*log((Sj^2+(2*A*Sj)+B)/B)) - C*(atan((Sj+A)/E) - atan(A/E)))

eq15= latex(V_i == V_inf*sin(beta_i) + symsum((lambda_j/(2*pi))*I_ij,j,1,n))
eq16= latex(lambda_j*I_ij == -2*pi*V_inf*cos(beta_i))

