clear
clc
%Buckling Criterion
syms L t SF B F E I K A Sigma Sigmacr
SF = 25;
B = 0.05;
F = 5668.53;
E = 69e9;
K = 1;
L = .18;
I = (B*(t^3)+2*t*((B-t)^3))/3;
A = 3*B*t - 2*t^2;
Pcrit = ((pi^2)*E*I/((K*L)^2));
Sigma = F/A;
Sigmacr = Pcrit/A

eqn = SF*Sigma == Sigmacr
solutions = vpa(solve(eqn,t))