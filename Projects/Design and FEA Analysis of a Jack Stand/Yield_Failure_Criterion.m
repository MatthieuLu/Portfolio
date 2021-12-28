clear
clc
%Yield Criterion
syms sigma sigmastar t B F A SF
sigmastar = 5.515e7
SF = 25
B = 0.05
F = 5668.53
A = 3*B*t - 2*t^2
sigma = F/A

eqn = SF*sigma == sigmastar
S = vpa(solve(eqn,t))