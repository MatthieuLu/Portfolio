clear
clc
% Minimum Mass
syms Density L A B t m_i;
Density = 2810;
L = 0.18;
t = 0.0005198894197323006737573663303929;
B = 0.05;
A = 3*B*t - 2*t^2;
m_i = Density*L*A;

TotalMinimumMass= 4*m_i
