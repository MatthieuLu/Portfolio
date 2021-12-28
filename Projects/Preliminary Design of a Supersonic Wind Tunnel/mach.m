function M = mach(A)
if A < 1
    error('What''s wrong with you?')
else
    syms M
    AMR = A == ((5 + M^2)^3)/(216*M);
    M = double(vpasolve(AMR,M,[0 Inf]));
end
end