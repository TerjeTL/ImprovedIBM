clc; clear; close all;

% T(r) = A*ln(r) + B

r_1 = 1;
r_2 = 2;

T_1 = 100;
T_2 = 200;

syms A B
eqns = [A*log(r_1) + B == T_1, A*log(r_2) + B == T_2];

S = vpasolve(eqns, [A B]);

A_sol = S.A;
B_sol = S.B;

fprintf("Coefficients:\n-------------\n" + ...
    "\tA: %f\n\tB: %f\n", A_sol, B_sol);


