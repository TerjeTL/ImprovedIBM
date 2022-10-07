clc; clear; close all;

% T(r) = A*ln(r) + B

r_1 = 0.15;
r_2 = 0.45;
c_x = 0.5;
c_y = 0.5;

T_1 = 10;
T_2 = 200;

inner_neumann = true;

syms A B

if inner_neumann
    eqns = [A/r_1 == T_1, A*log(r_2) + B == T_2];
else
    eqns = [A*log(r_1) + B == T_1, A*log(r_2) + B == T_2];
end

S = vpasolve(eqns, [A B]);

A_sol = S.A;
B_sol = S.B;

phi_data = load("data_export.csv");
% Remove unwanted datapoints
phi_data(phi_data==0) = nan;

h = 1.0/(width(phi_data)-1);

error = zeros(size(phi_data));

for i = 1:size(error,1)
    for j = 1:size(error,2)
        p_x = (i-1)*h - c_x;
        p_y = (j-1)*h - c_y;
        
        r = sqrt(p_x^2 + p_y^2);
        if (isnan(phi_data(i,j)) == false)
            error(i,j) = ( A_sol*log(r) + B_sol );
        end
    end
end
error(error==0) = nan;

hold on;
%surf(phi_data);
surf(abs(error));
hold off