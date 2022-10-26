clc; clear; close all;
addpath('./functions/');

FILE = "../export_data.h5";

%% Get Geometric Params
% radius, x, y, bc_phi
raw_data = h5read("../export_data.h5", "/geometry/inner/data");
r_inner = raw_data(1);
bc_phi_inner = raw_data(4);
bc_inner = h5read("../export_data.h5", "/geometry/inner/bc");

raw_data = h5read("../export_data.h5", "/geometry/outer/data");
r_outer = raw_data(1);
bc_phi_outer = raw_data(4);
bc_outer = h5read("../export_data.h5", "/geometry/outer/bc");

%% Get Solution to Analytical Eqns
[A, B] = analytical_solver(r_inner, bc_phi_inner, bc_inner, ...
    r_outer, bc_phi_outer, bc_outer);

%% Get Steady State Solution
fprintf("Loading Steady State Solutions...")
mesh_0 = load_steady_state_solution(FILE, "mesh_0", 0);
mesh_1 = load_steady_state_solution(FILE, "mesh_1", 0);
mesh_2 = load_steady_state_solution(FILE, "mesh_2", 0);
mesh_3 = load_steady_state_solution(FILE, "mesh_3", 0);
mesh_4 = load_steady_state_solution(FILE, "mesh_4", 0);
mesh_5 = load_steady_state_solution(FILE, "mesh_5", 0);
mesh_6 = load_steady_state_solution(FILE, "mesh_6", 0);
fprintf(" OK\n")

%% Richardson method
fprintf("Generating Richardson Extrapolated Solutions...")
r_0 = richardson_extrapolation(mesh_0, mesh_1);
r_1 = richardson_extrapolation(mesh_1, mesh_2);
r_2 = richardson_extrapolation(mesh_2, mesh_3);
fprintf(" OK\n")

%% Analytical Meshes
fprintf("Generating Analytical Solutions...")
analytical_mesh_0 = analytical_mesh(mesh_0, A, B, r_inner, r_outer, 0);
analytical_mesh_1 = analytical_mesh(mesh_1, A, B, r_inner, r_outer, 0);
analytical_mesh_2 = analytical_mesh(mesh_2, A, B, r_inner, r_outer, 0);
analytical_mesh_3 = analytical_mesh(mesh_3, A, B, r_inner, r_outer, 0);
analytical_mesh_4 = analytical_mesh(mesh_4, A, B, r_inner, r_outer, 0);
analytical_mesh_5 = analytical_mesh(mesh_5, A, B, r_inner, r_outer, 0);
analytical_mesh_6 = analytical_mesh(mesh_6, A, B, r_inner, r_outer, 0);
fprintf(" OK\n")

%% Error Meshes
fprintf("Calculating Error Meshes...")
error_0 = abs(mesh_0 - analytical_mesh_0);
error_1 = abs(mesh_1 - analytical_mesh_1);
error_2 = abs(mesh_2 - analytical_mesh_2);
error_3 = abs(mesh_3 - analytical_mesh_3);
error_4 = abs(mesh_4 - analytical_mesh_4);
error_5 = abs(mesh_5 - analytical_mesh_5);
error_6 = abs(mesh_6 - analytical_mesh_6);

error_0_r = abs(r_0 - analytical_mesh_0);
error_1_r = abs(r_1 - analytical_mesh_1);
error_2_r = abs(r_2 - analytical_mesh_2);
fprintf(" OK\n")

%% Error 2-Norm
fprintf("Calculating Error 2-Norm...")
error_0_norm = calculate_two_norm(error_0);
error_test = calculate_two_norm(error_0);

error_1_norm = calculate_two_norm(error_1);
error_2_norm = calculate_two_norm(error_2);
error_3_norm = calculate_two_norm(error_3);
error_4_norm = calculate_two_norm(error_4);
error_5_norm = calculate_two_norm(error_5);
error_6_norm = calculate_two_norm(error_6);

error_0_norm_r = calculate_two_norm(error_0_r);
error_1_norm_r = calculate_two_norm(error_1_r);
error_2_norm_r = calculate_two_norm(error_2_r);
fprintf(" OK\n")

h = [];
E = [];

h(1) = 1/11;
h(2) = 1/22;
h(3) = 1/44;
h(4) = 1/88;
h(5) = 1/176;
h(6) = 1/352;
h(7) = 1/704;

E(1) = error_0_norm;
E(2) = error_1_norm;
E(3) = error_2_norm;
E(4) = error_3_norm;
E(5) = error_4_norm;
E(6) = error_5_norm;
E(7) = error_6_norm;

E_r = [];
E_r(1) = error_0_norm_r;
E_r(2) = error_1_norm_r;
E_r(3) = error_2_norm_r;

first_order = 1.0*h.^(1.0);
second_order = 1.0*h.^(2.0);
fourth_order = 1.0*h.^(4.0);

hold on
loglog(h, E);
loglog(h(1:end-1), E_r);
loglog(h, first_order, '-k');
loglog(h, second_order, '--k');
loglog(h, fourth_order, '-.k');


p = polyfit(log(h),log(E),0);
% z = polyval(p,log(h));
%loglog(h,exp(z));
hold off
set(gca, 'XScale', 'log', 'YScale', 'log');

xlabel('h','interpreter', 'latex', 'FontSize', 24) 
ylabel('$|\!|\mathrm{T-T_{ex}}|\!|_2$', 'interpreter', 'latex', 'FontSize', 24);
% Create legend
legend1 = legend('Solution', '1. order','2. order', '4. order', 'FontSize', 24);
set(legend1,...
    'Position',[0.630952386221006 0.142460322001624 0.267857137588518 0.165476185934884]);

%% funcs

function s = plot_mesh_surface(mesh)
    [X, Y] = meshgrid_from_mesh(mesh);
    s = surf(X,Y,mesh,FaceColor='interp');
end

function err = calculate_two_norm(error_mesh)
    tmp = sum(error_mesh.^2, 'all') * (1/(size(error_mesh, 1) - 1))^2;
    err = tmp^(1/2);
end

function richardson_mesh = richardson_extrapolation(coarse_mesh, fine_mesh)
    richardson_mesh = fine_mesh(1:2:end, 1:2:end) + (fine_mesh(1:2:end, 1:2:end) - coarse_mesh)/3.0;
end