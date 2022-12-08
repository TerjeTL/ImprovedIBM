clc; clear; close all;
addpath('./functions/');

%FILE = "../seven_levels_ss.h5";
FILE = "../export_data_high_res.h5";

%% Get Geometric Params
% radius, x, y, bc_phi
raw_data = h5read(FILE, "/geometry/inner/data");
r_inner = raw_data(1);
bc_phi_inner = raw_data(4);
bc_inner = h5read(FILE, "/geometry/inner/bc");

raw_data = h5read(FILE, "/geometry/outer/data");
r_outer = raw_data(1);
bc_phi_outer = raw_data(4);
bc_outer = h5read(FILE, "/geometry/outer/bc");

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
mesh_7 = load_steady_state_solution(FILE, "mesh_7", 0);
fprintf(" OK\n")

%% Richardson method
fprintf("Generating Richardson Extrapolated Solutions...")

%mesh_0c = crop_boundaries(mesh_0, r_inner+0.1, r_outer-0.1);
%mesh_1c = crop_boundaries(mesh_1, r_inner+0.1, r_outer-0.1);
%mesh_2c = crop_boundaries(mesh_2, r_inner+0.1, r_outer-0.1);
%mesh_3c = crop_boundaries(mesh_3, r_inner+0.1, r_outer-0.1);
%mesh_4c = crop_boundaries(mesh_4, r_inner+0.1, r_outer-0.1);
%mesh_5c = crop_boundaries(mesh_5, r_inner+0.1, r_outer-0.1);
%mesh_6c = crop_boundaries(mesh_6, r_inner+0.1, r_outer-0.1);

r_0 = richardson_extrapolation(mesh_0, mesh_1);
r_1 = richardson_extrapolation(mesh_1, mesh_2);
r_2 = richardson_extrapolation(mesh_2, mesh_3);
r_3 = richardson_extrapolation(mesh_3, mesh_4);
r_4 = richardson_extrapolation(mesh_4, mesh_5);
r_5 = richardson_extrapolation(mesh_5, mesh_6);
r_6 = richardson_extrapolation(mesh_6, mesh_7);
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
analytical_mesh_7 = analytical_mesh(mesh_7, A, B, r_inner, r_outer, 0);

analytical_mesh_0c = crop_boundaries(analytical_mesh_0, r_inner+0.1, r_outer-0.1);
analytical_mesh_1c = crop_boundaries(analytical_mesh_1, r_inner+0.1, r_outer-0.1);
analytical_mesh_2c = crop_boundaries(analytical_mesh_2, r_inner+0.1, r_outer-0.1);
analytical_mesh_3c = crop_boundaries(analytical_mesh_3, r_inner+0.1, r_outer-0.1);
analytical_mesh_4c = crop_boundaries(analytical_mesh_4, r_inner+0.1, r_outer-0.1);
%analytical_mesh_5c = crop_boundaries(analytical_mesh_5, r_inner+0.1, r_outer-0.1);
%analytical_mesh_6c = crop_boundaries(analytical_mesh_6, r_inner+0.1, r_outer-0.1);
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
error_7 = abs(mesh_7 - analytical_mesh_7);

error_0_coarse = error_0;
error_1_coarse = error_1(1:2:end,1:2:end);
error_2_coarse = error_2(1:4:end,1:4:end);
error_3_coarse = error_3(1:8:end,1:8:end);
error_4_coarse = error_4(1:12:end,1:12:end);
error_5_coarse = error_5(1:24:end,1:24:end);
error_6_coarse = error_6(1:48:end,1:48:end);
error_7_coarse = error_7(1:96:end,1:96:end);

error_0_r = abs(r_0 - analytical_mesh_0);
error_1_r = abs(r_1 - analytical_mesh_1);
error_2_r = abs(r_2 - analytical_mesh_2);
error_3_r = abs(r_3 - analytical_mesh_3);
error_4_r = abs(r_4 - analytical_mesh_4);
error_5_r = abs(r_5 - analytical_mesh_5);
error_6_r = abs(r_6 - analytical_mesh_6);

error_0_r_coarse = error_0_r;
error_1_r_coarse = error_1_r(1:2:end, 1:2:end);
error_2_r_coarse = error_2_r(1:4:end, 1:4:end);
error_3_r_coarse = error_3_r(1:8:end, 1:8:end);
error_4_r_coarse = error_4_r(1:16:end, 1:16:end);
error_5_r_coarse = error_5_r(1:32:end, 1:32:end);
error_6_r_coarse = error_6_r(1:64:end, 1:64:end);

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
error_7_norm = calculate_two_norm(error_7);

error_0_norm_r = calculate_two_norm(error_0_r);
error_1_norm_r = calculate_two_norm(error_1_r);
error_2_norm_r = calculate_two_norm(error_2_r);
error_3_norm_r = calculate_two_norm(error_3_r);
error_4_norm_r = calculate_two_norm(error_4_r);
error_5_norm_r = calculate_two_norm(error_5_r);
error_6_norm_r = calculate_two_norm(error_6_r);

error_0_norm_rp = calculate_two_norm_point(error_0_r_coarse, 4, 9);
error_1_norm_rp = calculate_two_norm_point(error_1_r_coarse, 4, 9);
error_2_norm_rp = calculate_two_norm_point(error_2_r_coarse, 4, 9);
error_3_norm_rp = calculate_two_norm_point(error_3_r_coarse, 4, 9);
error_4_norm_rp = calculate_two_norm_point(error_4_r_coarse, 4, 9);
error_5_norm_rp = calculate_two_norm_point(error_5_r_coarse, 4, 9);
error_6_norm_rp = calculate_two_norm_point(error_6_r_coarse, 4, 9);

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
h(8) = 1/1408;

E(1) = error_0_norm;
E(2) = error_1_norm;
E(3) = error_2_norm;
E(4) = error_3_norm;
E(5) = error_4_norm;
E(6) = error_5_norm;
E(7) = error_6_norm;
E(8) = error_7_norm;

E_r = [];
E_r(1) = error_0_norm_rp;
E_r(2) = error_0_norm_rp;
E_r(3) = error_2_norm_rp;
E_r(4) = error_3_norm_rp;
E_r(5) = error_4_norm_rp;
E_r(6) = error_5_norm_rp;
E_r(7) = error_6_norm_rp;

%E_p = [];
%E_p(1) = error_0_coarse(7,3);
%E_p(2) = error_1_coarse(7,3);
%E_p(3) = error_2_coarse(7,3);
%E_p(4) = error_3_coarse(7,3);
%E_p(5) = error_4_coarse(7,3);
%E_p(6) = error_5_coarse(7,3);
%E_p(7) = error_6_coarse(7,3);
%E_p(8) = error_7_coarse(7,3);

first_order = 1.0*h.^(1.0);
second_order = 1.0*h.^(2.0);
third_order = 1.0*h.^(3.0);
fourth_order = 1.0*h.^(4.0);

hold on
loglog(h, E);
%loglog(h, E_p);
loglog(h(1:end-1), E_r*35.);
loglog(h, first_order, '-k');
loglog(h, second_order, '--k');
loglog(h, third_order, '--k');
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

function err = calculate_two_norm_point(mesh, i, j)
    val = mesh(i,j);
    tmp = val^2 * (1/(size(mesh, 1) - 1))^2;
    err = tmp^(1/2);
end

function richardson_mesh = richardson_extrapolation(coarse_mesh, fine_mesh)
    richardson_mesh = fine_mesh(1:2:end, 1:2:end) + (fine_mesh(1:2:end, 1:2:end) - coarse_mesh)/3.0;
    %fine_mesh = fine_mesh(1:2:end, 1:2:end);
    %richardson_mesh = fine_mesh + (fine_mesh - coarse_mesh)./(2.0^(2.0) - 1);
    %richardson_mesh = coarse_mesh + (fine_mesh - coarse_mesh).*2.0^(3.0)./(2.0^(3.0) - 1);
    %richardson_mesh = 4/3 * fine_mesh - 1/3 * coarse_mesh;
end