clc; clear; close all;
addpath('./functions/');

FILE = "../export_data.h5";

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
boundary_value_mesh = h5read(FILE, "/solutions/mesh_1/steady_state/boundary_values");
flags = h5read(FILE, "/solutions/mesh_0/steady_state/boundary_values");

[X, Y] = meshgrid_from_mesh(boundary_value_mesh);
analytical = sqrt((X-0.5).^2 + (Y-0.5).^2);
analytical = analytical_value(double(A),double(B),analytical);
analytical(boundary_value_mesh == 0) = 0;
analytical(boundary_value_mesh > 1.5) = 0;


mesh_1 = load_steady_state_solution(FILE, "mesh_3", nan);


[theta_sol, value_sol, rho_sol] = convert_to_array(boundary_value_mesh);
[theta_an, value_an, rho_an] = convert_to_array(analytical);

error = value_sol - value_an;

yyaxis left
plot(theta_sol, error)
yyaxis right
plot(theta_sol, rho_an-r_inner)

%% Analytical Meshes
analytical_mesh_0 = analytical_mesh(mesh_0, A, B, r_inner, r_outer, nan);
analytical_mesh_1 = analytical_mesh(mesh_1, A, B, r_inner, r_outer, nan);

%% Error Meshes
error_0 = abs(mesh_0 - analytical_mesh_0);
error_1 = abs(mesh_1 - analytical_mesh_1);

%% Plots
tiledlayout(2,3)

% level 0
nexttile
analytical_plot_0 = plot_mesh_surface(analytical_mesh_0);
zlim([0 2.5]);

nexttile
plot_0 = plot_mesh_surface(mesh_0);
zlim([0 2.5]);

nexttile
error_plot_0 = plot_mesh_surface(error_0);

% level 1
nexttile
analytical_plot_1 = plot_mesh_surface(analytical_mesh_1);
zlim([0 2.5]);

nexttile
plot_1 = plot_mesh_surface(mesh_1);
zlim([0 2.5]);

nexttile
error_plot_1 = plot_mesh_surface(error_1);
%zlim([0 0.45]);

set(gcf,'position',[get(0, 'Screensize')])

%% funcs
function s = plot_mesh_surface(mesh)
    [X, Y] = meshgrid_from_mesh(mesh);
    s = surf(X,Y,mesh,FaceColor='interp');
end

function [theta, value, rho] = convert_to_array(boundary_mat)
    [X, Y] = meshgrid_from_mesh(boundary_mat);
    r = sqrt((X-0.5).^2 + (Y-0.5).^2);

    [theta, rho] = cart2pol(X-0.5, Y-0.5);

    theta = theta .* 180/pi;
    
    boundary_mat(boundary_mat > 1.5) = 0;

    mesh = boundary_mat;
    mesh_to_array = mesh(mesh > 0);
    theta_to_array = theta(mesh > 0);
    rho_to_array = rho(mesh > 0);

    [theta, idx] = sort(theta_to_array);
    value = mesh_to_array(idx);
    rho = rho_to_array(idx);
end