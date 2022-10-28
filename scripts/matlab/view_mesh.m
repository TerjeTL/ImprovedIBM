clc; clear; close all;
addpath('./functions/');

FILE = "../seven_levels_ss.h5";
%FILE = "../export_data.h5";

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

mesh_0 = load_steady_state_solution(FILE, "mesh_1", nan);

mesh_1 = load_steady_state_solution(FILE, "mesh_2", nan);
mesh_1 = mesh_1(1:2:end, 1:2:end);

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