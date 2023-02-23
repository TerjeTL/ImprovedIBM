clc; clear; close all;
addpath('./functions/');

%FILE = "../seven_levels_ss.h5";
%FILE = "../export_data_high_res.h5";
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

mesh_0 = load_steady_state_solution(FILE, "mesh_2", nan);
mesh_1 = load_steady_state_solution(FILE, "mesh_3", nan);
%mesh_1 = mesh_1(1:8:end, 1:8:end);

r_1 = richardson_extrapolation(mesh_0, mesh_1);

%% Analytical Meshes
analytical_mesh_0 = analytical_mesh(mesh_0, A, B, r_inner, r_outer, nan);
analytical_mesh_1 = analytical_mesh(mesh_1, A, B, r_inner, r_outer, nan);

%% Error Meshes
error_0 = mesh_0 - analytical_mesh_0;
error_1 = mesh_1 - analytical_mesh_1;
error_r = r_1 - analytical_mesh_1(1:2:end, 1:2:end);

%% Plots
tiledlayout(2,3)

% level 0
nexttile
analytical_plot_0 = plot_mesh_surface(analytical_mesh_0);
zlim([0 2]);

nexttile
plot_0 = plot_mesh_surface(mesh_0);
zlim([0 2]);

nexttile
error_plot_0 = plot_mesh_surface(error_1);

% level 1
nexttile
analytical_plot_1 = plot_mesh_surface(analytical_mesh_1);
zlim([0 2]);

nexttile
plot_1 = plot_mesh_surface(mesh_1);
zlim([0 2]);

nexttile
error_plot_1 = plot_mesh_surface(error_r);
%zlim([0 0.45]);

set(gcf,'position',[get(0, 'Screensize')])

%% funcs
function s = plot_mesh_surface(mesh)
    [X, Y] = meshgrid_from_mesh(mesh);
    s = surf(X,Y,mesh,FaceColor='interp');
end

function richardson_mesh = richardson_extrapolation(coarse_mesh, fine_mesh)
    richardson_mesh = fine_mesh(1:2:end, 1:2:end) + (fine_mesh(1:2:end, 1:2:end) - coarse_mesh)/3.0;
    %fine_mesh = fine_mesh(1:2:end, 1:2:end);
    %richardson_mesh = fine_mesh + (fine_mesh - coarse_mesh)./(2.0^(2.0) - 1);
    %richardson_mesh = coarse_mesh + (fine_mesh - coarse_mesh).*2.0^(3.0)./(2.0^(3.0) - 1);
    %richardson_mesh = 4/3 * fine_mesh - 1/3 * coarse_mesh;
end

% function richardson_mesh = richardson_extrapolation(coarse_mesh, fine_mesh)
%     % Lets try the "completed richardson extrapolation"
%     richardson_mesh = fine_mesh;
% 
%     c_0 = 1.0/3.0 * (fine_mesh(1:2:end, 1:2:end) - coarse_mesh);
%     c_2 = 1.0/3.0 * (fine_mesh(3:2:end, 3:2:end) - coarse_mesh(2:end, 2:end));
%     c_1 = 0.5 * (c_0(2:end, 2:end) + c_2);
% 
%     richardson_mesh(1:2:end, 1:2:end) = fine_mesh(1:2:end, 1:2:end) + c_0;
%     richardson_mesh(2:2:end, 2:2:end) = fine_mesh(2:2:end, 2:2:end) + c_1;
% end