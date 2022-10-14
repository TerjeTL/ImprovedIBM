clc; clear; close all;

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

mesh_0 = steady_state_solution("mesh_0");
mesh_1 = steady_state_solution("mesh_1");

%% Analytical Meshes
analytical_mesh_0 = analytical_mesh(mesh_0, A, B, r_inner, r_outer);
analytical_mesh_1 = analytical_mesh(mesh_1, A, B, r_inner, r_outer);

%% Error Meshes
error_0 = abs(mesh_0 - analytical_mesh_0);
error_1 = abs(mesh_1 - analytical_mesh_1);

%% Plots
tiledlayout(2,3)

% level 0
nexttile
analytical_plot_0 = plot_mesh_surface(analytical_mesh_0);
zlim([0 250]);

nexttile
plot_0 = plot_mesh_surface(mesh_0);
zlim([0 250]);

nexttile
error_plot_0 = plot_mesh_surface(error_0);

% level 1
nexttile
analytical_plot_1 = plot_mesh_surface(analytical_mesh_1);
zlim([0 250]);

nexttile
plot_1 = plot_mesh_surface(mesh_1);
zlim([0 250]);

nexttile
error_plot_1 = plot_mesh_surface(error_1);
zlim([0 0.45]);

set(gcf,'position',[get(0, 'Screensize')])

%% funcs
function mesh = steady_state_solution(mesh_level_str)
    dir = "/solutions/" + mesh_level_str + "/time_dict";
    time_list = h5read("../export_data.h5", dir);
    str_time_list = compose('%0.6f', time_list);
    dir = "/solutions/" + mesh_level_str + "/time_data/" + str_time_list(end);
    mesh = h5read("../export_data.h5", dir);
    mesh(mesh==0) = nan;
end

function mesh = analytical_mesh(mesh_size, A, B, r_inner, r_outer)
    [X, Y] = meshgrid_from_mesh(mesh_size);
    r = sqrt((X-0.5).^2 + (Y-0.5).^2);

    mesh = zeros(size(r));
    for i = 1:size(mesh,1)
        for j = 1:size(mesh,2)
            if (r(i,j) > r_inner && r(i,j) < r_outer)
                mesh(i,j) = analytical_value(A, B, r(i,j));
            end
        end
    end
    mesh(mesh==0) = nan;
end

function [X, Y] = meshgrid_from_mesh(mesh)
    x = linspace(0, 1, size(mesh, 1));
    y = linspace(0, 1, size(mesh, 2));
    
    [X, Y] = meshgrid(x, y);
end

function s = plot_mesh_surface(mesh)
    [X, Y] = meshgrid_from_mesh(mesh);
    s = surf(X,Y,mesh,FaceColor='interp');
end