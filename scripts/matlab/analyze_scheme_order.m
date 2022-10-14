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
fprintf("Loading Steady State Solutions ...")
mesh_0 = steady_state_solution("mesh_0", 0);
mesh_1 = steady_state_solution("mesh_1", 0);
mesh_2 = steady_state_solution("mesh_2", 0);
mesh_3 = steady_state_solution("mesh_3", 0);
fprintf(" OK\n")

%% Analytical Meshes
fprintf("Generating Analytical Solutions ...")
analytical_mesh_0 = analytical_mesh(mesh_0, A, B, r_inner, r_outer, 0);
analytical_mesh_1 = analytical_mesh(mesh_1, A, B, r_inner, r_outer, 0);
analytical_mesh_2 = analytical_mesh(mesh_2, A, B, r_inner, r_outer, 0);
analytical_mesh_3 = analytical_mesh(mesh_3, A, B, r_inner, r_outer, 0);
fprintf(" OK\n")

%% Error Meshes
fprintf("Calculating Error ...")
error_0 = abs(mesh_0 - analytical_mesh_0);
error_1 = abs(mesh_1 - analytical_mesh_1);
error_2 = abs(mesh_2 - analytical_mesh_2);
error_3 = abs(mesh_3 - analytical_mesh_3);
fprintf(" OK\n")

%% Error 2-Norm
fprintf("Calculating Error 2-Norm...")
error_0_norm = norm(error_0);
error_test = calculate_two_norm(error_0);

error_1_norm = norm(error_1);
error_2_norm = norm(error_2);
error_3_norm = norm(error_3);
fprintf(" OK\n")

h = [];
E = [];

h(1) = 1/41;
h(2) = 1/82;
h(3) = 1/164;
h(4) = 1/328;

E(1) = error_0_norm;
E(2) = error_1_norm;
E(3) = error_2_norm;
E(4) = error_3_norm;

second_order = 1.0*h.^(2.0);

hold on
loglog(h, E);
%loglog(h, second_order);


 p = polyfit(log(h),log(E),1);
 z = polyval(p,log(h));
 loglog(h,exp(z))

 hold off

%% funcs
function mesh = steady_state_solution(mesh_level_str, replace_0)
    dir = "/solutions/" + mesh_level_str + "/time_dict";
    time_list = h5read("../export_data.h5", dir);
    str_time_list = compose('%0.6f', time_list);
    dir = "/solutions/" + mesh_level_str + "/time_data/" + str_time_list(end);
    mesh = h5read("../export_data.h5", dir);
    mesh(mesh==0) = replace_0;
end

function mesh = analytical_mesh(mesh_size, A, B, r_inner, r_outer, replace_0)
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
    mesh(mesh==0) = replace_0;
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

function err = calculate_two_norm(error_mesh)
    tmp = sum(error_mesh.^2, 'all') * 1/(size(error_mesh, 1) - 1);
    err = tmp^(1/2);
end