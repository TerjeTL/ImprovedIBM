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

%% Example grid
x = linspace(0, 1, 100);
y = linspace(0, 1, 100);

[X, Y] = meshgrid(x, y);
r = sqrt((X-0.5).^2 + (Y-0.5).^2);

%% Calculate Z-Values for this Mesh
analytical_mesh = zeros(size(r));
for i = 1:size(analytical_mesh,1)
    for j = 1:size(analytical_mesh,2)
        if (r(i,j) > r_inner && r(i,j) < r_outer)
            analytical_mesh(i,j) = analytical_value(A, B, r(i,j));
        end
    end
end
analytical_mesh(analytical_mesh==0) = nan;

%% continues...
tiledlayout(1,2)
nexttile
surf(X,Y,analytical_mesh);
zlim([0 250]);
nexttile
surf(X,Y,analytical_mesh);
zlim([0 250]);
