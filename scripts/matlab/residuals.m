clc; clear; close all;

%% Get Residual Data
residual_0 = load_residual("mesh_0");
residual_1 = load_residual("mesh_1");
residual_2 = load_residual("mesh_2");
residual_3 = load_residual("mesh_3");

%residual_0(1,2) = residual_0(1,2);
%residual_1(1,2) = residual_1(1,2)/2;
%residual_2(1,2) = residual_2(1,2)/4;
%residual_3(1,2) = residual_3(1,2)/8;

residual_0(2:end,2) = residual_0(2:end,2);
residual_1(2:end,2) = residual_1(2:end,2);
residual_2(2:end,2) = residual_2(2:end,2);
residual_3(2:end,2) = residual_3(2:end,2);

%residual_0(:,2) = residual_0(:,2)/residual_0(1,2);
%residual_1(:,2) = residual_1(:,2)/residual_1(1,2);
%residual_2(:,2) = residual_2(:,2)/residual_2(1,2);
%residual_3(:,2) = residual_3(:,2)/residual_3(1,2);

hold on
plot(residual_0(:,1), residual_0(:,2));
plot(residual_1(:,1), residual_1(:,2));
plot(residual_2(:,1), residual_2(:,2));
plot(residual_3(:,1), residual_3(:,2));
hold off
set(gca, 'XScale', 'linear', 'YScale', 'log');

xlabel('h','interpreter', 'latex', 'FontSize', 24) 
ylabel('$|\!|\mathrm{T-T_{ex}}|\!|_2$', 'interpreter', 'latex', 'FontSize', 24);
% Create legend
%legend1 = legend('Solution', '1. order','2. order', '4. order', 'FontSize', 24);
%set(legend1, 'Position', ...
%    [0.630952386221006 0.142460322001624 0.267857137588518 0.165476185934884]);


function residual = load_residual(mesh_level_str)
    dir = "/solutions/" + mesh_level_str + "/time_dict";
    time_list = h5read("../export_data.h5", dir);
    dir = "/solutions/" + mesh_level_str + "/euclidian_norm";
    res_list = h5read("../export_data.h5", dir);
    residual = cat(2, time_list, res_list);
end

function plot_residual(mesh_level_str)
    residual = load_residual(mesh_level_str);
    plot(residual(:,1), residual(:,2));
end