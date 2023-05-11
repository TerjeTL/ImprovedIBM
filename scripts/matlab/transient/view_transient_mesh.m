clc; clear; close all;
addpath('./../functions/');



%FILE = "../seven_levels_ss.h5";
%FILE = "../export_data_high_res.h5";
FILE = "./../../export_data.h5";
GROUP_SIZE = "89";
TIME = 3;

alpha_x = 2.0; %order in space
alpha_t = 1.0; %order in time

richardson_extrp_dir = "/solutions/richardson_extrp/transient/";
solution_dir = richardson_extrp_dir + "group_size_" + GROUP_SIZE + "/";

%% Get Geometric Params
% radius, x, y, bc_phi
%raw_data = h5read(FILE, "/geometry/Inner Circle/data");
%r_inner = raw_data(1);
%bc_phi_inner = raw_data(4);
%bc_inner = h5read(FILE, "/geometry/Inner Circle/bc");

raw_data = h5read(FILE, "/geometry/Outer Circle/data");
r_outer = raw_data(1);
bc_phi_outer = raw_data(4);
bc_outer = h5read(FILE, "/geometry/Outer Circle/bc");

% Get List of Logged Time Levels
times = get_hdf5_subdirs(FILE, "solutions/richardson_extrp/transient/group_size_" + GROUP_SIZE + "/");

% Get Solutions of First Time Level
mesh_0 = h5read(FILE, solution_dir + times(TIME) + "/coarse");
mesh_1 = h5read(FILE, solution_dir + times(TIME) + "/fine");

analytical_mesh = analytical_transient_mesh(mesh_0, str2double(times(TIME)), r_outer, nan);

%mesh_exact = h5read(FILE, richardson_extrp_dir + "analytical_size_45/" + times(TIME));
mesh_re = h5read(FILE, solution_dir + times(TIME) + "/richardson_extrp");

mesh_0(mesh_0==0) = nan;
mesh_1(mesh_1==0) = nan;
%mesh_exact(mesh_exact==0) = nan;
mesh_re(mesh_re==0) = nan;


completed_re = zeros(size(mesh_1));
for i = 1 : (size(mesh_1, 1)-2)
  for j = 1 : (size(mesh_1, 2)-2)
    I = (i+1)/2;
    J = (j+1)/2;
    
    if (rem(i+1,2) == 0 && rem(j, 2) ~= 0)  
        % even i+1, odd j
        c1 = mesh_re(I,J) - mesh_1(i,j);
        c2 = mesh_re(I+1,J) - mesh_1(i+2,j);
        c = 1.0/2.0 * (c1 + c2);
    
        completed_re(i+1, j) = mesh_1(i+1, j) + c;
    end
    if (rem(i,2) ~= 0 && rem(j+1, 2) == 0)
        % even j+1, odd i
        c1 = mesh_re(I,J) - mesh_1(i,j);
        c2 = mesh_re(I,J+1) - mesh_1(i,j+2);
        c = 1.0/2.0 * (c1 + c2);
    
        completed_re(i, j+1) = mesh_1(i, j+1) + c;
    end
    if (rem(i+1,2) == 0 && rem(j+1, 2) == 0)
        % statements even i+1, j+1
        c1 = mesh_re(I,J) - mesh_1(i,j);
        c2 = mesh_re(I+1,J) - mesh_1(i+2,j);
        c3 = mesh_re(I,J+1) - mesh_1(i,j+2);
        c4 = mesh_re(I+1,J+1) - mesh_1(i+2,j+2);

        c = 1.0/4.0 * (c1 + c2 + c3 + c4);

        completed_re(i+1, j+1) = mesh_1(i+1, j+1) + c;
    end
    if (rem(i,2) ~= 0 && rem(j, 2) ~= 0)
        completed_re(i, j) = mesh_re(I, J);
    end 
  end
end

ones = mesh_1;
ones(ones ~= 0) = 1.0;
mesh_re_fine = h5read(FILE, solution_dir + times(TIME) + "/richardson_extrp_fine");
mesh_re_fine(mesh_re_fine==0) = nan;
%mesh_re_fine = mesh_re_fine.*ones;

analytical_1 = analytical_transient_mesh(mesh_1, str2double(times(TIME)), r_outer, nan);

error_0 = mesh_0 - analytical_mesh;
error_1 = mesh_1 - analytical_1;
error_re = mesh_re - analytical_mesh;
error_completed = completed_re - analytical_1;


%% Plots
tiledlayout(2,3)

% level 0
nexttile
plot_0 = plot_mesh_surface(mesh_0);
%zlim([0 2]);

nexttile
plot_1 = plot_mesh_surface(mesh_1);
%zlim([0 2]);

nexttile
plot_re = plot_mesh_surface( error_completed);

nexttile
plot_err_0 = plot_mesh_surface(error_0);

nexttile
plot_err_1 = plot_mesh_surface(error_1);

nexttile
plot_err_re = plot_mesh_surface(error_re);

set(gcf,'position',[get(0, 'Screensize')])

%% funcs
function s = plot_mesh_surface(mesh)
    [X, Y] = meshgrid_from_mesh(mesh);
    s = surf(X,Y,mesh,FaceColor='interp');
end

function directories = get_hdf5_subdirs(file, dir)
    s = struct;
    s.datasets = {};
    s.groups = {}; 
    fid = H5F.open(file, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    anotherid = H5O.open(fid, dir, 'H5P_DEFAULT');
    [status, sOut] = H5O.visit(anotherid, 'H5_INDEX_NAME', 'H5_ITER_NATIVE', @getObjInfo, s);
    
    directories = sOut.groups(2:end);
end

function [status, sOut] = getObjInfo(rootID, name, sIn)
    objID = H5O.open(rootID, name, 'H5P_DEFAULT');
    obj_info=H5O.get_info(objID);
    H5O.close(objID);

    switch(obj_info.type)
        case H5ML.get_constant_value('H5G_GROUP')
            sIn.groups{end+1} = name;
    end
    status = 0; % keep iterating
    sOut = sIn; 
end

function val = fine_grid_re(re, fine, i, j)
    I = (i+1)/2;
    J = (j+1)/2;
    if (rem(i+1,2) == 0 && rem(j, 2) ~= 0)  
        % even i+1, odd j
        c1 = re(I,J) - fine(i,j);
        c2 = re(I+1,J) - fine(i+2,j);
        c = 1.0/2.0 * (c1 + c2);

        mat(i+1, j) = fine(i+1, j) + c;
    elseif (rem(i,2) ~= 0 && rem(j+1, 2) == 0)
        % even j+1, odd i
        c1 = re(I,J) - fine(i,j);
        c2 = re(I,J+1) - fine(i,j+2);
        c = 1.0/2.0 * (c1 + c2);

        mat(i, j+1) = fine(i, j+1) + c;
    elseif (rem(i+1,2) == 0 && rem(j+1, 2) == 0)
        % statements even i, j
        c1 = re(I,J) - fine(i,j);
        c2 = re(I+1,J) - fine(i+2,j);
        c3 = re(I,J+1) - fine(i,j+2);
        c4 = re(I+1,J+1) - fine(i+2,j+2);

        c = 1.0/4.0 * (c1 + c2 + c3 + c4 + c5);

        mat(i+1, j+1) = fine(i+1, j+1) + c;
    elseif (rem(I,2) == 0 && rem(J, 2) == 0)
        mat(i, j) = re(I, J);
    end 
end
