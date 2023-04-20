clc; clear; close all;
addpath('./../functions/');



%FILE = "../seven_levels_ss.h5";
%FILE = "../export_data_high_res.h5";
FILE = "./../../export_data.h5";
GROUP_SIZE = "89";
TIME = 3;

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
mesh_exact = h5read(FILE, richardson_extrp_dir + "analytical_size_89/" + times(TIME));
mesh_re = h5read(FILE, solution_dir + times(TIME) + "/richardson_extrp");

mesh_0(mesh_0==0) = nan;
mesh_1(mesh_1==0) = nan;
mesh_exact(mesh_exact==0) = nan;
mesh_re(mesh_re==0) = nan;


error = mesh_1 - mesh_exact;

%% Plots
tiledlayout(1,3)

% level 0
nexttile
plot_0 = plot_mesh_surface(mesh_1);
zlim([0 2]);

nexttile
plot_1 = plot_mesh_surface(mesh_exact);
zlim([0 2]);

nexttile
plot_re = plot_mesh_surface(error);

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