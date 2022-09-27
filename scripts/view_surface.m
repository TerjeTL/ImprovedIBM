clc; clear; close all;

phi_data = load("converged_result.csv");

% Remove unwanted datapoints
phi_data(phi_data==0) = nan;

surf(phi_data);