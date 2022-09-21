clc; clear; close all;

load("first_data.csv");

surf(min(first_data, 200.0));