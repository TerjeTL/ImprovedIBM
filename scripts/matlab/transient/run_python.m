%pyenv('Version', "C:/Users/forme/AppData/Local/Programs/Python/Python39/python.exe")

clc; clear; close all;
%res = pyrunfile("bessel_stuff_matlab.py","phi",time=0.1, r=0.8, thermal_conductivity=1.0, r_outer=1.0);

xs = linspace(0, 1, 1000);

ys = arrayfun(@(x) analytical_solution(0.5, x, 1.0, 1.0), xs);

plot(xs, ys);

%print(res)



% 
% print(res);

%res = pyrunfile("addac.py","z",x=3,y=2);

function val = analytical_solution(t, r, k, r_max)
    val = pyrunfile("bessel_stuff_matlab.py","phi", time=t, r_loc=r, thermal_conductivity=k, r_outer=r_max);
end