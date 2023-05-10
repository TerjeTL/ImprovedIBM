function val = analytical_solution(t, r, k, r_max)
    val = pyrunfile("bessel_stuff_matlab.py","phi", time=t, r_loc=r, thermal_conductivity=k, r_outer=r_max);
end