function mesh = analytical_mesh(mesh_size, A, B, r_inner, r_outer, replace_0)
    [X, Y] = meshgrid_from_mesh(mesh_size);
    r = sqrt((X-0.5).^2 + (Y-0.5).^2);
    r(r > r_outer) = replace_0;
    r(r < r_inner) = replace_0;
    r(r ~= replace_0) = analytical_value(double(A),double(B),r(r ~= replace_0));
    mesh = r;
end