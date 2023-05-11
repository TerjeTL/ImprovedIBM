function mesh = analytical_transient_mesh(mesh_size, time, r_outer, replace_0)
    [X, Y] = meshgrid_from_mesh(mesh_size);
    r = sqrt((X-0.5).^2 + (Y-0.5).^2);
    r(r > r_outer) = replace_0;
    r(r ~= replace_0) = arrayfun(@(r) analytical_transient_solution(time, r, 2.0, r_outer), r(r ~= replace_0));
    mesh = r;
end