function cropped_mesh = crop_boundaries(mesh, r_inner, r_outer)
    [X, Y] = meshgrid_from_mesh(mesh);
    r = sqrt((X-0.5).^2 + (Y-0.5).^2);
    r(r > r_outer) = 0;
    r(r < r_inner) = 0;
    r(r ~= 0) = 1.0;
    cropped_mesh = mesh .* r;
end