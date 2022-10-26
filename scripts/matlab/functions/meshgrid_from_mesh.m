function [X, Y] = meshgrid_from_mesh(mesh)
    x = linspace(0, 1, size(mesh, 1));
    y = linspace(0, 1, size(mesh, 2));
    
    [X, Y] = meshgrid(x, y);
end