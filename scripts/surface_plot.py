# import numpy as np

# import pyvista as pv
# from pyvista import examples

# # Make data
# x = np.arange(-10, 10, 0.25)
# y = np.arange(-10, 10, 0.25)
# x, y = np.meshgrid(x, y)
# r = np.sqrt(x**2 + y**2)
# z = np.sin(r)

# # Create and plot structured grid
# grid = pv.StructuredGrid(x, y, z)
# grid.plot()

import numpy as np
import pandas as pd
import h5py
import pyvista as pv

with h5py.File('export_data.h5', 'r') as hf:
    t = hf['solutions/mesh_0/time_dict'][:]

x = np.arange(-10, 10, 0.5)
y = np.arange(-10, 10, 0.5)
x, y = np.meshgrid(x, y)
r = np.sqrt(x**2 + y**2)
z = np.sin(r)*10

# Create and structured surface
grid = pv.StructuredGrid(x, y, z)

# Create a plotter object and set the scalars to the Z height
plotter = pv.Plotter(notebook=False, off_screen=False)
plotter.add_mesh(
    grid,
    scalars=z.ravel(),
    lighting=False,
    show_edges=True,
    scalar_bar_args={"title": "Height"},
    clim=[-5, 5]
)

# Open a gif
#plotter.open_gif("wave.gif")

pts = grid.points.copy()

plotter.show_grid()
plotter.show(interactive_update=True)

# Update Z and write a frame for each updated position
nframe = 15000
for phase in np.linspace(0, 2 * np.pi, nframe + 1)[:nframe]:
    z = np.sin(r + phase*15)*5
    pts[:, -1] = z.ravel()
    plotter.update_coordinates(pts, render=False)
    plotter.update_scalars(z.ravel(), render=False)

    # Write a frame. This triggers a render.
    plotter.update()
    #plotter.write_frame()

# Closes and finalizes movie
plotter.close()