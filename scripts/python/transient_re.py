import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import h5py
import pandas as pd

file = h5py.File("./../export_data.h5", 'r')


re_group_dset = file['solutions/richardson_extrp/transient/group_size_45']
times = [x for x in re_group_dset.keys()]



print(re_group_dset)
print(times)

fine_grid = file['solutions/richardson_extrp/transient/group_size_45/' + times[0] + '/fine'][...]

print(fine_grid)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

xs = np.linspace(0, 1, np.shape(fine_grid)[0])
ys = np.linspace(0, 1, np.shape(fine_grid)[1])
X, Y = np.meshgrid(xs, ys)

surf = ax.plot_surface(X, Y, fine_grid, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()