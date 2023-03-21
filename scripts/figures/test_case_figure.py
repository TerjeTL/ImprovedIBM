import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

plt.figure(figsize=(8,8))

SMALL_SIZE = 18
MEDIUM_SIZE = 24
BIGGER_SIZE = 28

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

n_cell_x = 22
n_cell_y = n_cell_x

x, y = np.meshgrid(np.linspace(0,1, n_cell_x+1), np.linspace(0, 1, n_cell_y+1))

segs1 = np.stack((x,y), axis=2)
segs2 = segs1.transpose(1,0,2)
plt.gca().add_collection(LineCollection(segs1, colors='black', zorder=0))
plt.gca().add_collection(LineCollection(segs2, colors='black', zorder=0))

plt.scatter(x, y, s=30, c='b', zorder=1)

# geometry
r_inner = 0.15
r_outer = 0.45

# the two circles
thetas = np.linspace(0,2*np.pi, 200)
# you don't need r = np.one(len(thetas))

x_unitcirc = np.cos(thetas)
y_unitcirc = np.sin(thetas)

x_outer = np.cos(thetas)*r_outer
y_outer = np.sin(thetas)*r_outer

x_inner = x_unitcirc * r_inner
y_inner = y_unitcirc * r_inner

xs = np.linspace(-1,1, 201)
ys = np.linspace(-1,1, 201)

# mesh for contours
xv,yv = np.meshgrid(xs,ys)

# generate the level map
r = xv**2 + yv**2

# transpose
xv += 0.5
yv += 0.5

x_outer += 0.5
y_outer += 0.5

x_inner += 0.5
y_inner += 0.5

# plot the contours with two levels only
# notice the xv, yv parameters
plt.contourf(xv, yv, r, levels=[r_outer**2, 2**2], colors=('white'), hatches=['/'])

plt.contourf(xv, yv, r, levels=[0.0**2, r_inner**2], hatches=['/'], colors='white')

# plot the two circles
plt.plot(x_outer, y_outer, color='r', linewidth=3)
plt.plot(x_inner, y_inner, color='r', linewidth=3)

plt.text(0.375, 0.505, '$T = 1$', fontsize = 24, rotation=45, 
        bbox = dict(facecolor = 'white', edgecolor='white'))

plt.text(0.10, 0.805, '$T = 2$', fontsize = 24, rotation=45, 
        bbox = dict(facecolor = 'white', edgecolor='white'))

plt.xlim([0, 1])
plt.ylim([0, 1])

plt.xlabel('x')
plt.ylabel('y')

plt.show()