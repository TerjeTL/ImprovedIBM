import matplotlib.pyplot as plt

c = (0.5, 0.5)
r_1 = 0.2
r_2 = 0.5

circle1 = plt.Circle(c, r_1, color='r', hatch='//////')
circle2 = plt.Circle(c, r_2, color='blue')

fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()

ax.add_patch(circle1)
ax.add_patch(circle2)

plt.show()