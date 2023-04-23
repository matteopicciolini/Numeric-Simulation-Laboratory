#matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#numpy
import numpy as np

fig = plt.figure(figsize=plt.figaspect(0.5))

n, x, y, z = np.loadtxt('../Data/02.2_distr_non_unif_sphere.dat', unpack=True, usecols=(0, 1, 2, 3))
ax1= fig.add_subplot(1, 2, 1, projection='3d')

u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:20j]
xx = np.cos(u)*np.sin(v)
yy = np.sin(u)*np.sin(v)
zz = np.cos(v)
ax1.plot_wireframe(xx, yy, zz, lw=0.3, color = "gray")

ax1.scatter(x,y,z, s=1)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("z")
ax1.set_box_aspect([1,1,1])

n, x, y, z = np.loadtxt('../Data/02.2_distr_unif_sphere.dat', unpack=True, usecols=(0, 1, 2, 3))
ax2 = fig.add_subplot(1, 2, 2, projection='3d')

ax2.plot_wireframe(xx, yy, zz, lw=0.3, color = "gray")
ax2.scatter(x,y,z, s=1)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("z")
ax2.set_box_aspect([1,1,1])
plt.show()