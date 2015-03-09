from aux import set_plt_params
import numpy as np
from kl_expansion import kl_covariance
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

a = 0.5
c = 1
x = np.linspace(-a, a, 100)
x1, x2, K5 = kl_covariance(x=x, m=5, a=a)
x1, x2, K10 = kl_covariance(x=x, m=10, a=a)
K = np.exp(-c * np.abs(x1-x2))

set_plt_params(relative_fig_width=0.99)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(x1, x2, K, rstride=2, cstride=2, alpha=0.3)

cset = ax.contour(x1, x2, K, zdir='z', offset=-0.3, cmap=cm.coolwarm)

ax.set_zlim(-0.3, 1) 
ax.set_zticks(np.linspace(0,1,5)) 
ax.set_xlabel("$x$") 
ax.set_ylabel("$y$") 

plt.savefig("covariance.pdf", bbox_inches='tight', transparent=True)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(x1, x2, K5, rstride=2, cstride=2, alpha=0.3)

cset = ax.contour(x1, x2, K, zdir='z', offset=-0.3, cmap=cm.coolwarm)

ax.set_zlim(-0.3, 1) 
ax.set_zticks(np.linspace(0,1,5)) 
ax.set_xlabel("$x$") 
ax.set_ylabel("$y$") 

plt.savefig("covarianceN5.pdf", bbox_inches='tight', transparent=True)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(x1, x2, K10, rstride=2, cstride=2, alpha=0.3)

cset = ax.contour(x1, x2, K, zdir='z', offset=-0.3, cmap=cm.coolwarm)

ax.set_zlim(-0.3, 1) 
ax.set_zticks(np.linspace(0,1,5)) 
ax.set_xlabel("$x$") 
ax.set_ylabel("$y$") 

plt.savefig("covarianceN10.pdf", bbox_inches='tight', transparent=True)
