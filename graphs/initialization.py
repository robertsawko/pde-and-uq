from aux import set_plt_params
import numpy as np
import matplotlib.pyplot as plt


N = 8
pol_orders = np.arange(3, N + 1)
x = np.linspace(0, 2*np.pi, 1000)
m = 4
f1 = 0.0
for k in np.arange(1, m + 1):
    f1 = f1 + np.random.randn() * np.sin(k*x) + np.random.randn() * np.cos(k*x)


f1 = f1/np.sqrt(m);

m = 40
f2 = 0.0

for k in np.arange(1, m + 1):
    f2 = f2 + np.random.randn() * np.sin(k*x) + np.random.randn() * np.cos(k*x)

f2 = f2/np.sqrt(m);

set_plt_params(relative_fig_width=0.8)

fig = plt.figure()
ax = fig.gca()
ax.set_xticks((0, 0.5 * np.pi, np.pi, 1.5*np.pi, 2*np.pi))
ax.set_xticklabels(("0", "$\\frac{\pi}{2}$", "$\pi$", "$\\frac{3\pi}{2}$", "$2\pi$"))
ax.set_xlabel("$x$")
ax.set_ylabel("$u_0$")
l1, = plt.plot(x, f1)
l2, = plt.plot(x, f2)
fig.savefig("initialization.pgf", bbox_inches='tight', transparent=True)

figlegend = plt.figure(figsize=(0.8, 1.1))
l = figlegend.legend(
    [l1, l2], ["m=4", "m=40"], ncol=1, loc='upper right')
l.get_frame().set_facecolor('none')
figlegend.savefig(
    "initialization_legend.pgf", transparent=True)

