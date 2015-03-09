from aux import set_plt_params
import numpy as np
import matplotlib.pyplot as plt
import pickle
from itertools import cycle


def plot_mean_and_std_comparison():
    mesh_sizes = [50, 75, 100]
    dgp1 = dict((
        m, np.genfromtxt("pdf_dgp1_{0:03d}.csv".format(m),
        delimiter=',', skip_header=1, skip_footer=1)
        ) for m in mesh_sizes
    )
    dgp2 = dict((
        m, np.genfromtxt("pdf_dgp2_{0:03d}.csv".format(m),
        delimiter=',', skip_header=1, skip_footer=1)
        ) for m in mesh_sizes
    )

    def opacity_fun(x):
        return np.max([0.2, x])

    set_plt_params(relative_fig_width=0.41)

    fig = plt.figure()
    fig.suptitle("DG Quadratic", fontsize=7)
    ax = fig.gca()
    ax.set_xlabel("$a$")
    ax.set_ylabel("$p(a, x=\pi, t=1)$")
    lines = []
    labels = []
    line_styles = [":", "--", "-"]
    linestyle_cycler = cycle(line_styles)
    for m in mesh_sizes:
        l, = plt.plot(dgp2[m][:,9] - 4, dgp2[m][:, 1],
            linestyle= next(linestyle_cycler),
            alpha=0.5)
        lines.append(l)
        labels.append('Mesh {0}'.format(m))

    plt.subplots_adjust(top=0.8)
    fig.savefig("pdf_mesh_convergence_dgp2.pgf", bbox_inches='tight', transparent=True)

    fig = plt.figure()
    fig.suptitle("DG Linear", fontsize=7)
    ax = fig.gca()
    ax.set_xlabel("$a$")
    ax.set_ylabel("$p(a, x=\pi, t=1)$")
    for m in mesh_sizes:
        plt.plot(
            dgp1[m][:,9] - 4, dgp1[m][:, 1],
            linestyle= next(linestyle_cycler),
            alpha=0.5)

    plt.subplots_adjust(top=0.8)
    fig.savefig("pdf_mesh_convergence_dgp1.pgf", bbox_inches='tight', transparent=True)

    figlegend = plt.figure(figsize=(3.5, 0.6))
    l = figlegend.legend(
        lines, labels, ncol=3, loc='upper right')
    l.get_frame().set_facecolor('none')
    figlegend.savefig(
        "pdf_mesh_convergence-legend.pgf", transparent=True)

plot_mean_and_std_comparison()
