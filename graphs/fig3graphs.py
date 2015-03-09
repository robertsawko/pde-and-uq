from aux import set_plt_params
import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.cm as cm


def kinetic_plot_histogram(x, p):
    set_plt_params(relative_fig_width=0.42)
    fig = plt.figure()
    fig.suptitle("PDF from joint-PDF equation", fontsize=7)
    ax = fig.gca()
    ax.set_xticks((0, np.pi, 2 * np.pi))
    ax.set_xticklabels(("0", "$\pi$", "$2\pi$"))
    ax.set_xlabel("$x$")
    ax.set_ylabel("$a$")
    plt.tricontourf(
        x[:, 0], x[:, 1] - 4, p,
        cmap=cm.coolwarm, levels=np.linspace(-0.2,1.2,40))
    plt.subplots_adjust(top=0.8)
    fig.savefig("histogram-kin.pgf", bbox_inches='tight', transparent=True)


def kinetic_extract_mean_and_std(x, p):
    Xb = np.sort(np.unique(x[:, 0]))
    mean = np.zeros(len(Xb))
    std = np.zeros(len(Xb))

    for i, xb in enumerate(Xb):
        ind = np.argwhere(xb == x[:, 0])
        a = x[ind, 1] - 4
        pdf = p[ind]
        # print np.concatenate((a, pdf), axis=1)
        mean[i] = np.inner(a.T, pdf.T) / np.sum(pdf)
        centralised_squarred_a = (a - mean[i])**2
        std[i] = np.inner(centralised_squarred_a.T, pdf.T) / np.sum(pdf)

    return Xb, mean, std


def mc_plot_histogram(X, A, H):
    set_plt_params(relative_fig_width=0.42)
    fig = plt.figure()
    fig.suptitle("PDF from Monte-Carlo", fontsize=7)
    ax = fig.gca()
    ax.set_xticks((0, np.pi, 2*np.pi))
    ax.set_xticklabels(("0", "$\pi$", "$2\pi$"))
    ax.set_xlabel("$x$")
    ax.set_ylabel("$a$")
    plt.contourf(
        X, A, H,
        cmap=cm.coolwarm)
    plt.subplots_adjust(top=0.8)
    fig.savefig("histogram-mc.pgf", bbox_inches='tight', transparent=True)


def mc_extract_mean_and_std(X, A, H):
    Xb = np.sort(np.unique(X[0, :]))
    mean = np.zeros(len(Xb))
    std = np.zeros(len(Xb))

    for i, xb in enumerate(Xb):
        a = A[:, i]
        pdf = H[:, i]
        # print np.concatenate((a, pdf), axis=1)
        mean[i] = np.inner(a.T, pdf.T) / np.sum(pdf)
        centralised_squarred_a = (a - mean[i])**2
        std[i] = np.inner(centralised_squarred_a.T, pdf.T) / np.sum(pdf)

    return Xb, mean, std


def plot_mean_and_std_comparison(X, A, H, x, p):
    mc_x, mc_mean, mc_std = mc_extract_mean_and_std(X, A, H)
    kin_x, kin_mean, kin_std = kinetic_extract_mean_and_std(x, p)
    with open('zero.pickle') as f:
        zero_u, zero_x = pickle.load(f)

    set_plt_params(relative_fig_width=0.41)

    fig = plt.figure()
    fig.suptitle("Mean", fontsize=7)
    ax = fig.gca()
    ax.set_xticks((0, np.pi, 2*np.pi))
    ax.set_xticklabels(("0", "$\pi$", "$2\pi$"))
    ax.set_xlabel("$x$")
    ax.set_ylabel("$a$")
    ax.set_ylim([-1.1, 1.1])
    l1, = plt.plot(mc_x, mc_mean, '--')
    l2, = plt.plot(kin_x, kin_mean)
    l3, = plt.plot(zero_x, zero_u)
    plt.subplots_adjust(top=0.8)
    fig.savefig("compare_mean.pgf", bbox_inches='tight', transparent=True)

    fig = plt.figure()
    fig.suptitle("Standard deviation", fontsize=7)
    ax = fig.gca()
    ax.set_xticks((0, np.pi, 2*np.pi))
    ax.set_xticklabels(("0", "$\pi$", "$2\pi$"))
    ax.set_xlabel("$x$")
    ax.set_ylabel("$a$")
    plt.plot(mc_x, mc_std, '--')
    plt.plot(kin_x, kin_std)
    plt.subplots_adjust(top=0.8)
    fig.savefig("compare_std.pgf", bbox_inches='tight', transparent=True)

    figlegend = plt.figure(figsize=(3.5, 0.5))
    l = figlegend.legend(
        [l1, l2, l3],
        ['MC', 'Joint PDF', 'Naive solution'],
        ncol=3, loc='upper right')
    l.get_frame().set_facecolor('none')
    figlegend.savefig(
        "compare_legend.pgf", transparent=True)

with open('mc_histogram.pickle') as f:
    X, A, H = pickle.load(f)

with open('xp.pickle') as f:
    x, p = pickle.load(f)

mc_plot_histogram(X, A, H)
kinetic_plot_histogram(x, p)

plot_mean_and_std_comparison(X, A, H, x, p)
