from aux import set_plt_params
import numpy as np
import matplotlib.pyplot as plt
import pickle


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


def plot_mean_and_std_comparison(X, A, H):
    mc_x, mc_mean, mc_std = mc_extract_mean_and_std(X, A, H)
    N = 8
    pol_orders = np.arange(3, N + 1)
    gpc = dict(
        (
            k,
            np.genfromtxt(
                "gpcN{0}.csv".format(k),
                delimiter=',', skip_header=1, skip_footer=1)
        ) for k in pol_orders)

    def opacity_fun(x):
        return 0.95 * x**3 + 0.05

    set_plt_params(relative_fig_width=0.41)

    fig = plt.figure()
    fig.suptitle("Mean", fontsize=7)
    ax = fig.gca()
    ax.set_xticks((0, np.pi, 2*np.pi))
    ax.set_xticklabels(("0", "$\pi$", "$2\pi$"))
    ax.set_xlabel("$x$")
    ax.set_ylabel("$a$")
    plt.plot(mc_x, mc_mean, '--')
    for n in pol_orders:
        plt.plot(
            gpc[n][:, 6 + n], gpc[n][:, 0],
            color='#348ABD',
            alpha=opacity_fun(n/float(N)))
    plt.subplots_adjust(top=0.8)
    fig.savefig("gpc_mean.pgf", bbox_inches='tight', transparent=True)

    fig = plt.figure()
    fig.suptitle("Standard deviation", fontsize=7)
    ax = fig.gca()
    ax.set_xticks((0, np.pi, 2*np.pi))
    ax.set_xticklabels(("0", "$\pi$", "$2\pi$"))
    ax.set_xlabel("$x$")
    ax.set_ylabel("$a$")
    l1, = plt.plot(mc_x, mc_std, '--', label='MC')
    lines = [l1]
    labels = ['MC']
    for n in pol_orders:
        l, = plt.plot(
            gpc[n][:, 6 + n], np.sum(gpc[n][:,0:n]**2, axis=1),
            color='#348ABD',
            alpha=opacity_fun(float(n) / N)
        )
        lines.append(l)
        labels.append('gPC N={0}'.format(n))

    plt.subplots_adjust(top=0.8)
    fig.savefig("gpc_std.pgf", bbox_inches='tight', transparent=True)

    figlegend = plt.figure(figsize=(2.3, 0.6))
    l = figlegend.legend(
        lines, labels, ncol=2, loc='upper right')
    l.get_frame().set_facecolor('none')
    figlegend.savefig(
        "gpc_legend.pgf", transparent=True)

with open('mc_histogram.pickle') as f:
    X, A, H = pickle.load(f)

plot_mean_and_std_comparison(X, A, H)
