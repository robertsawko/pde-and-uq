from numpy import linspace, logspace, pi, cos, sqrt, zeros_like
from numpy.random import seed, rand
from matplotlib.pyplot import figure
from scipy.signal import welch


def ftransform_analytical(omega, b):
    return 1 / pi * b / (b**2 + omega**2)


def shinozuka(x, b, Phi, delta_w=1 / 2 / pi):
    '''
    Shinozuka process-generation method

    :param x: spatial position
    :param b: coefficient in the correlation exponent
    :param Phi: random phases
    :param delta_w: frequency increment and smallest resolved frequency

    Based on

        Shinozuka, M. and Deodatis, G. "Simulation of stochastic process by
        spectral representation", 1991 Applied Mechanics Reviews, pp 191-204

    Note: If the smallest frequency is larger than the domain length the
    periodicity may be visible.
    '''
    M = len(Phi)
    f0 = zeros_like(x)
    for k in range(M):
        omega = (k + 1) * delta_w
        Sff = 1 / pi * b / (b**2 + omega**2)
        f0 += sqrt(4 * Sff * delta_w) * cos(omega * x + Phi[k])
    return f0


if __name__ == '__main__':
    L = 2 * pi
    N = 20000
    x = linspace(0, L, N)

    fig_corr = figure()
    fig_sig = figure()
    axs = fig_sig.gca()
    axc = fig_corr.gca()
    seed(123)
    b = 100
    maxM = 8192
    Phi = rand(maxM) * 2 * pi

    for m in [maxM // 4, maxM // 2, maxM]:
        y = shinozuka(x, b=b, delta_w=1 / L, Phi=Phi[0:m])
        # larger nperseg makes smaller frequencies visible
        f, pxx = welch(y, fs=1 / (x[1] - x[0]), nperseg=N // 10)

        axc.loglog(
            2 * pi * f, pxx / 4 / pi,
            label='Empirical m={0:n}'.format(m))
        axs.plot(x, y, label='m={0:n}'.format(m))

    omega = logspace(-1, 4, 10000)
    a = ftransform_analytical(omega, b=b)
    axc.loglog(omega, a, label='Analytical')
    axs.legend()
    axc.legend()
    fig_sig.show()
    fig_corr.show()
