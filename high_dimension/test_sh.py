'''
Testing Shinozuka method for signal generation.
'''
from numpy import logspace, linspace, pi
from numpy.random import seed, rand
from matplotlib.pyplot import figure
from scipy.signal import periodogram
from scipy.integrate import trapz
from synthesis import ftransform_analytical, shinozuka

if __name__ == '__main__':
    L = 2 * pi
    N = 20000
    x = linspace(0, L, N)
    repeats = 100

    fig_corr = figure()
    fig_sig = figure()
    axs = fig_sig.gca()
    axc = fig_corr.gca()
    seed(123)
    b = 100
    maxM = 2**12
    target_energy = 0.5

    for m in [maxM // 4, maxM // 2, maxM]:
        # larger nperseg makes smaller frequencies visible
        Sxxs = []
        for n in range(repeats):
            Phi = rand(m) * 2 * pi
            y = shinozuka(x, b=b, delta_w=1 / L, Phi=Phi[0:m])
            f, pxx = periodogram(y, fs=1 / (x[1] - x[0]))
            omega = 2 * pi * f
            Sxxs.append(pxx / 4 / pi)
        Sxx = sum(Sxxs) / repeats
        axc.loglog(omega, Sxx, label='Empirical m={0:n}'.format(m))
        axs.plot(x, y, label='m={0:n}'.format(m))
        print('Captured energy percentage with {1:d} modes: {0:0.1f}%'.format(
            trapz(Sxx, x=omega) / target_energy * 100, m))

    omega = logspace(-1, 4, 10000)
    a = ftransform_analytical(omega, b=b)
    axc.set_ylim([10**-6, 5 * 10**-3])
    axc.loglog(omega, a, label='Analytical')
    axs.legend()
    axc.legend()
    fig_sig.show()
    fig_corr.show()
