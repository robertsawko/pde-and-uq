from numpy import linspace, pi, sin, var
from synthesis import shinozuka
import pickle
from numpy.random import seed, rand, randn


if __name__ == '__main__':
    L = 2 * pi
    N = 40000
    x = linspace(0, L, N)

    seed(123)
    b = 100
    M = 16384
    # Phi = rand(maxM) * 2 * pi
    target_energy = 0.5

    eta0 = randn(1)
    Phi = rand(M) * 2 * pi
    noise = shinozuka(x, b=b, delta_w=1 / L, Phi=Phi[0:M])
    y = 0.5 * (1.0 + eta0 + sin(x)) + noise
    with open('u0.pickle', 'wb') as f:
        pickle.dump([x, y], f)
