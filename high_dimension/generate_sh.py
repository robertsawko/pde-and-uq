from numpy import linspace, pi
from synthesis import shinozuka
import pickle
from numpy.random import seed, rand


if __name__ == '__main__':
    L = 2 * pi
    N = 40000
    x = linspace(0, L, N)

    seed(123)
    b = 100
    M = 16384
    # Phi = rand(maxM) * 2 * pi
    target_energy = 0.5

    Phi = rand(M) * 2 * pi
    y = shinozuka(x, b=b, delta_w=1 / L, Phi=Phi[0:M])

    with open('u0.pickle', 'wb') as f:
        pickle.dump([x, y], f)
