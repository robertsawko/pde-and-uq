"""
Functions necessary to perform Karhunen-Loeve expansion of

.. math:: C(x, y) = \exp(-c|x-y|)

on interval :math: `(-a, a)`. Based on the solution given in

Ghanem, R. and Spanos, P. "Stochastic Finite Elements: A Spectral Approach",
1991 Springer, pp 29-33
"""
import numpy as np
from scipy.optimize import brentq


def kl_omega(n, a, c, debug=False):
    """Return frequency of even KL eigen functions.

    The frequency is the solution to

    .. math::
        c - \omega \mathrm {tan} (\omega a) = 0

    :param n: the index of eigenfunction
    :param a: the domain limit is :math: `(-a, a)`
    :param c: is the constant in the exponent
    :param debug: debug flag, prints values of the roots (default False)
    :returns: :math:`\omega` that solves the equation
    :rtype: double
    """
    f = lambda omega: c - omega * np.tan(omega * a)
    T = np.pi / a
    s = (n - 0.5) * T + 1e-9
    if n == 0:
        s = 1e-9
    e = (n + 0.5) * T - 1e-9
    # Invoke root finding algorithm Brent (1973)
    omega0 = brentq(f, s, e)
    if debug:
        print "root no {0} is {1} found in ({2}, {3})".format(n, omega0, s, e)
    return omega0


def kl_omega_star(n, a, c, debug=False):
    """
    Return frequency of odd KL eigen functions

    The frequency is the solution

    .. math:: \omega + c \mathrm{tan}(\omega a) = 0

    :param n: the index of eigenfunction
    :param a: the domain limit is (-a, a)
    :param c: is the constant in the exponent
    :param debug: debug flag, prints values of the roots
    :returns: :math:`\omega` that solves the equation
    :rtype: double
    """
    f = lambda omega: omega + c * np.tan(omega * a)
    T = np.pi / a
    s = (n - 0.5) * T + 1e-9
    e = (n + 0.5) * T - 1e-9
    # Invoke root finding algorithm Brent (1973)
    omega0 = brentq(f, s, e)
    if debug:
        print "root no {0} is {1} found in ({2}, {3})".format(n, omega0, s, e)
    return omega0


def kl_eigenpair(n, a, c):
    """
    Return and eigenvalue and and an eigenfunction corresponding to nth element
    in KL expansion.

    :param n: the index of eigenfunction
    :param a: the domain limit is (-a, a)
    :param c: is the constant in the exponent
    :returns: :math:`(\lambda_n, \Psi_n(x))`
    :rtype: Tuple(double, lambda x:)
    """
    if n % 2 == 0:
        omega_n = kl_omega(n / 2, a, c, debug=False)
        return 2 * c / (omega_n**2 + c**2),\
            lambda x: np.cos(omega_n * x)\
            / np.sqrt(a + np.sin(2*omega_n*a) / (2*omega_n))
    else:
        omega_n = kl_omega_star((n + 1) / 2, a, c, debug=False)
        return 2 * c / (omega_n**2 + c**2),\
            lambda x: np.sin(omega_n*x)\
            / np.sqrt(a - np.sin(2 * omega_n*a) / (2*omega_n))


def kl_coefficients(n, a, c):
    """
    Returns the scaling coefficients and the frequencies of the sin and cos
    function from the KL expansion.

    :param n: the index of eigenfunction
    :param a: the domain limit is (-a, a)
    :param c: is the constant in the exponent
    :returns: :math:`(\sqrt{\lambda_n}/A, \omega)` where :math:`A =(1\pm\mathrm{sin}(2\omega a)/(2\omega))`
    :rtype: Tuple(double, double)
    """
    if n % 2 == 0:
        omega_n = kl_omega(n / 2, a, c, debug=False)
        lambda_n = 2 * c / (omega_n**2 + c**2)
        psi_coeff_n = 1 / np.sqrt(a + np.sin(2*omega_n*a) / (2*omega_n))
        return np.sqrt(lambda_n) * psi_coeff_n, omega_n
    else:
        omega_n = kl_omega_star((n + 1) / 2, a, c, debug=False)
        lambda_n = 2 * c / (omega_n**2 + c**2)
        psi_coeff_n = 1 / np.sqrt(a - np.sin(2*omega_n*a) / (2*omega_n))
        return np.sqrt(lambda_n) * psi_coeff_n, omega_n


def kl_expansion(x, m, a=0.5, c=1, eta=None):
    """
    Returns the function object representing the KL expansion.
    """
    w = 0
    if eta is None:
        eta = np.random.randn(m)
    for n in range(m):
        lambda_n, fn = kl_eigenpair(n, a, c)
        w = w + eta[n] * np.sqrt(lambda_n) * fn(x)
    return w


def kl_covariance(x, m, a=0.5, c=1):
    """
    Returns the covariance reconstructed from KL expansion. Added for testing
    purposes.
    """
    x1 = np.outer(x, np.ones(x.size))
    x2 = x1.T
    K = np.zeros((x.size, x.size))
    for n in range(m):
        lambda_n, fn = kl_eigenpair(n, a, c)
        K = K + lambda_n * fn(x1) * fn(x2)
    return x1, x2, K


def kl_expansion_fourier_projected(x, m, n, a=np.pi, c=1, eta=None, phase=0):
    """
    Returns a KL expansion projected on the Fourier series for a given x.
    """
    w = 0
    if eta is None:
        eta = np.random.randn(m)
    for n in range(m):
        psi_n, omega_n = kl_coefficients(n, a, c)
        for k in range(1, n):
            if n % 2 == 0:
                kl_term = eta[n] * psi_n \
                    * (2 * k * np.sin(a*k) * np.cos(a*omega_n)
                        - 2 * omega_n * np.cos(a*k) * np.sin(a*omega_n))\
                    / (k**2 - omega_n**2)
                norm = a + np.sin(2*a*k) / (2*k)
                w = w + kl_term * np.cos(k * x) / norm
            else:
                kl_term = eta[n] * psi_n \
                    * (2 * omega_n * np.sin(a*k) * np.cos(a*omega_n)
                        - 2 * k * np.cos(a*k) * np.sin(a*omega_n))\
                    / (k**2 - omega_n**2)
                norm = a - np.sin(2*a*k) / (2*k)
                w = w + kl_term * np.sin(k * x - phase) / norm
    return w
