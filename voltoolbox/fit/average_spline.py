import numpy as np
import scipy.integrate as integrate
from voltoolbox import AverageSplineCurve


def avg(f, z):
    if abs(z) < 1.0e-15:
        return f(0.0)
    x = integrate.romberg(f, 0.0, z)
    return x / z


def put_wing_func(a, b, c):
    assert (a < b) & (b < c)

    def f(z):
        res = (max(0, c - z)**3 - max(0, b - z)**3) / (c - b)
        res -= (max(0, b - z)**3 -max(0, a - z)**3) / (b - a)
        return res / 6.0

    return f


def call_wing_func(a, b, c):
    assert (a < b) & (b < c)

    def f(z):
        res = (max(0, z - a)**3 -max(0, z - b)**3) / (b - a)
        res -= (max(0, z - b)**3 -max(0, z - c)**3) / (c - b)
        return res / 6.0
    return f


def convex_func(a, b, c):
    assert (a < b) & (b < c)

    def f(z):
        res = (max(0, z - a)**3 - max(0, z - b)**3) / (b - a)
        res -= (max(0, z - b)**3 - max(0, z - c)**3) / (c - b)
        res -= (b - a)**2

        res += (max(0, c - z)**3 - max(0, b - z)**3) / (c - b)
        res -= (max(0, b - z)**3 - max(0, a - z)**3) / (b - a)
        res -= (c-b)**2
        return res / 12.0
    return f


def avg_spline_crv(nodes, f):
    return AverageSplineCurve(nodes, [f(x) for x in nodes])


def avg_convex_func(a, b, c):
    f = convex_func(a, b, c)
    return avg_spline_crv((a, b, c), f)


def avg_put_wing_func(a, b, c):
    f = put_wing_func(a, b, c)
    return avg_spline_crv((a, b, c), f)


def avg_call_wing_func(a, b, c):
    f = call_wing_func(a, b, c)
    return avg_spline_crv((a, b, c), f)


class AverageSpline(object):

    def __init__(self, basis_funcs, basis_penalty):
        assert len(basis_funcs) == len(basis_penalty)
        self.basis_funcs = basis_funcs
        self.basis_penalty = basis_penalty

    @classmethod
    def average_spline(cls, put_zs=None, call_zs=None):

        basis_funcs = [
            np.vectorize(lambda z: 1.0),
            np.vectorize(lambda z: z),
            np.vectorize(avg_convex_func(*(-2.5, 0.0, 2.5)).vol),
        ]

        basis_penalty = [0.0, 0.0, 0.0]

        if put_zs:
            for i in range(len(put_zs) - 1):
                if i == 0:
                    a, b, c = (-put_zs[1], -put_zs[0], 0.0)
                else:
                    a, b, c = (-put_zs[i+1], -put_zs[i], -put_zs[i-1])

                basis_funcs += [np.vectorize(avg_put_wing_func(*(a, b, c)).vol)]
                basis_penalty += [1.0]

        if call_zs:
            for i in range(len(call_zs) - 1): 
                if i==0:
                    a, b, c = (0.0, call_zs[0], call_zs[1])
                else:
                    a, b, c = (call_zs[i-1], call_zs[i], call_zs[i+1])
                basis_funcs += [np.vectorize(avg_call_wing_func(*(a, b, c)).vol)]
                basis_penalty += [1.0]
        return cls(basis_funcs, basis_penalty)

    def sample_basis(self, xs: np.array):
        return np.column_stack( [f(xs) for f in self.basis_funcs])

    def fit_coeffs(self, xs, vals, errs=None,
                    *, smoothing=1.0e-12):
        if errs is None:
            errs = np.array([1.0] * xs.shape[0])

        basis_vals = self.sample_basis(xs)
        basis_vals /= np.column_stack([errs])
        target = vals / errs

        var = np.matmul(basis_vals.T, basis_vals)

        basis_penalty = np.array(self.basis_penalty)
        var_regul = np.diag(basis_penalty) * var.trace()
        sum_pen = basis_penalty.sum()
        if abs(sum_pen) > 0.0:
            var_regul /= sum_pen

        var += smoothing * var_regul

        cov = basis_vals.T.dot(target)
        return np.linalg.solve(var, cov)
