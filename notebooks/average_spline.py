import numpy as np
import scipy.integrate as integrate

def avg(f, z):
    if abs(z) < 1.0e-15:
        return f(0.0)
    x, _ = integrate.quad(f, 0.0, z)
    return x / z

class AverageSpline(object):
    
    def __init__(self, put_zs, call_zs):
        assert(len(put_zs) > 1 and len(call_zs) > 1)
        self.put_zs = tuple(put_zs)
        self.call_zs = tuple(call_zs)
    
    def put_wing(self, i, z):
        assert(0 <= i < len(self.put_zs) - 1)   
        if i == 0:
            a, b, c = (-self.put_zs[1], -self.put_zs[0], 0.0)
        else:
            a, b, c = (-self.put_zs[i+1], -self.put_zs[i], -self.put_zs[i-1])
        
        res = (max(0, c - z)**3 - max(0, b - z)**3) / (c - b)
        res -= (max(0, b - z)**3 -max(0, a - z)**3) / (b - a)
        return res / 6.0
    
    def call_wing(self, i, z):
        assert(0 <= i < len(self.call_zs) - 1)       
        if i==0:
            a, b, c = (0.0, self.call_zs[0], self.call_zs[1])
        else:
            a, b, c = (self.call_zs[i-1], self.call_zs[i], self.call_zs[i+1])

        res = (max(0, z - a)**3 -max(0, z - b)**3) / (b - a)
        res -= (max(0, z - b)**3 -max(0, z - c)**3) / (c - b)
        return res / 6.0
    
    def convex(self, z):
        a, b, c = (-2.0, 0.0, 2.0)

        res = (max(0, z - a)**3 - max(0, z - b)**3) / (b - a)
        res -= (max(0, z - b)**3 - max(0, z - c)**3) / (c - b)
        res -= (b - a)**2
            
        res += (max(0, c - z)**3 - max(0, b - z)**3) / (c - b)
        res -= (max(0, b - z)**3 - max(0, a - z)**3) / (b - a)
        res -= (c-b)**2
        return res / 12.0  
        
    def slope(self, z):
        return z
    
    def put_wing_avg(self, i, z):
        return avg(lambda x: self.put_wing(i, x), z)
    
    def call_wing_avg(self, i, z):
        return avg(lambda x: self.call_wing(i, x), z)

    def convex_avg(self, z):
        return avg(self.convex, z)
    
    def slope_avg(self, z):
        return avg(self.slope, z)
    
    def compute_avg_vals(self, xs):
        vals = [np.array([1.0] * xs.shape[0]),
                np.vectorize(self.slope_avg)(xs), 
                np.vectorize(self.convex_avg)(xs)]
        for i in range(0, len(self.call_zs) - 1):
            vals += [np.vectorize(lambda x: self.call_wing_avg(i, x))(xs)]
        for i in range(0, len(self.put_zs) - 1):
            vals += [np.vectorize(lambda x: self.put_wing_avg(i, x))(xs)]
        return np.column_stack(vals)
    
    def compute_vals(self, xs):
        vals = [np.array([1.0] * xs.shape[0]),
                np.vectorize(self.slope)(xs), 
                np.vectorize(self.convex)(xs)]
        for i in range(0, len(self.call_zs) - 1):
            vals += [np.vectorize(lambda x: self.call_wing(i, x))(xs)]
        for i in range(0, len(self.put_zs) - 1):
            vals += [np.vectorize(lambda x: self.put_wing(i, x))(xs)]
        return np.column_stack(vals)

    def fit_coeffs(self, xs, vals, errs=None, *, smoothing=1.0e-12):        
        if errs is None:
            errs = np.array([1.0] * xs.shape[0])        

        basis_vals = self.compute_avg_vals(xs)
        basis_vals /= np.column_stack([errs])
        target = vals / errs

        var = np.matmul(basis_vals.T, basis_vals)
        var_regul = smoothing * var.trace() * np.identity(var.shape[0])
        var_regul[0,0] = 0.0  # dont need to smooth level
        var_regul[1,1] = 0.0  # dont need to smooth slope
        var_regul[2,2] = 0.0  # dont need to smooth convexity       
        var += var_regul

        cov = basis_vals.T.dot(target)
        return np.linalg.solve(var, cov)
        

