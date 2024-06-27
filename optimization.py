import numpy as np
from scipy.optimize import minimize
from scipy.integrate import solve_ivp, quad


def f_in_b(b, c):
    ''' Computes optimal reel-in factor f_in. '''

    # Coefficients of derivative of cycle power with regard to reel-in factor
    a2 = -1.5*c["f_out"] - np.cos(b)
    a1 = 2 * c["f_out"] * np.cos(b)
    a0 = 0.5 * c["force_factor_out"]/c["force_factor_in"] * c["f_out"] * (np.cos(c["elevation_angle_out"])-c["f_out"])**2 - 0.5*c["f_out"]*np.cos(b)**2
    
    # Variable subsitution to produce a depressed cubic
    p = a1 - a2**2 / 3
    q = 2/27 * a2**3 - a2*a1/3 + a0

    # Calculation of the real root, the two other roots are complex
    t = None
    if p < 0 and 4*p**3 + 27*q**2 > 0:
        t = -2 * np.sign(q) * np.sqrt(-p/3) * np.cosh(np.arccosh(-1.5*np.abs(q)/p * np.sqrt(-3/p))/3)
    elif p > 0:
        t = -2 * np.sqrt(p/3) * np.sinh(np.arcsinh(1.5*q/p * np.sqrt(3/p))/3)
    else:
        raise RuntimeError("multiple real roots")

    return t - a2/3


def dbdl(l, b, c):
    ''' Computes derivative of elevation angle with regard to tether length. '''
    f = f_in_b(*b, c)
    return (c["E_in"]*(np.cos(b)-f) - np.sin(b)) / (f*l)


def beta_in_l(c):
    ''' Integrates reel-in trajectory. '''
    b_sol = solve_ivp(dbdl, (c["tether_length_max"], c["tether_length_min"]), [c["elevation_angle_out"]], dense_output=True, args=[c]).sol
    return lambda l : b_sol(l)[0]


def reel_out_power_and_time(c):
    ''' Computes reel-out power and reel-out duration. '''
    p = c["force_factor_out"] * (np.cos(c["elevation_angle_out"])-c["f_out"])**2 * c["f_out"]
    t = (c["tether_length_max"] - c["tether_length_min"])/c["f_out"]
    return p, t


def average_reel_in_power_and_time(c):
    ''' Computes average reel-in power and reel-in duration. '''
    b_l         = beta_in_l(c)
    f_in_l      = lambda l : f_in_b(b_l(l), c)
    p_int       = lambda l : c["force_factor_in"] * (np.cos(b_l(l))-f_in_l(l))**2
    t_int       = lambda l : 1/f_in_l(l)

    En       = quad(p_int, c["tether_length_max"], c["tether_length_min"], epsrel=1e-3)[0]  # Reel-in energy
    t        = quad(t_int, c["tether_length_max"], c["tether_length_min"], epsrel=1e-3)[0]  # Reel-in duration

    return En/t, t


def negative_cycle_power(x, c):
    ''' Computes the negative cycle power. ''' 
    c["f_out"]  = x[0]

    P_out, t_out = reel_out_power_and_time(c)
    P_in_avg, t_in = average_reel_in_power_and_time(c)

    return -(P_out*t_out + P_in_avg*t_in) / (t_out + t_in)


def optimal_f_out(c):
    ''' Optimizes the reel-out factor. '''
    return minimize(negative_cycle_power, 2e-1, args=(c), method='powell', bounds=[(1e-3, np.cos(c["elevation_angle_out"]))]).x[0]
