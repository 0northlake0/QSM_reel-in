import numpy as np
from scipy.optimize import minimize
from scipy.integrate import solve_ivp, romb


def f_in_b(b, c):
    ''' Computes optimal reel-in factor f_in. '''
    # Coefficients of derivative of cycle power with regard to reel-in factor
    # a3 = 1.0
    a2 = -1.5*c["f_out"] - np.cos(b)
    a1 = 2 * c["f_out"] * np.cos(b)
    a0 = 0.5 * c["force_factor_out"]/c["force_factor_in"] * c["f_out"] * (np.cos(c["elevation_angle_out"])-c["f_out"])**2 - 0.5*c["f_out"]*np.cos(b)**2
    
    # Variable subsitution to produce a depressed cubic
    p = a1 - a2**2 / 3
    q = 2/27 * a2**3 - a2*a1/3 + a0

    # Check for the number of real roots
    if np.any(4*p**3 + 27*q**2 < 0):
        raise RuntimeError("multiple real roots")  # usually caused by unrealistic input values, f.ex f_out too small or big.
    else:
        pass
    
    # Calculation of the real root, the two other roots are complex
    p1, q1 = p[p<0], q[p<0]
    p2, q2 = p[p>0], q[p>0]

    t = np.zeros_like(b)
    t[p<0] = -2 * np.sign(q1) * np.sqrt(-p1/3) * np.cosh(np.arccosh(-1.5*np.abs(q1)/p1 * np.sqrt(-3/p1))/3)  # Case 1
    t[p>0] = -2 * np.sqrt(p2/3) * np.sinh(np.arcsinh(1.5*q2/p2 * np.sqrt(3/p2))/3)  # Case 2

    # Variable resubstitution
    f = t - a2/3

    # Correct for maximum reel-in speed
    f = np.maximum(f, c["min_reel_in_speed"]/c["wind_speed"])
    # Correct for maximum reel-in power
    f = np.maximum(f, np.cos(b) - np.sqrt(-c["nominal_generator_power_in"] / (c["force_factor_in"] * c["kite_planform_area"] * 0.5 * c["atmosphere_density"] * c["wind_speed"]**3)))

    return f


def dbdl(l, b, c):
    ''' Computes derivative of elevation angle with regard to tether length. '''
    f = f_in_b(b, c)
    return (c["E_in"]*(np.cos(b)-f) - np.sin(b)) / (f*l)


def beta_in_l(c):
    ''' Integrates reel-in trajectory. '''
    return solve_ivp(dbdl, (c["tether_length_max"], c["tether_length_min"]), [c["elevation_angle_out"]], t_eval=c["l_range"], args=[c]).y[0]


def reel_out_power_and_time(c):
    ''' Computes reel-out power and reel-out duration. '''
    p = c["force_factor_out"] * (np.cos(c["elevation_angle_out"])-c["f_out"])**2 * c["f_out"]
    t = (c["tether_length_max"] - c["tether_length_min"])/c["f_out"]
    return p, t


def average_reel_in_power_and_time(c):
    ''' Computes average reel-in power and reel-in duration. '''
    b  = beta_in_l(c)
    f  = f_in_b(b, c)

    En = romb(-c["force_factor_in"] * (np.cos(b)-f)**2, dx=c["sample_spacing"])  # Reel-in energy
    t  = romb(-1./f, dx=c["sample_spacing"])  # Reel-in duration

    return En/t, t


def objective_function(x, c):
    ''' Computes the negative cycle power. ''' 
    c["f_out"]  = x[0]

    P_out, t_out = reel_out_power_and_time(c)
    P_in_avg, t_in = average_reel_in_power_and_time(c)

    return -(P_out*t_out + P_in_avg*t_in) / (t_out + t_in)


def optimal_f_out(c):
    ''' Optimizes the reel-out factor. '''
    return minimize(objective_function, 2e-1, args=(c), bounds=[(1e-2, 0.7)]).x[0]
