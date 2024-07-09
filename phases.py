import numpy as np

from optimization import *


def phase_velocities(c):
    ''' Computes the wind speeds of the boundaries between the operational regimes. '''
    # regime 1 -> regime 2
    v_nf = np.sqrt(2. * c["nominal_tether_force"] / (c["atmosphere_density"] * c["kite_planform_area"] * c["force_factor_out"] * (np.cos(c["elevation_angle_out"]) - c["f_out_unconstrained"])**2))

    # regime 2 -> regime 3
    v_np = v_nf * (1 + (c["nominal_generator_power_out"]/(c["nominal_tether_force"] * v_nf) - c["f_out_unconstrained"]) / np.cos(c["elevation_angle_out"]))

    return v_nf, v_np

def f_out_vw(v_w, c):
    ''' Computes the reel-out factor for all wind speeds. '''
    v_nf, v_np = phase_velocities(c)
    
    # regime 1
    result = np.full_like(v_w, c["f_out_unconstrained"])

    # regime 2
    result[v_w > v_nf] = v_nf/v_w[v_w > v_nf] * (np.cos(c["elevation_angle_out"]) * (v_w[v_w > v_nf]/v_nf - 1) + c["f_out_unconstrained"])
    
    # regime 3
    f_out_np = v_nf/v_np * (np.cos(c["elevation_angle_out"]) * (v_np/v_nf - 1) + c["f_out_unconstrained"])
    result[v_w > v_np] = f_out_np * v_np/v_w[v_w > v_np]
    return result

def adapted_elevation_angle_out(v_w, c):
    ''' Computes the reel-out elevation angle for all wind speeds. Reel-out elevation angle is increased in regime 3 to decrease the efficiency, thus keeping the power below/at the limit. '''
    _, v_np = phase_velocities(c)

    # regime 1 & 2
    result = np.full_like(v_w, c["elevation_angle_out"])
    
    # regime 3
    result[v_w > v_np] = np.arccos(np.sqrt(2 * c["nominal_tether_force"] / (c["atmosphere_density"] * v_w[v_w > v_np]**2 * c["kite_planform_area"] * c["force_factor_out"])) + f_out_vw(v_w, c)[v_w > v_np])
    return result


def power_curve(v_w, c):
    ''' Computes the power curve and other relevant quantities for all wind speeds. '''
    c["f_out_unconstrained"] = optimal_f_out(c)
    velocity_bounds = phase_velocities(c)

    reeling_factor_out = f_out_vw(v_w, c)
    reel_out_angles = adapted_elevation_angle_out(v_w, c)
    power_out = []
    power_in = []
    cycle_power = []
    time_out = []
    time_in = []

    
    for v, f, b in zip(v_w, reeling_factor_out, reel_out_angles):
        c["elevation_angle_out"] = b
        c["f_out"] = f
        c["wind_speed"] = v

        p_in_avg, t_in = average_reel_in_power_and_time(c)
        power_in.append(p_in_avg * c["kite_planform_area"] * 0.5 * c["atmosphere_density"] * v**3)
        time_in.append(t_in / v)

        p_out, t_out = reel_out_power_and_time(c)
        power_out.append(p_out * c["kite_planform_area"] * 0.5 * c["atmosphere_density"] * v**3)
        time_out.append(t_out / v)

        cycle_power.append((p_out*t_out + p_in_avg*t_in) / (t_out + t_in) * c["kite_planform_area"] * 0.5 * c["atmosphere_density"] * v**3)
        pass
    
    return velocity_bounds, reeling_factor_out, reel_out_angles, np.array(power_out), np.array(power_in), np.array(cycle_power), np.array(time_out), np.array(time_in)


def reel_in_trajectory(vw, c):
    ''' Computes relevant quantities for the reel-in trajectory at wind speed vw. '''
    c["f_out_unconstrained"] = optimal_f_out(c)
    c["f_out"] = f_out_vw(vw, c)[0]
    c["elevation_angle_out"] = adapted_elevation_angle_out(vw, c)[0]
    c["wind_speed"] = vw[0]

    b = beta_in_l(c)
    f = f_in_b(b, c)
    p = -c["force_factor_in"] * (np.cos(b)-f)**2 * c["kite_planform_area"] * 0.5 * c["atmosphere_density"] * vw**3

    return b, f, p