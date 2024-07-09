import numpy as np
import matplotlib.pyplot as plt

from phases import power_curve, reel_in_trajectory


def plot_power_curve(velocity_range, inputs):
    ''' Plotting the power curve for all velocities in velocity_range, with the system characteristics as given in inputs. '''
    velocity_bounds, reeling_factor_out, reel_out_angles, power_out, power_in, cycle_power, time_out, time_in = power_curve(velocity_range, inputs)

    # Plot power curve
    plt.figure(dpi=300)
    plt.plot(velocity_range, cycle_power*1e-3, label=r"$P_c$")
    plt.plot(velocity_range, power_out*1e-3, label=r"$P_o$")
    plt.plot(velocity_range, -power_in*1e-3, label="$-P_i$")
    plt.vlines(velocity_bounds, 0, 25, colors='black', linestyles='dotted')
    plt.ylim(0, 25)
    plt.xlabel("Wind speed (m/s)")
    plt.ylabel("Mechanical power (kW)")
    plt.legend(framealpha=0)
    plt.savefig("Power_curve.svg")

    # Plotting other computed quantities has to be implemented manually if desired
    pass


def plot_single_reel_in_trajectory(vw, inputs):
    ''' Plotting the trajectory and the evolution of the power and reeling factor for a single reel-in at one wind speed, with the system characteristics as given in inputs. '''

    b, f, p = reel_in_trajectory(np.array([vw]), inputs.copy())

    # Plot reel-in trajectory
    plt.figure(dpi=300)
    plt.plot(inputs["l_range"]*np.cos(b[0]), inputs["l_range"]*np.sin(b[0]), "k", label="Reel-out")
    plt.plot(inputs["l_range"]*np.cos(b), inputs["l_range"]*np.sin(b), label="Reel-in")
    plt.plot(np.linspace(-inputs["tether_length_min"], inputs["tether_length_min"]), np.sqrt(inputs["tether_length_min"]**2 - np.linspace(-inputs["tether_length_min"], inputs["tether_length_min"])**2), "k:", label="Minimum tether length")
    plt.axis("scaled")
    plt.title(f"Wind speed: {vw} m/s")
    plt.xlabel("Horizontal distance (m)")
    plt.ylabel("Height (m)")
    plt.legend(framealpha=0)
    plt.savefig("Reel-in_trajectory.svg")

    # Plot reel-in power & factor
    fig, ax1 = plt.subplots(dpi=300)
    ax1.plot(inputs["l_range"], p*1e-3, label=r"$P_i$")
    ax1.yaxis.set_inverted(True)
    ax1.xaxis.set_inverted(True)
    ax1.set(xlabel="Tether length (m)", ylabel="Reel-in power (kW)", title=f"Wind speed: {vw} m/s")
    ax2 = ax1.twinx()
    ax2.plot(inputs["l_range"], f, "k", label=r"$f_i$")
    ax2.yaxis.set_inverted(True)
    ax2.xaxis.set_inverted(True)
    ax2.set(ylabel="Reel-in factor (-)")
    fig.legend(loc="upper right", bbox_to_anchor=(1, 0.8), bbox_transform=ax1.transAxes, framealpha=0)
    fig.tight_layout()
    fig.savefig("Reel-in_power_and_factor.svg")
    pass


# Environmental properties
atmosphere_density          =  1.225      # kg/m**3

# Kite properties
kite_planform_area          =  16.7       # m**2
kite_lift_coefficient_out   =  1.0        # -
kite_drag_coefficient_out   =  0.2        # -
kite_lift_coefficient_in    =  0.14       # -
kite_drag_coefficient_in    =  0.07       # -

# Tether properties
nominal_tether_force        =  5000       # N
tether_drag_coefficient     =  1.1        # -
tether_diameter             =  0.00484    # m
tether_length_max           =  375        # m
tether_length_min           =  200        # m

# Generator properties
nominal_generator_power_out =  20000      # W
nominal_generator_power_in  =  -20000      # W

# Operational parameters
elevation_angle_out         =  25         # deg
min_reel_in_speed           =  -8         # m/s

# Computation parameters
tether_length_samples       = 2**8 + 1

# Derived properties
rm_out = 0.5 * (tether_length_min + tether_length_max)
CD_out = kite_drag_coefficient_out + 0.25 * tether_drag_coefficient \
         * tether_diameter * rm_out / kite_planform_area
E2out  = (kite_lift_coefficient_out / CD_out)**2
E_in   = (kite_lift_coefficient_in  / kite_drag_coefficient_in)
E2in   = E_in**2
force_factor_out = kite_lift_coefficient_out * np.sqrt(1+1/E2out) * (1+E2out)
force_factor_in  = kite_lift_coefficient_in  * np.sqrt(1+1/E2in)  * (1+E2in)
elevation_angle_out = np.deg2rad(elevation_angle_out)
l_range = np.linspace(tether_length_max, tether_length_min, tether_length_samples)
sample_spacing = (tether_length_max - tether_length_min)/(tether_length_samples-1)


# all inputs / constants
inputs = {"elevation_angle_out" : elevation_angle_out, "force_factor_out" : force_factor_out, "force_factor_in" : force_factor_in, "E_in" : E_in, "tether_length_max" : tether_length_max, "tether_length_min" : tether_length_min, "nominal_tether_force" : nominal_tether_force, "kite_planform_area" : kite_planform_area, "atmosphere_density" : atmosphere_density, "nominal_generator_power_out" : nominal_generator_power_out, "nominal_generator_power_in" : nominal_generator_power_in, "min_reel_in_speed" : min_reel_in_speed, "wind_speed": None, "l_range" : l_range, "sample_spacing" : sample_spacing}


# Power curve
velocity_range = np.linspace(1, 20, 191)  # Wind speed range of interest (m/s)
plot_power_curve(velocity_range, inputs.copy())

# Reel-in trajectory
wind_speed = 10.  # m/s
plot_single_reel_in_trajectory(wind_speed, inputs.copy())
