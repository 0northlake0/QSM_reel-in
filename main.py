import numpy as np
import matplotlib.pyplot as plt

from phases import power_curve


# Environmental properties
atmosphere_density        =  1.225      # kg/**3

# Kite properties
kite_planform_area        =  16.7       # m**2
kite_lift_coefficient_out =  1.0        # -
kite_drag_coefficient_out =  0.2        # -
kite_lift_coefficient_in  =  0.14       # -
kite_drag_coefficient_in  =  0.07       # -

# Tether properties
nominal_tether_force      =  5000       # N
tether_drag_coefficient   =  1.1        # -
tether_diameter           =  0.00484    # m
tether_length_max         =  375        # m
tether_length_min         =  200        # m

# Generator properties
nominal_generator_power   =  20000      # W

# Operational parameters
elevation_angle_out       =  25         # deg

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

# all inputs / constants
inputs = {"elevation_angle_out" : elevation_angle_out, "force_factor_out" : force_factor_out, "force_factor_in" : force_factor_in, "E_in" : E_in, "tether_length_max" : tether_length_max, "tether_length_min" : tether_length_min, "nominal_tether_force" : nominal_tether_force, "kite_planform_area" : kite_planform_area, "atmosphere_density" : atmosphere_density, "nominal_generator_power" : nominal_generator_power}


# Wind speed range of interest
velocity_range = np.linspace(1, 20, 1901)

# Computation
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
plt.legend()
plt.savefig("Power_curve.svg")
