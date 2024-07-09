# QSM reel-in
## Implementation of a quasi-steady reel-in for power curve computation of Airborne Wind Energy.

This project extends the power computations performed in https://github.com/awecourse/resources to include a QSM model of the reel-in phase instead of an asymptotic reel-in trajectory.


### Theoretical background

The performance analysis of an airborne wind energy system is performed as described in [1], using the control strategy first described in [2], specifically implementing a constant reel-in glide ratio, thus suitable for soft kites. The simulation of the reel-in is extended using the approach described in [3].

Constraints on the tether force and the maximum reel-out generator power are implemented, leading to a 3-regime behaviour, as described in theory [2]. Additionnally, constraints on the minimum reel-in speed and the minimum reel-in generator power are implemented, minimum meaning the maximal absolute value here, the two quantities being negative. A maximum reel-out speed is not implemented, because the 3-regime strategy as described in [2] assumes that it can be adapted as needed to satify the tether force constraint in regime 2.

The optimal reel-out factor in the unconstrained wind speed regime is numerically optimized, where for each candidate reel-out factor the reel-in trajectory is computed (explained below). The reel-out factors when constrained by the tether force and/or the maximum power are computed according to the theory presented in [2], specifically, the elevation angle is increased when entering the power-limited regime.

The reel-in trajectory is integrated using the expression for a quasi-steady reel-in presented in [3, p.254], which is moved into the "tether-length" domain, in order to have clearly defined start- and endpoints, allowing for the use of a numerical ODE solver, implementing more advanced integration schemes. The optimal reel-in factor (see below) is recomputed for each reel-in elevation angle.

To solve for the reel-in factor, the general formulation of the cycle power, assuming a constant reel-in elevation angle, is maximized for each point along the trajectory. Setting the derivative equal to zero, a cubic equation with two complex roots and one real root is obtained, thus the real root corresponds to the reel-in factor that maximizes the cycle power, which thus also only has a single extremum. Roots can also be found numerically, but this approach significiantly increases computational time. For certain system parameters, more than one real root exist, these cases seem to be constrained to extreme/unrealistic inputs, f.ex. a very big reel-out factor. If more than one real root exists, a RuntimeError is thrown, and the user might need to manually investigate possible causes.

### References

[1] Van der Vlugt, R., Peschel, J., Schmehl, R.: Design and Experimental Characterization of a Pumping Kite Power System. In: Ahrens, U., Diehl, M., Schmehl, R. (eds) Airborne Wind Energy. Green Energy and Technology. Springer, Berlin, Heidelberg. Chap. 23, pp. 403-425, 2013. https://doi.org/10.1007/978-3-642-39965-7_23

[2] Luchsinger, R.H.: Pumping Cycle Kite Power. In: Ahrens, U., Diehl, M., and Schmehl, R. (eds) Airborne Wind Energy. Green Energy and Technology. Springer, Berlin, Heidelberg. Chap. 3, pp. 47-64, 2013. https://doi.org/10.1007/978-3-642-39965-7_3

[3] Fechner, U., Schmehl, R.: Model-Based Efficiency Analysis of Wind Power Conversion by a
Pumping Kite Power System. In: Ahrens, U., Diehl, M., and Schmehl, R. (eds) Airborne Wind Energy. Green Energy and Technology. Springer, Berlin, Heidelberg. Chap. 14, pp. 249-270, 2013. https://doi.org/10.1007/978-3-642-39965-7_14
