from copy import deepcopy
import numpy as np


def passbyval(func):
    # Function wrapper that performs a deep copy on input arrays so that the input array is not modified
    def new(*args):
        cargs = [deepcopy(arg) for arg in args]
        return func(*cargs)
    return new


@passbyval
def transport(koc, org_carbon, bulk_density, degradation_aqueous, soil_water_m_all, delta_x, kflag, runoff,
              leaching, runoff_effic, pesticide_mass_soil):

    # Returns an array in which a value from 'a' is added and 'b' is multiplied on each iteration
    def cumulative_multiply_and_add(a, b):
        return (np.hstack((b[::-1].cumprod()[::-1], 1)) * np.hstack((0, a * b))).cumsum()[1:] / b[::-1].cumprod()[::-1]

    runoff *= runoff_effic
    leaching_dates = np.where(leaching > 0.0)[0]

    # Get retardation and deg_total
    kd = koc * org_carbon if kflag else koc
    retardation = (soil_water_m_all / delta_x) + (bulk_density * kd)
    deg_total = degradation_aqueous + ((runoff + leaching) / (delta_x * retardation))

    # Get degradation rate for each day
    degradation_rate = np.full(pesticide_mass_soil.size, np.exp(-degradation_aqueous))  # Non-leaching days
    degradation_rate[leaching_dates] =  np.exp(-deg_total[leaching_dates])  # Leaching days

    # Get total mass by accumulating pesticide_mass_soil and degrading by degradation rate
    total_mass = cumulative_multiply_and_add(pesticide_mass_soil, degradation_rate)

    runoff_mass = np.zeros_like(runoff)
    runoff_mass[leaching_dates] = \
        (((total_mass / retardation / delta_x) / deg_total) * (1 - degradation_rate) * runoff)[leaching_dates]

    # @@@ - not sure exactly the source of this problem, seems to be related to the -2 indexing in read.scenario
    return np.hstack([0.0, runoff_mass[:-1]])


def waterbody_concentration(q, xc, total_runoff, total_runoff_mass):

    # Develop hydrograph
    baseflow = np.ones_like(total_runoff) * 1.e-6
    average_runoff = np.average(total_runoff) # m3/d
    volume = (xc * 40.)
    baseflow[q >= average_runoff] = q[q >= average_runoff] - average_runoff
    total_flow = baseflow + total_runoff

    # Compute daily concentration from runoff
    conc_days = ((total_runoff_mass > 0.0) & (volume > 0.0))
    daily_concentration = np.zeros_like(total_runoff)
    daily_concentration[conc_days] = total_runoff_mass[conc_days] / volume[conc_days]

    # Modify concentration @@@ - this doesn't seem to do anything.  exp kadj always 0 or 1
    k_adj = total_flow / volume
    exp_k = np.exp(-k_adj)
    avgconc_adj = np.zeros_like(total_runoff)
    conc = 0
    for d in range(1, daily_concentration.size):
        conc += daily_concentration[d]
        avgconc_adj[d] = conc / k_adj[d] * (1 - exp_k[d])
        conc *= exp_k[d]

    # Adjust for change in units
    runoff_conc = np.nan_to_num((total_runoff_mass * 1000000.) / total_runoff)
    avgconc_adj *= 1000000.

    return total_flow, baseflow, avgconc_adj, runoff_conc

