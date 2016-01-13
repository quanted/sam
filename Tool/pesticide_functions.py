import numpy as np


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
    baseflow = np.ones_like(total_runoff) * 1.e-6  # 1.e-6 is minimum baseflow
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


def pesticide_applications(plant_factor, stageflag, distribflag, appnumrec_global, appmass_global,
                           appmethod_init, app_windows, stagedays, cropstage, rain, soil_2cm, covmax, foliar_deg,
                           washoff_coeff):

    def application_dates():
        cumulative_maturity = np.int16((plant_factor == 1).cumsum())
        overcount = np.maximum.accumulate(cumulative_maturity * (plant_factor == 0))
        count_mature = np.concatenate(([0], (cumulative_maturity - overcount)[:-1]))
        n_plant = np.concatenate(([0], ((plant_factor == 0) & (count_mature > 0)).cumsum()))

        harvest_dates = np.where(n_plant[:-1] != n_plant[1:])[0]
        plant_dates = harvest_dates - (2 * count_mature[harvest_dates])
        emerg_dates = plant_dates + 7
        mature_dates = harvest_dates - count_mature[harvest_dates]

        return {1: plant_dates, 2: emerg_dates, 3: mature_dates, 4: harvest_dates}

    def details():

        twindow1, pct1, twindow2, pct2 = app_windows
        appmass_init = appmass[0]
        appcount = 0

        for date in application_dates():

            # @@@ - if (cropstage, distribflag, stageflag) == (1, 2, 1): there are gaps (e.g. at 840)
            # @@@ with start and end date offsets (tp1 = 0 v tp1 = 1, etc)
            if (cropstage, distribflag) in [(1, 1), (2, 1), (3, 2), (3, 3), (4, 2), (4, 3)]: # @@@
                appcount = 0

            if distribflag == 1:  # uniform application
                appmass[appcount:appcount+twindow1] = appmass_init / twindow1

            elif distribflag == 2:  # step application
                appmass[appcount:appcount+twindow1] = (pct1 / 100 * appmass_init) / twindow1
                appmass[appcount+twindow1:appcount+twindow1+twindow2] = (pct2 / 100 * appmass_init) / twindow2

            elif distribflag == 3:  # triangular application
                thalf = twindow1 / 2
                leg = ((appmass_init / twindow1) / thalf) * ((np.arange(thalf) + 1) ** 2 - (np.arange(thalf) ** 2))
                appmass[appcount:appcount+twindow1] = np.hstack((leg, leg[::-1]))

            full_window = twindow1 + twindow2 if distribflag == 2 else twindow1
            multiplier = 1 if stageflag == 1 else -1
            appnumrec[appcount:appcount+full_window] = np.arange(full_window) + date + (multiplier * stagedays) - 1
            appcount += full_window

        applications = np.zeros((2, plant_factor.size))
        applications[0][appnumrec] = appmass
        applications[1][appnumrec] = appmethod_init  # @@@ - how does any other application method get used?

        return applications


    def process():

        def days_since(a):
            return np.arange(len(a)) - np.hstack(([0.0], np.maximum.accumulate(np.arange(len(a)) * (a > 0))[:-1]))

        applied_to_soil = applications[0] * (applications[1] == 1)
        applied_to_canopy = applications[0] * (applications[1] == 2)
        retained_by_canopy = applied_to_canopy * plant_factor * covmax
        canopy_to_soil = applied_to_canopy - retained_by_canopy

        pesticide_mass_soil = (applied_to_soil + canopy_to_soil) * soil_2cm
        degradation = np.exp(-days_since(rain + canopy_to_soil) * foliar_deg)
        washoff = np.exp(-rain * washoff_coeff)

        canopy_to_soil_days = np.where(canopy_to_soil + rain > 0)[0]
        canopy_mass = 0.0
        for day in canopy_to_soil_days:
            canopy_mass = (canopy_mass + retained_by_canopy[day]) * degradation[day]
            pesticide_mass_soil[day] += canopy_mass - (canopy_mass * washoff[day])
            canopy_mass *= washoff[day]

        return pesticide_mass_soil

    # Create local copies of application data arrays
    appmass = np.copy(appmass_global)
    appnumrec = np.copy(appnumrec_global)

    applications = details()

    pesticide_mass_soil = process()

    return pesticide_mass_soil