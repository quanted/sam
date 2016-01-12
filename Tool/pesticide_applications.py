import numpy as np
from pesticide_functions import passbyval


@passbyval
def details(plant_factor, stageflag, distribflag, appnumrec, appmass, appmethod, app_windows, stagedays, cropstage):

    def application_dates():
        cumulative_maturity = np.int16((plant_factor == 1).cumsum())
        overcount = np.maximum.accumulate(cumulative_maturity * (plant_factor == 0))
        count_mature = np.concatenate(([0], (cumulative_maturity - overcount)[:-1]))
        n_plant = np.concatenate(([0], ((plant_factor == 0) & (count_mature > 0)).cumsum()))

        harvest_dates = np.where(n_plant[:-1] != n_plant[1:])[0]
        plant_dates = harvest_dates - (2 * count_mature[harvest_dates])
        emerg_dates = plant_dates + 7
        mature_dates = harvest_dates - count_mature[harvest_dates]

        return {1: plant_dates, 2: emerg_dates, 3: mature_dates, 4: harvest_dates}[cropstage]

    twindow1, pct1, twindow2, pct2 = app_windows
    appmass_init = appmass[0]
    appcount = 0

    for date in application_dates():

        # @@@ - if (cropstage, distribflag, stageflag) == (1, 2, 1): there are gaps (e.g. at 840)
        # @@@ with start and end date offsets (tp1 = 0 v tp1 = 1, etc)
        if (cropstage, distribflag) in [(1, 1), (2, 1), (3, 2), (3, 3), (4, 2), (4, 3)]: # @@@
            appcount = 0

        if distribflag == 1: # uniform application
            appmass[appcount:appcount+twindow1] = appmass_init / twindow1

        elif distribflag == 2: # step application
            appmass[appcount:appcount+twindow1] = (pct1 / 100 * appmass_init) / twindow1
            appmass[appcount+twindow1:appcount+twindow1+twindow2] = (pct2 / 100 * appmass_init) / twindow2

        elif distribflag == 3: # triangular application
            thalf = twindow1 / 2
            leg = ((appmass_init / twindow1) / thalf) * ((np.arange(thalf) + 1) ** 2 - (np.arange(thalf) ** 2))
            appmass[appcount:appcount+twindow1] = np.hstack((leg, leg[::-1]))

        full_window = twindow1 + twindow2 if distribflag == 2 else twindow1
        multiplier = 1 if stageflag == 1 else -1
        appnumrec[appcount:appcount+full_window] = np.arange(full_window) + date + (multiplier * stagedays) - 1
        appcount += full_window

    pesticide_applications = np.zeros((2, plant_factor.size))
    pesticide_applications[0][appnumrec] = appmass
    pesticide_applications[1][appnumrec] = appmethod

    return pesticide_applications


@passbyval
def process(pesticide_applications, plant_factor, rain, soil_distribution_2cm, covmax, foliar_degradation, washoff_coeff):

    def days_since(a):
        return np.arange(len(a)) - np.hstack(([0.0], np.maximum.accumulate(np.arange(len(a)) * (a > 0))[:-1]))

    applied_to_soil = pesticide_applications[0] * (pesticide_applications[1] == 1)
    applied_to_canopy = pesticide_applications[0] * (pesticide_applications[1] == 2)
    retained_by_canopy = applied_to_canopy * plant_factor * covmax
    canopy_to_soil = applied_to_canopy - retained_by_canopy

    pesticide_mass_soil = (applied_to_soil + canopy_to_soil) * soil_distribution_2cm
    degradation = np.exp(-days_since(rain + canopy_to_soil) * foliar_degradation)
    washoff = np.exp(-rain * washoff_coeff)

    canopy_to_soil_days = np.where(canopy_to_soil + rain > 0)[0]
    canopy_mass = 0.0
    for day in canopy_to_soil_days:
        canopy_mass = (canopy_mass + retained_by_canopy[day]) * degradation[day]
        pesticide_mass_soil[day] += canopy_mass - (canopy_mass * washoff[day])
        canopy_mass *= washoff[day]

    return pesticide_mass_soil


if __name__ == "__main__":
    print("This is a library. Run pesticide_calculator.py")