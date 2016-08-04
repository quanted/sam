import math

import numpy as np

class ConvolutionArray:
    def __init__(self):
        self.counter = 0
        self.initialized = False

        self.array = np.array([])
        self.lookup_dict = None

    def initialize(self, n_recipes, n_dates):
        self.array = np.zeros((3, n_recipes + 1, n_dates))
        self.initialized = True

    def update(self, recipe_id, mass, runoff, baseflow):
        self.array[0, self.counter] = mass
        self.array[1, self.counter] = runoff
        self.array[2, self.counter] = baseflow
        self.counter += 1

    def write_to_file(self, outfile_path):

        np.savez_compressed(outfile_path, mass=self.array[0], runoff=self.array[1], baseflow=self.array[2])


def applications(i, scenario):
    """
    Creates a time-series array of pesticide applications and transfers those to the soil
    """

    from Tool.parameters import soil, plant

    def application_dates():

        # Predicts the dates of pesticide application based on crop stage
        cumulative_maturity = np.int16((scenario.plant_factor == 1).cumsum())
        overcount = np.maximum.accumulate(cumulative_maturity * (scenario.plant_factor == 0))
        count_mature = np.concatenate(([0], (cumulative_maturity - overcount)[:-1]))
        n_plant = np.concatenate(([0], ((scenario.plant_factor == 0) & (count_mature > 0)).cumsum()))

        harvest_dates = np.where(n_plant[:-1] != n_plant[1:])[0]
        plant_dates = harvest_dates - (2 * count_mature[harvest_dates])
        emergence_dates = plant_dates + 7
        mature_dates = harvest_dates - count_mature[harvest_dates]

        return {1: plant_dates, 2: emergence_dates, 3: mature_dates, 4: harvest_dates}[i.cropstage]

    count = 0
    application_mass = np.zeros(len(i.dates))  # Longer than it needs to be, but won't be totally filled anyway
    day_of_application = np.zeros(len(i.dates), dtype=np.int16)
    for date in application_dates():
        if i.refine == "uniform_step":  # uniform application
            application_mass[count:count + i.refine_time_window1] = i.application_rate / i.refine_time_window1
        elif i.refine == "step":  # step application
            application_mass[count:count + i.refine_time_window1] = \
                (i.refine_percent_applied1 / 100 * i.application_rate) / i.refine_time_window1
            application_mass[count + i.refine_time_window1:count + i.refine_time_window1 + i.refine_time_window2] = \
                (i.refine_percent_applied2 / 100 * i.application_rate) / i.refine_time_window2
        elif i.refine == "triangular":  # triangular application
            thalf = i.refine_time_window1 / 2
            leg = ((i.application_rate / i.refine_time_window1) / thalf) * \
                  ((np.arange(thalf) + 1) ** 2 - (np.arange(thalf) ** 2))
            application_mass[count:count + i.refine_time_window1] = np.hstack((leg, leg[::-1]))

        full_window = i.refine_time_window1 + i.refine_time_window2 if i.refine == 2 else i.refine_time_window1
        multiplier = 1 if i.stageflag == 1 else -1
        day_of_application[count:count + full_window] = np.arange(full_window) + date + (multiplier * i.stagedays) - 1
        count += full_window

    daily_mass, daily_method = np.zeros(scenario.size), np.zeros(scenario.size)
    daily_mass[day_of_application] = application_mass

    if i.application_method == 1:  # Soil application
        pesticide_mass_soil = daily_mass * soil.distrib_2cm

    elif i.application_method == 2:  # Canopy application
        retained_by_canopy = daily_mass * scenario.plant_factor * scenario.covmax
        canopy_to_soil = daily_mass - retained_by_canopy
        canopy_to_soil_days = np.where(canopy_to_soil + scenario.rain > 0)[0]
        pesticide_mass_soil = canopy_to_soil * soil.distrib_2cm
        additions = scenario.rain + canopy_to_soil
        days_since_addition = \
            np.arange(len(additions)) - \
            np.hstack(([0.0], np.maximum.accumulate(np.arange(len(additions)) * (additions > 0))[:-1]))
        degradation = np.exp(-days_since_addition * plant.foliar_degradation)
        washoff = np.exp(-scenario.rain * plant.washoff_coeff)
        canopy_mass = 0.0
        for day in canopy_to_soil_days:
            canopy_mass = (canopy_mass + retained_by_canopy[day]) * degradation[day]
            pesticide_mass_soil[day] += canopy_mass - (canopy_mass * washoff[day])
            canopy_mass *= washoff[day]

    return pesticide_mass_soil


def transport(pesticide_mass_soil, scenario, i):
    """
    Simulate transport of pesticide through runoff and erosion
    """
    from Tool.parameters import soil

    # Returns an array in which a value from 'a' is added and 'b' is multiplied on each iteration
    def cumulative_multiply_and_add(a, b):
        return (np.hstack((b[::-1].cumprod()[::-1], 1)) * np.hstack((0, a * b))).cumsum()[1:] / b[::-1].cumprod()[::-1]

    # Initialize output
    out_array = np.zeros((2, scenario.size))
    runoff_mass, erosion_mass = out_array

    runoff = scenario.runoff * soil.runoff_effic
    leach_dates = np.where(scenario.leaching > 0.0)[0]

    # Get retardation and deg_total
    retardation = (scenario.soil_water / soil.delta_x) + (scenario.bulk_density * scenario.kd)
    deg_total = i.degradation_aqueous + ((runoff + scenario.leaching) / (soil.delta_x * retardation))

    # Get degradation rate for each day
    degradation_rate = np.full(pesticide_mass_soil.size, np.exp(-i.degradation_aqueous))  # Non-leaching days
    degradation_rate[leach_dates] = np.exp(-deg_total[leach_dates])  # Leaching days

    # Get total mass by accumulating pesticide_mass_soil and degrading by degradation rate
    total_mass = cumulative_multiply_and_add(pesticide_mass_soil, degradation_rate)

    # Compute the mass of pesticide in runoff
    runoff_mass[leach_dates] = \
        (((total_mass / retardation / soil.delta_x) / deg_total) * (1 - degradation_rate) * runoff)[leach_dates]  # conc[kg/m3]*[m] = kg/m2

    # Compute the mass of pesticide from erosion
    if i.process_erosion:
        log_erosion = np.zeros(scenario.erosion.shape)
        log_erosion[scenario.erosion > 0.0] = scenario.erosion[scenario.erosion > 0.0]
        enrich = np.exp(2.0 - (0.2 * log_erosion))
        erosion_intensity = soil.erosion_effic / soil.soil_depth  # Assume uniform extraction, no decline, MMF
        enriched_eroded_mass = (scenario.erosion / 100000.) * enrich  # kg/ha -> g/cm2 (kg/ha*(ha/10000 m2)(1000 g/kg)(m2/10000 cm2) = 1/100000
        erosion = enriched_eroded_mass * scenario.kd * erosion_intensity * 10000.  # g/cm2 *(m3/g)*10000cm2/m2 -> [m]  !Enrich based on PRZM, kd: kd sediment may differ from kd erosion
        erosion_mass[leach_dates] = \
            (((total_mass / retardation / soil.delta_x) / deg_total) * (1 - degradation_rate) * erosion)[leach_dates]  # conc[kg/m3]*[m] = kg/m2

    # JCH - Need to offset by one day for some reason. Slowing things down?
    out_array = np.hstack(([[0.0], [0.0]], out_array[:,:-1]))

    return out_array.T


<<<<<<< HEAD
=======
def waterbody_concentration(q, v, length, total_runoff, runoff_mass, erosion_mass,
                            process_benthic=True, deg_photolysis=None, deg_hydrolysis=None, deg_wc = None, deg_benthic = None, koc=None):

    # MMF - we need area_wb, daily_depth, depth_0 for benthic calcs to work
>>>>>>> mfry

def waterbody_concentration(flow, runoff_and_erosion, runoff_and_erosion_mass,
                            process_benthic=True, degradation_aqueous=None, koc=None):

    def solute_holding_capacity():
        # Calculates Solute Holding capacities and mass transfer between water column and benthic regions

        # Aqueous volumes in each region
        vol1 = depth * surface_area     # total volume in water column, approximately equal to water volume alone
        vol2a = benthic.depth * surface_area  # total benthic volume
        vol2 = vol2a * benthic.porosity       # total benthic pore water volume

        # Default EXAMS conditions for partitioning
        kow = koc / .35      # DEFAULT EXAMS CONDITION ON Kow  p.35
        kpdoc1 = kow * .074  # DEFAULT RELATION IN EXAMS (LITTORAL)
        kpdoc2 = koc         # DEFAULT RELATION IN EXAMS (BENTHIC) p.16 of EXAMS 2.98 (or is it Kow*.46 ?)
        xkpb = 0.436 * kow ** .907  # DEFAULT RELATION IN EXAMS

        # mass in littoral region
        vol1a = depth[0] * surface_area  # initial volume corresponding with suspended matter reference
        m_sed_1 = water_column.sused * vol1a * .001  # SEDIMENT MASS LITTORAL
        m_bio_1 = water_column.plmas * vol1a * .001  # BIOLOGICAL MASS LITTORAL
        m_doc_1 = water_column.doc * vol1a * .001  # DOC MASS LITTORAL

        # partitioning coefficients of individual media
        kd_sed_1 = koc * water_column.froc * .001  # Kd of sediment in littoral [m3/kg]
        kd_sed_2 = koc * benthic.froc * .001  # Kd of sediment in benthic
        kd_bio = xkpb / 1000.  # Kd of biological organisms
        kd_doc_1 = kpdoc1 / 1000.  # Kd of DOC in littoral region
        kd_doc_2 = kpdoc2 / 1000.  # Kd of DOC in benthic region

        # mass in benthic region
        m_sed_2 = benthic.bulk_density * vol2a * 1000.  # as defined by EXAMS Tool.parameters m_sed_2 = BULKD/PCTWA*VOL2*100000.
        m_bio_2 = benthic.bnmas * surface_area * .001
        m_doc_2 = benthic.doc * vol2 * .001

        # solute holding capacity in region 1
        capacity_1 = kd_sed_1 * m_sed_1 + kd_bio * m_bio_1 + kd_doc_1 * m_doc_1 + vol1

        # solute holding capacity in region 2
        capacity_2 = kd_sed_2 * m_sed_2 + kd_bio * m_bio_2 + kd_doc_2 * m_doc_2 + vol2

        # Fraction going to water column and benthic
        fw1 = vol1 / capacity_1  # fw1 is daily, vol1 is daily
        fw2 = vol2 / capacity_2

        theta = capacity_2 / capacity_1

        sed_conv_factor = vol2 / fw2 / m_sed_2  # converts pore water to [Total Conc normalized to sed mass]

        # Omega mass transfer - Calculates littoral to benthic mass transfer coefficient
        omega = benthic.d_over_dx / benthic.depth  # (m3/hr)/(3600 s/hr)

        return fw1, fw2, theta, sed_conv_factor, omega

    def simultaneous_diffeq(gamma1, gamma2, omega, theta, daily_aq_peak):
        """
        ANALYTICAL SOLUTION FOR THE TWO SIMULTANEOUS DIFFERENTIAL EQNS:
                  dm1/dt = Am1 + Bm2
                  dm2/dt = Em1 + Fm2
        WITH INITIAL VALUES m1 AND m2 FOR m1 AND m2
        mn1 IS OUTPUT VALUE FOR m1 AFTER TIME T
        mn2 IS OUTPUT VALUE FOR m2 AFTER TIME T
        mavg1 IS AVERAGE VALUE OF m1 OVER TIME T
        """

        t_end = 86400.  # seconds, time step of ONE DAY
        m1, m2 = daily_aq_peak

        # Calculate constants for simultaneous_diffeq: A,B,E,F
        # This reduces the model equivalent parameters to the coefficients needed for solving simultaneous_diffeq
        a = -gamma1 - omega * theta
        b = omega * theta
        e = omega
        f = -gamma2 - omega

        af = a + f
        dif = 4 * ((f * a) - (b * e))
        bbb = np.sqrt(af * af - dif)

        root1 = (af + bbb) / 2.
        root2 = (af - bbb) / 2.

        dd = (root1 - a) / b
        ee = (root2 - a) / b
        ff = ee - dd
        x1 = (ee * m1 - m2) / ff
        y1 = (m2 - dd * m1) / ff

        # Calculate new concentrations for next step
        rt1 = root1 * t_end
        rt2 = root2 * t_end
        exrt1 = np.exp(rt1)
        exrt2 = np.exp(rt2)
        ccc = x1 * exrt1
        ddd = y1 * exrt2

        #values for m1 and m2 after time step t_end
        mn = np.array([(ccc + ddd),              # Water column
                       (dd * ccc + ee * ddd)])   # Benthic

        # AVERAGE DAILY CONCENTRATION SOLUTION: set up for daily average, but can be changed by adjusting time step
        gx = x1 / root1
        hx = y1 / root2

        term1 = gx * exrt1  # term3 = -X1/root1*exp(root1*T1)
        term2 = hx * exrt2  # term4 = -Y1/root2*exp(root2*T1
        term3 = -gx
        term4 = -hx

        mavg = np.array([(term1 + term2 + term3 + term4),                       # Water column
                         (term1 * dd + term2 * ee + term3 * dd + term4 * ee)])  # Benthic

        mavg /= t_end

        return mn, mavg

    from Tool.parameters import benthic, soil, water_column, stream_channel

    total_runoff, total_erosion = runoff_and_erosion
    runoff_mass, erosion_mass = runoff_and_erosion_mass
    n_dates = total_runoff.size

    # Develop hydrograph
    average_runoff = np.average(total_runoff)  # m3/d
    baseflow = np.ones_like(flow.q) * 1.e-6  # 1.e-6 is minimum baseflow
    baseflow[flow.q >= average_runoff] = flow.q[flow.q >= average_runoff] - average_runoff
    total_flow = baseflow + total_runoff

    # Compute stream channel properties
    cross_section = total_flow / flow.v
    volume = cross_section * 40.

    # Compute daily concentration from runoff
    conc_days = ((runoff_mass > 0.0) & (volume > 0.0))
    daily_concentration = np.zeros(n_dates)
    daily_concentration[conc_days] = runoff_mass[conc_days] / volume[conc_days]

    #Overall littoral deg rate (k_adj). Note: deg_wc = aerobic aquatic metabolism rate, deg_benthic = aerobic sorbed metabolism rate
    k_adj = (total_flow / volume) + (deg_photolysis + deg_hydrolysis)*fw1 + (deg_wc*fw1) + deg_benthic*(1-fw1)
    exp_k = np.exp(-k_adj)
    aqconc_avg_wb = np.zeros_like(total_runoff)

    #Overall benthic deg rate (k_adj2). Note: aerobic & anaerobic metabolism rates are equated, so deg_wc = anaer aq met rate, deg_benthic = anaer sorb met rate
    k_adj2 = deg_wc * fw2 + deg_benthic * (1 - fw2) + deg_hydrolysis * fw2

    if process_benthic:

        width = stream_channel.a * np.power(cross_section, stream_channel.b)
        depth = cross_section / width
        surface_area = width * flow.l

        # Compute benthic solute holding capacity
        fw1, fw2, theta, sed_conv_factor, omega = solute_holding_capacity()

        mass_input = np.vstack([runoff_mass + ((1. - soil.prben) * erosion_mass),  # Water Column
                                soil.prben * erosion_mass]).T                      # Benthic

        # Beginning day aquatic concentrations, considered Peak Aqueous Daily Conc in Water Column
        daily_peak = np.zeros((n_dates, 2))
        daily_avg = np.zeros((n_dates, 2))
        fw = np.array([fw1, fw2]).T

    else:
        daily_avg = daily_peak = np.zeros([1, 2])

    # Reset starting values
    aqconc_wb = 0
    mn = np.array([0, 0])

    for d in range(daily_concentration.size):

        # Compute daily average concentration in the water body - when no Benthic layer considered
        aqconc_wb += daily_concentration[d]    #initial water body concentration for current time step

        # Daily avg aq conc in water body, area under curve/t = Ci/k*(1-e^-k), NO benthic
        aqconc_avg_wb[d] = aqconc_wb / k_adj[d] * (1 - exp_k[d])

        #initial water body concentration for next time step
        aqconc_wb *= exp_k[d]

        if process_benthic:  # Addition of Benthic region, assume VVWM benthic properties

            # Add mass input to antecedent mass
            daily_mass = mn + mass_input[d]

            # Get daily volume of water column and benthic zones  # JCH - different from volume because vol 40.  why?
            daily_vol = np.array([(depth[d] * surface_area[d]),                             # Water column
                                  (benthic.depth * surface_area[d] * benthic.porosity)])    # Benthic zone

            # Convert to aqueous concentrations (peak) at beginning of day
            daily_peak[d] = daily_mass * fw[d] / daily_vol

<<<<<<< HEAD
            # Solve simultaneous differential equation
            new_aqconc, daily_avg[d] = \
                simultaneous_diffeq(k_adj[d], degradation_aqueous, omega, theta[d], daily_peak[d])
=======
            #For simul diffeq soln: mn1,mn2,mavg1,mavg2 = new_aqconc1, new_aqconc2, aqconc_avg1[d], aqconc_avg2[d]
            #Note: aqconc_avg1 and aqconc_avg2 are outputted - Daily avg aq conc in WC and Benthic regions
            new_aqconc1, new_aqconc2, aqconc_avg1[d], aqconc_avg2[d] = \
                simultaneous_diffeq(k_adj[d], k_adj2[d], omega, theta, aq_conc1[d], aq_conc2[d])
>>>>>>> mfry

            # Masses m1 and m2 after time step, t_end
            mn = new_aqconc / fw[d] * daily_vol

            # Daily average mass in WC - not currently outputted, but calculated in VVWM
            #mavg1_store = aqconc_avg[0][d] / fw1[d] * depth[d] * surface_area[d]

    with np.errstate(divide='ignore', invalid='ignore'):
        runoff_conc = np.nan_to_num(runoff_mass / total_runoff)
        erosion_conc = np.nan_to_num(erosion_mass / total_runoff)

    # Adjust for change in units - for concentrations: 1 kg/m3 to 1000000. ug/L
    conv = 1000000.  # kg/m3 -> ug/L
    aqconc_avg_wb, daily_avg, daily_peak, runoff_conc = \
        map(lambda x: x * conv, (aqconc_avg_wb, daily_avg, daily_peak, runoff_conc))

    # Diagnostics
    aqconc_avg_wb, daily_avg, daily_peak = \
        tuple(map(lambda x: x.T, (aqconc_avg_wb, daily_avg, daily_peak)))

    return total_flow, baseflow, runoff_conc, aqconc_avg_wb, daily_avg, daily_peak

if __name__ == "__main__":
    import pesticide_calculator
    pesticide_calculator.main()
