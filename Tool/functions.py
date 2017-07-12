import math
import numpy as np

from numba import njit, guvectorize


@njit
def canopy_loop(n_dates, application_mass, plant_factor, covmax, soil_2cm, foliar_degradation, rain, washoff_coeff):
    canopy_mass = 0
    canopy_to_soil = np.zeros(n_dates)  # Cumulative
    last_application = 0
    for day in range(n_dates):
        if application_mass[day] > 0:
            canopy_pesticide_additions = application_mass[day] * plant_factor[day] * covmax
            canopy_to_soil[day] = (application_mass[day] - canopy_pesticide_additions) * soil_2cm
            canopy_mass = canopy_pesticide_additions + canopy_mass * np.exp(
                (day - last_application) * foliar_degradation)
            last_application = day
        if rain[day] > 0:
            canopy_mass *= np.exp((day - last_application) * foliar_degradation)
            pesticide_remaining = canopy_mass * np.exp(-rain[day] * washoff_coeff)
            canopy_to_soil[day] += canopy_mass - pesticide_remaining
            last_application = day
    return canopy_to_soil


@njit
def concentration_loop(n_dates, daily_concentration, k_adj, daily_volume, mass_input, fw1, fw2, omega, theta, deg_aq):
    # Beginning day aquatic concentrations, considered Peak Aqueous Daily Conc in Water Column
    daily_peak = np.zeros((2, n_dates))
    daily_avg = np.zeros((2, n_dates))
    aqconc_avg_wb = np.zeros(n_dates)

    # Reset starting values
    exp_k = np.exp(-k_adj)
    aqconc_wb = 0
    antecedent_mass = np.zeros(2)  # mn

    for day in range(daily_concentration.size):
        # Add mass input to antecedent mass
        daily_mass = antecedent_mass + mass_input[day]

        # Convert to aqueous concentrations (peak) at beginning of day
        # JCH - fw comes from solute_holding_capacity. Fraction going into each section. Should fw[0] + fw[1] = 1?
        daily_peak[0, day] = daily_mass[0] * fw1[day] / daily_volume[day, 0]
        daily_peak[1, day] = daily_mass[1] * fw2[day] / daily_volume[day, 1]

        # Compute daily average concentration in the water body - when no Benthic layer considered
        aqconc_wb += daily_concentration[day]  # initial water body concentration for current time step

        # Daily avg aq conc in water body, area under curve/t = Ci/k*(1-e^-k), NO benthic
        aqconc_avg_wb[day] = aqconc_wb / k_adj[day] * (1 - exp_k[day])

        # initial water body concentration for next time step
        aqconc_wb *= exp_k[day]

        # For simul diffeq soln: mn1,mn2,mavg1,mavg2 = new_aqconc1, new_aqconc2, aqconc_avg1[d], aqconc_avg2[d]
        # Note: aqconc_avg1 and aqconc_avg2 are outputted - Daily avg aq conc in WC and Benthic regions
        new_aqconc, wc_avg, benthic_avg = simultaneous_diffeq(k_adj[day], deg_aq, omega, theta[day], daily_peak[:, day])
        daily_avg[0, day] = wc_avg
        daily_avg[1, day] = benthic_avg

        # Masses m1 and m2 after time step, t_end
        antecedent_mass[0] = new_aqconc[0] / fw1[day] * daily_volume[day, 0]
        antecedent_mass[1] = new_aqconc[1] / fw2[day] * daily_volume[day, 1]

    return aqconc_avg_wb, daily_avg, daily_peak


@guvectorize(['void(float64[:], float64[:], float64[:])'], '(n),(o)->(n)', nopython=True)
def cumulative_multiply_and_add(a, b, res):
    res[0] = a[0]
    for i in range(1, a.size):
        res[i] = res[i - 1] * b[i - 1] + a[i]


def impulse_response_function(alpha, beta, length):
    def gamma_distribution(t, a, b):
        a, b = map(float, (a, b))
        tau = a * b
        return ((t ** (a - 1)) / (((tau / a) ** a) * math.gamma(a))) * math.exp(-(a / tau) * t)

    return np.array([gamma_distribution(i, alpha, beta) for i in range(length)])


@guvectorize(['void(float64[:], int16[:], int16[:], int16[:], int16[:], float64[:])'], '(p),(o),(o),(p),(n)->(o)',
             nopython=True)
def moving_window(time_series, window_sizes, endpoints, years_since_start, year_sizes, res):
    # Count the number of times the concentration exceeds the test threshold in each year
    counts = np.zeros((year_sizes.size, window_sizes.size))
    for test_number in range(window_sizes.size):
        window_size = window_sizes[test_number]
        threshold = endpoints[test_number]
        window_sum = np.sum(time_series[:window_size])
        for day in range(window_size, len(time_series)):
            year = years_since_start[day]
            window_sum += time_series[day] - time_series[day - window_size]
            avg = window_sum / window_size
            if avg > threshold:
                counts[year, test_number] += 1

    # Average the number of yearly exceedances for each test
    for test_number in range(window_sizes.size):
        exceedance = 0
        for year in range(year_sizes.size):
            exceedance += counts[year, test_number] / year_sizes[year]
        res[test_number] = exceedance / year_sizes.size


@njit
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

    # values for m1 and m2 after time step t_end
    mn = np.zeros(2)
    mn[0] = ccc + ddd  # Water column
    mn[1] = dd * ccc + ee * ddd  # Benthic

    # AVERAGE DAILY CONCENTRATION SOLUTION: set up for daily average, but can be changed by adjusting time step
    gx = x1 / root1
    hx = y1 / root2

    term1 = gx * exrt1  # term3 = -X1/root1*exp(root1*T1)
    term2 = hx * exrt2  # term4 = -Y1/root2*exp(root2*T1
    term3 = -gx
    term4 = -hx

    mavg_wc = (term1 + term2 + term3 + term4) / t_end  # Water column
    mavg_ben = (term1 * dd + term2 * ee + term3 * dd + term4 * ee) / t_end  # Benthic

    return mn, mavg_wc, mavg_ben


def solute_holding_capacity(depth, surface_area, koc):
    """Calculates Solute Holding capacities and mass transfer between water column and benthic regions"""

    from .parameters import benthic, water_column

    # Aqueous volumes in each region
    vol1 = depth * surface_area  # total volume in water column, approximately equal to water volume alone
    vol2a = benthic.depth * surface_area  # total benthic volume
    vol2 = vol2a * benthic.porosity  # total benthic pore water volume

    # Default EXAMS conditions for partitioning
    kow = koc / .35  # DEFAULT EXAMS CONDITION ON Kow  p.35
    kpdoc1 = kow * .074  # DEFAULT RELATION IN EXAMS (LITTORAL)
    kpdoc2 = koc  # DEFAULT RELATION IN EXAMS (BENTHIC) p.16 of EXAMS 2.98 (or is it Kow*.46 ?)
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
    m_sed_2 = benthic.bulk_density * vol2a * 1000.  # as defined by EXAMS parameters m_sed_2 = BULKD/PCTWA*VOL2*100000.
    m_bio_2 = benthic.bnmas * surface_area * .001
    m_doc_2 = benthic.doc * vol2 * .001

    # solute holding capacity in regions 1 and 2
    capacity_1 = kd_sed_1 * m_sed_1 + kd_bio * m_bio_1 + kd_doc_1 * m_doc_1 + vol1
    capacity_2 = kd_sed_2 * m_sed_2 + kd_bio * m_bio_2 + kd_doc_2 * m_doc_2 + vol2

    # Fraction going to water column and benthic
    fw1 = vol1 / capacity_1  # fw1 is daily, vol1 is daily
    fw2 = vol2 / capacity_2

    theta = capacity_2 / capacity_1

    sed_conv_factor = vol2 / fw2 / m_sed_2  # converts pore water to [Total Conc normalized to sed mass]

    # Omega mass transfer - Calculates littoral to benthic mass transfer coefficient
    omega = benthic.d_over_dx / benthic.depth  # (m3/hr)/(3600 s/hr)

    return fw1, fw2, theta, sed_conv_factor, omega


"""Scenario stuff"""


def compute_concentration(transported_mass, runoff, n_dates, q):
    """ Concentration function for time of travel """

    """
    JCH - the idea here is that baseflow is estimated as the difference between total predicted q (erom) and
    mean modeled runoff. Predicted baseflow isn't event specific and is not sensitive to runoff itself apart
    from the mean.  This balances the sum of the modeled Q with the sum of the predicted Q
    """
    mean_runoff = runoff.mean()  # m3/d
    baseflow = np.subtract(q, mean_runoff, out=np.zeros(n_dates), where=(q > mean_runoff))
    total_flow = runoff + baseflow
    concentration = np.divide(transported_mass, total_flow, out=np.zeros(n_dates),
                              where=(total_flow != 0))
    runoff_concentration = np.divide(transported_mass, runoff, out=np.zeros(n_dates),
                                     where=(runoff != 0))

    return total_flow, map(lambda x: x * 1000000., (concentration, runoff_concentration))  # kg/m3 -> ug/L


def partition_benthic(reach, runoff, runoff_mass, erosion_mass):
    from .parameters import soil, stream_channel, benthic

    try:
        reach = self.region.flow_file.fetch(reach)
        q, v, l = reach.q, reach.v, reach.l
    except AttributeError:
        return None, None, (None, None)

    mean_runoff = runoff.mean()  # m3/d
    baseflow = np.subtract(q, mean_runoff, out=np.zeros(self.i.n_dates), where=(q > mean_runoff))
    total_flow = runoff + baseflow
    mixing_cell = 40.  # meters
    cross_section = total_flow / v
    width = stream_channel.a * np.power(cross_section, stream_channel.b)
    depth = cross_section / width
    surface_area = width * l
    volume = np.array([(depth * surface_area),  # Water column
                       (benthic.depth * surface_area * benthic.porosity)])  # Benthic zone

    # Compute concentration in runoff of runoff mass and erosion mass
    runoff_conc = np.divide(runoff_mass, runoff, out=np.zeros(self.i.n_dates), where=(runoff != 0))
    daily_conc = np.divide(runoff_mass + erosion_mass, mixing_cell, out=np.zeros(self.i.n_dates),
                           where=(runoff_mass + erosion_mass > 0.0) & (mixing_cell > 0.0))

    # Divide mass loading between water column and benthic zones
    mass_input = np.vstack([runoff_mass + ((1. - soil.prben) * erosion_mass),  # Water Column
                            soil.prben * erosion_mass]).T  # Benthic
    # Partition concentration into benthic and water column concentrations
    # This needs to be cleaned up
    # Compute benthic solute holding capacity
    fw1, fw2, theta, sed_conv_factor, omega = solute_holding_capacity(depth, surface_area, self.i.koc)

    k_adj = np.array((total_flow / mixing_cell) + (self.i.deg_photolysis + self.i.deg_hydrolysis) * fw1 + \
                     (self.i.deg_wc * fw1) + self.i.deg_benthic * (1 - fw1))

    aqconc_avg_wb, daily_avg, daily_peak = \
        concentration_loop(self.i.n_dates, daily_conc, k_adj, volume,
                              mass_input, fw1, fw2, omega, theta, self.i.deg_aqueous)

    return map(lambda x: x * 1000000., (runoff_conc, aqconc_avg_wb, daily_avg, daily_peak))


def pesticide_applications(n_dates, applications, new_year, active_crops, event_dates, plant_factor, rain, covmax):
    from .parameters import soil, plant

    scenario_applications = set(filter(lambda x: x.crop in active_crops, applications))
    application_mass = np.zeros((2, n_dates))
    canopy_applications = False

    # Determine how much pesticide is applied and when
    for app in scenario_applications:
        index = ['groud', 'foliar'].index(app.method)
        if index:
            canopy_applications = True
        start_dates = np.int16(new_year + event_dates[app.event] + app.offset)
        first_window = \
            np.repeat(start_dates, app.window1) + np.tile(np.arange(app.window1), len(start_dates))
        application_mass[index, first_window] += (app.rate * (app.pct1 / 100.)) / app.window1
        if app.refine == 'step':
            second_window = \
                np.repeat(start_dates + app.window1, app.window2) + \
                np.tile(np.arange(app.window2), len(start_dates))
            application_mass[index, second_window] += (app.rate * (app.pct2 / 100.)) / app.window2

    pesticide_mass_soil = application_mass[0] * soil.distrib_2cm
    if canopy_applications:
        pesticide_mass_soil += \
            canopy_loop(n_dates, application_mass[1], np.array(plant_factor), covmax,
                           soil.distrib_2cm, plant.foliar_degradation, np.array(rain), plant.washoff_coeff)

    return pesticide_mass_soil


def pesticide_transport(n_dates, pesticide_mass_soil, runoff, erosion, leaching, org_carbon, bulk_density, soil_water,
                        koc, kd_flag, deg_aqueous):
    """ Simulate transport of pesticide through runoff and erosion """
    from .parameters import soil

    # Initialize output
    runoff_mass = np.zeros(n_dates)
    erosion_mass = np.zeros(n_dates)
    runoff = runoff * soil.runoff_effic
    leach_dates = np.where(leaching > 0.0)[0]

    kd = koc * org_carbon if kd_flag else koc

    # Get retardation and deg_total
    retardation = (soil_water / soil.delta_x) + (bulk_density * kd)
    deg_total = deg_aqueous + ((runoff + leaching) / (soil.delta_x * retardation))

    # Get degradation rate for each day
    degradation_rate = np.full(pesticide_mass_soil.size, np.exp(-deg_aqueous))  # Non-leaching days
    degradation_rate[leach_dates] = np.exp(-deg_total[leach_dates])  # Leaching days

    # Get total mass by accumulating pesticide_mass_soil and degrading by degradation rate
    total_mass = cumulative_multiply_and_add(pesticide_mass_soil, degradation_rate)

    # Compute conc
    average_conc = ((total_mass / retardation / soil.delta_x) / deg_total) * (1 - degradation_rate)

    # Compute the mass of pesticide in runoff
    runoff_mass[leach_dates] = average_conc[leach_dates] * runoff[leach_dates]  # conc[kg/m3]*[m] = kg/m2

    # Compute the mass of pesticide from erosion
    erosion_dates = np.where((erosion > 0) & (leaching > 0))[0]
    erosion_intensity = soil.erosion_effic / soil.soil_depth  # Assume uniform extraction, no decline, MMF
    enrich = np.exp(2.0 - (0.2 * np.log10(erosion[erosion_dates])))
    enriched_eroded_mass = erosion[erosion_dates] * enrich * kd * erosion_intensity * 0.1
    erosion_mass[erosion_dates] = average_conc[erosion_dates] * enriched_eroded_mass

    return runoff_mass, erosion_mass
