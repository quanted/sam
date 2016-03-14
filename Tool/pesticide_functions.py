import math
import numpy as np


def applications(input, scenario):
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

        return {1: plant_dates, 2: emergence_dates, 3: mature_dates, 4: harvest_dates}[input.cropstage]

    def days_since(a):
        # Creates a 1d array that counts the number of days since the last positive value
        return np.arange(len(a)) - np.hstack(([0.0], np.maximum.accumulate(np.arange(len(a)) * (a > 0))[:-1]))

    app = input.applications  # Shortcut to application timing details

    # Create local copies of application data arrays
    appmass = np.copy(input.appmass)
    appnumrec = np.copy(input.appnumrec)

    # Initialize variables
    appmass_init = appmass[0]
    appcount = 0

    for date in application_dates():

        if input.distribflag == 1:  # uniform application
            appmass[appcount:appcount + app.twindow1] = appmass_init / app.twindow1

        elif input.distribflag == 2:  # step application
            appmass[appcount:appcount + app.twindow1] = (app.pct1 / 100 * appmass_init) / app.twindow1
            appmass[appcount + app.twindow1:appcount + app.twindow1 + app.twindow2] = \
                (app.pct2 / 100 * appmass_init) / app.twindow2

        elif input.distribflag == 3:  # triangular application
            thalf = app.twindow1 / 2
            leg = ((appmass_init / app.twindow1) / thalf) * ((np.arange(thalf) + 1) ** 2 - (np.arange(thalf) ** 2))
            appmass[appcount:appcount + app.twindow1] = np.hstack((leg, leg[::-1]))

        full_window = app.twindow1 + app.twindow2 if input.distribflag == 2 else app.twindow1
        multiplier = 1 if input.stageflag == 1 else -1
        appnumrec[appcount:appcount + full_window] = np.arange(full_window) + date + (multiplier * input.stagedays) - 1
        appcount += full_window

    daily_mass, daily_method = np.zeros_like(scenario.plant_factor),  np.zeros_like(scenario.plant_factor)
    daily_mass[appnumrec] = appmass
    daily_method[appnumrec] = input.appmethod_init  # MMF - as of now we only offer the user to specify one app method for all dates

    applied_to_soil = daily_mass * (daily_method == 1)
    applied_to_canopy = daily_mass * (daily_method == 2)
    retained_by_canopy = applied_to_canopy * scenario.plant_factor * scenario.covmax
    canopy_to_soil = applied_to_canopy - retained_by_canopy
    canopy_to_soil_days = np.where(canopy_to_soil + scenario.rain > 0)[0]

    pesticide_mass_soil = (applied_to_soil + canopy_to_soil) * soil.distrib_2cm
    degradation = np.exp(-days_since(scenario.rain + canopy_to_soil) * plant.foliar_degradation)
    washoff = np.exp(-scenario.rain * plant.washoff_coeff)

    canopy_mass = 0.0
    for day in canopy_to_soil_days:
        canopy_mass = (canopy_mass + retained_by_canopy[day]) * degradation[day]
        pesticide_mass_soil[day] += canopy_mass - (canopy_mass * washoff[day])
        canopy_mass *= washoff[day]

    return pesticide_mass_soil


def transport(pesticide_mass_soil, scenario, input, process_erosion):
    """
    Simulate transport of pesticide through runoff and erosion
    """
    from Tool.parameters import soil

    # Returns an array in which a value from 'a' is added and 'b' is multiplied on each iteration
    def cumulative_multiply_and_add(a, b):
        return (np.hstack((b[::-1].cumprod()[::-1], 1)) * np.hstack((0, a * b))).cumsum()[1:] / b[::-1].cumprod()[::-1]

    scenario.runoff *= soil.runoff_effic
    leach_dates = np.where(scenario.leaching > 0.0)[0]

    # Erosion efficiency and enrichment
    if process_erosion:
        enrich = math.exp(2.0 - (0.2 * math.log10(scenario.erosion)))
        erosion_intensity = soil.erosion_effic / soil.soil_depth  # Assume uniform extraction, no decline, MMF
        enriched_eroded_mass = (scenario.erosion / 100000.) * enrich  # kg/ha -> g/cm2 (kg/ha*(ha/10000 m2)(1000 g/kg)(m2/10000 cm2) = 1/100000
        scenario.erosion = enriched_eroded_mass * scenario.kd * erosion_intensity * 10000.  # g/cm2 *(m3/g)*10000cm2/m2 -> [m]  !Enrich based on PRZM, kd: kd sediment may differ from kd erosion

    # Get retardation and deg_total
    retardation = (scenario.soil_water / soil.delta_x) + (scenario.bulk_density * scenario.kd)
    deg_total = input.degradation_aqueous + ((scenario.runoff + scenario.leaching) / (soil.delta_x * retardation))

    # Get degradation rate for each day
    degradation_rate = np.full(pesticide_mass_soil.size, np.exp(-input.degradation_aqueous))  # Non-leaching days
    degradation_rate[leach_dates] = np.exp(-deg_total[leach_dates])  # Leaching days

    # Get total mass by accumulating pesticide_mass_soil and degrading by degradation rate
    total_mass = cumulative_multiply_and_add(pesticide_mass_soil, degradation_rate)

    # Compute the mass of pesticide in runoff
    runoff_mass = np.zeros_like(scenario.runoff)
    runoff_mass[leach_dates] = \
        (((total_mass / retardation / soil.delta_x) / deg_total) * (1 - degradation_rate) * scenario.runoff)[
            leach_dates]  # conc[kg/m3]*[m] = kg/m2

    # Compute the mass of pesticide from erosion
    if process_erosion:
        erosion_mass = np.zeros_like(scenario.erosion)
        erosion_mass[leach_dates] = \
            (((total_mass / retardation / soil.delta_x) / deg_total) * (1 - degradation_rate) * scenario.erosion)[
                leach_dates]  # conc[kg/m3]*[m] = kg/m2

    # Offset by one day (JCH - Not sure why I have to offset by a day here, but it works)
    runoff_mass = np.hstack([0.0, runoff_mass[:-1]])
    if process_erosion:
        erosion_mass = np.hstack([0.0, erosion_mass[:-1]])
    else:
        erosion_mass = None

    return runoff_mass, erosion_mass


def waterbody_concentration(q, v, length, total_runoff, runoff_mass, erosion_mass,
                            process_benthic=True, degradation_aqueous=None, koc=None):  # MMF - we need area_wb, daily_depth, depth_0 for benthic calcs to work

    def solute_holding_capacity():
        '''
        Calculates Solute Holding capacities and mass transfer between water column and benthic regions
        '''

        from Tool.parameters import stream_channel

        width = stream_channel.a * np.power(xc, stream_channel.b)
        depth = (xc - width) / 2.0
        surface_area = width * length

        # Aqueous volumes in each region
        vol1 = daily_depth * sa  # total volume in water column, approximately equal to water volume alone
        vol2a = benthic.depth * sa  # total benthic volume
        vol2 = vol2a * benthic.porosity  # with EXAMS Tool.parameters   v2  = VOL2*BULKD*(1.-100./PCTWA)

        # Default EXAMS conditions for partitioning
        kow = koc / .35  # DEFAULT EXAMS CONDITION ON Kow  p.35
        kpdoc1 = kow * .074  # DEFAULT RELATION IN EXAMS (LITTORAL)
        kpdoc2 = koc  # DEFAULT RELATION IN EXAMS (BENTHIC) p.16 of EXAMS 2.98 (or is it Kow*.46 ?)
        xkpb = 0.436 * kow ** .907  # DEFAULT RELATION IN EXAMS

        # mass in littoral region
        vol1a = depth_0 * area_wb  # initial volume corresponding with suspended matter reference
        m_sed_1 = wc.sused * vol1a * .001  # SEDIMENT MASS LITTORAL
        m_bio_1 = wc.plmas * vol1a * .001  # BIOLOGICAL MASS LITTORAL
        m_doc_1 = wc.doc * vol1a * .001  # DOC MASS LITTORAL

        # partitioning coefficients of individual media
        kd_sed_1 = koc * wc.froc * .001  # Kd of sediment in littoral [m3/kg]
        kd_sed_2 = koc * benthic.froc * .001  # Kd of sediment in benthic
        kd_bio = xkpb / 1000.  # Kd of biological organisms
        kd_doc_1 = kpdoc1 / 1000.  # Kd of DOC in littoral region
        kd_doc_2 = kpdoc2 / 1000.  # Kd of DOC in benthic region

        # mass in benthic region
        m_sed_2 = benthic.bulk_density * vol2a * 1000.  # as defined by EXAMS Tool.parameters m_sed_2 = BULKD/PCTWA*VOL2*100000.
        m_bio_2 = benthic.bnmas * sa * .001
        m_doc_2 = benthic.doc * vol2 * .001

        # solute holding capacity in region 1
        capacity_1 = kd_sed_1 * m_sed_1 + kd_bio * m_bio_1 + kd_doc_1 * m_doc_1 + vol1

        # solute holding capacity in region 2
        capacity_2 = kd_sed_2 * m_sed_2 + kd_bio * m_bio_2 + kd_doc_2 * m_doc_2 + vol2

        fw1 = vol1 / capacity_1  # MMF - fw1 will have dimensions of [d], once daily_depth is function of [d] and vol1 is function of [d]
        fw2 = vol2 / capacity_2

        theta = capacity_2 / capacity_1

        sed_conv_factor = vol2 / fw2 / m_sed_2  # converts pore water to Total Conc normalized to sed mass

        # Omega mass transfer
        omega = benthic.d_over_dx / benthic.depth  # (m3/hr)/(3600 s/hr)

        return capacity_1, capacity_2, fw1, fw2, theta, sed_conv_factor, omega

    def simultaneous_diffeq(gamma1, gamma2, omega, theta, m1, m2):
        """
        ANALYTICAL SOLUTION FOR THE TWO SIMULTANEOUS DIFFERNTIAL EQNS:
                  dm1/dt = Am1 + Bm2
                  dm2/dt = Em1 + Fm2
        WITH INITIAL VALUES m1 AND m2 FOR m1 AND m2
        mn1 IS OUTPUT VALUE FOR m1 AFTER TIME T
        mn2 IS OUTPUT VALUE FOR m2 AFTER TIME T
        mavg1 IS AVERAGE VALUE OF m1 OVER TIME T
        """

        t_end = 86400.  # seconds, time step ONE DAY

        # MMF  Reducer - calculate constants for Simuldiff2: A,B,E,F
        a = -gamma1 - omega * theta
        b = omega * theta
        e = omega
        f = -gamma2 - omega

        af = a + f
        fxa = f * a
        bxe = b * e
        dif = 4 * (fxa - bxe)
        bbb = math.sqrt(af * af - dif)

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
        exrt1 = math.exp(rt1)
        exrt2 = math.exp(rt2)
        ccc = x1 * exrt1
        ddd = y1 * exrt2

        mn1 = ccc + ddd
        mn2 = dd * ccc + ee * ddd

        # AVERAGE DAILY CONCENTRATION: SET UP FOR DAILY AVERAGE, BUT CAN BE CHANGED BY CHANGING T1 AND T2
        gx = x1 / root1
        hx = y1 / root2

        term1 = gx * exrt1  # term3 = -X1/root1*exp(root1*T1)
        term2 = hx * exrt2  # term4 = -Y1/root2*exp(root2*T1
        term3 = -gx
        term4 = -hx

        mavg1 = (term1 + term2 + term3 + term4) / t_end  # mavg1=(term1+term2+term3+term4)/(T2-T1)

        mavg2 = (term1 * dd + term2 * ee + term3 * dd + term4 * ee) / t_end  # average compartment 2 concentration

        return mn1, mn2, mavg1, mavg2

    from Tool.parameters import benthic
    from Tool.parameters import soil
    from Tool.parameters import water_column as wc

    n_dates = total_runoff.size

    # Develop hydrograph
    baseflow = np.ones_like(total_runoff) * 1.e-6  # 1.e-6 is minimum baseflow
    average_runoff = np.average(total_runoff)  # m3/d
    baseflow[q >= average_runoff] = q[q >= average_runoff] - average_runoff
    total_flow = baseflow + total_runoff
    total_flow[v == 0] = 0.0  # JCH - We'll need to revisit this
    xc = total_flow / v
    volume = (xc * 40.)

    # Compute daily concentration from runoff
    conc_days = ((runoff_mass > 0.0) & (volume > 0.0))
    daily_concentration = np.zeros_like(total_runoff)
    daily_concentration[conc_days] = runoff_mass[conc_days] / volume[conc_days]

    k_adj = (total_flow / volume) + degradation_aqueous  # MMF - updated with deg aqueous rate, Include within daily loop below because total_flow & volume vary by day?
    exp_k = np.exp(-k_adj)                               # JCH - total flow and volume are already time-series here, so k_adj becomes a time series as well. If neither total flow or volume are
    avgconc_adj = np.zeros_like(total_runoff)            #       modified in the loop then the k_adj calculation doesn't need to be in the loop

    if process_benthic:

        # Compute benthic solute holding capacity
        fw1, fw2, theta, sed_conv_factor, omega = solute_holding_capacity()

        m1_input = runoff_mass + (1. - soil.prben) * erosion_mass
        m2_input = soil.prben * erosion_mass

        # Beginning day aquatic concentrations, considered Peak Aqueous Daily Conc in Water Column
        aq_conc1, aq_conc2 = np.zeros(n_dates), np.zeros(n_dates)
        aqconc_avg1, aqconc_avg2 = np.zeros(n_dates), np.zeros(n_dates)

    else:
        aqconc_avg1 = aqconc_avg2 = aq_conc1 = np.array([])

    # Reset starting values
    conc = mn1 = mn2 = 0

    for d in range(daily_concentration.size):

        # Compute concentration
        conc += daily_concentration[d]
        avgconc_adj[d] = conc / k_adj[d] * (1 - exp_k[d])
        conc *= exp_k[d]

        if process_benthic:

            # MMF  Addition of Benthic region, assume VVWM benthic properties
            m1 = mn1 + m1_input[d]  # daily peak mass
            m2 = mn2 + m2_input[d]

            # JCH - What are these variables for?  They aren't doing anything right now
            m1_store = m1
            m2_store = m2

            # MMF Convert to aqueous concentration
            aq_conc1[d] = m1 * fw1[d] / daily_depth[d] / sa
            aq_conc2[d] = m2 * fw2 / (benthic.depth * sa * benthic.porosity)

            new_aqconc1, new_aqconc2, aqconc_avg1[d], aqconc_avg2[d] = \
                simultaneous_diffeq(k_adj[d], degradation_aqueous, omega, theta, aq_conc1[d], aq_conc2[d])

            # Convert back to masses
            mn1 = new_aqconc1 / fw1[d] * daily_depth[d] * sa
            mn2 = new_aqconc2 / fw2 * benthic.depth * sa * benthic.porosity

            # JCH - This is also not currently being used.
            mavg1_store = aqconc_avg1[d] / fw1[d] * daily_depth[d] * sa

    # Adjust for change in units
    runoff_conc = np.nan_to_num((runoff_mass * 1000000.) / total_runoff)
    avgconc_adj *= 1000000.

    return total_flow, baseflow, avgconc_adj, runoff_conc, aqconc_avg1, aqconc_avg2, aq_conc1
