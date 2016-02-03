import math
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
    degradation_rate[leaching_dates] = np.exp(-deg_total[leaching_dates])  # Leaching days

    # Get total mass by accumulating pesticide_mass_soil and degrading by degradation rate
    total_mass = cumulative_multiply_and_add(pesticide_mass_soil, degradation_rate)

    runoff_mass = np.zeros_like(runoff)
    runoff_mass[leaching_dates] = \
        (((total_mass / retardation / delta_x) / deg_total) * (1 - degradation_rate) * runoff)[leaching_dates]

    # @@@ - not sure exactly the source of this problem, seems to be related to the -2 indexing in read.scenario
    return np.hstack([0.0, runoff_mass[:-1]])


def waterbody_concentration(q, xc, total_runoff, total_runoff_mass, degradation_aqueous, omega, theta):

    # Develop hydrograph
    baseflow = np.ones_like(total_runoff) * 1.e-6  # 1.e-6 is minimum baseflow
    average_runoff = np.average(total_runoff)  # m3/d
    volume = (xc * 40.)
    baseflow[q >= average_runoff] = q[q >= average_runoff] - average_runoff
    total_flow = baseflow + total_runoff

    # Compute daily concentration from runoff
    conc_days = ((total_runoff_mass > 0.0) & (volume > 0.0))
    daily_concentration = np.zeros_like(total_runoff)
    daily_concentration[conc_days] = total_runoff_mass[conc_days] / volume[conc_days]

    # Modify concentration @@@ - this doesn't seem to do anything.  exp kadj always 0 or 1
    # MMF - We should make sure flow, volume, and k_adj calcs are within daily loop
    k_adj = total_flow / volume  # + degradation_aqueous  #MMF k_adj[d] = total_flow[d]/volume[d]

    # MMF - Updates to water body - to include water column and static benthic regions
    # Addition of gamma1 and gamma2 degradation rates
    gamma1 = k_adj[d]
    gamma2 = degradation_aqueous

    # to avoid numerical difficulties
    gamma1 = min((gamma1, 1e-18))
    gamma2 = min((gamma2, 1e-18))

    exp_k = np.exp(-k_adj)
    avgconc_adj = np.zeros_like(total_runoff)
    conc = 0

    # MMF - reset starting values
    m1 = 0
    m2 = 0
    mn1 = 0
    mn2 = 0
    d elt = 86400.  # seconds, time step ONE DAY

    # MMF  Reducer - calculate constants for Simuldiff2: A,B,E,F
    a = -gamma1 - omega * theta
    b = omega * theta
    e = omega
    f = -gamma2 - omega

    for d in range(daily_concentration.size):
        conc += daily_concentration[d]
        avgconc_adj[d] = conc / k_adj[d] * (1 - exp_k[d])
        conc *= exp_k[d]

        # MMF  Addition of Benthic region, assume VVWM benthic properties
        # *** Not sure if we need "[d]" to be included here for variables or not?
        m1 = mn1 + m1_input[d]  # daily peak mass
        m2 = mn2 + m2_input[d]
        m1_store[d] = m1
        m2_store[d] = m2

        # MMF Convert to aqueous concentration
        aqconc1 = m1 * fw1[d] / daily_depth[d] / area_wb
        aqconc2 = m2 * fw2 / (benthic_depth * area_wb * porosity )

        # ******************************************************
        # store these beginning day aquatic concentrations
        # these variables are only used for delivery to output routines
        aq1_store[d] = aqconc1
        aq2_store[d] = aqconc2
        # ******************************************************

        # MMF - Use Analytical Solution for 2 Simultaneous Differential Equations from VVWM
        # dm1/dt = Am1 + Bm2
        # dm2/dt = Em1 + Fm2
        # WITH INITIAL VALUES m1 AND m2 FOR m1 AND m2
        # mn1 IS OUTPUT VALUE FOR m1 AFTER TIME T
        # mn2 IS OUTPUT VALUE FOR m2 AFTER TIME T
        # mavg1 IS AVERAGE VALUE OF m1 OVER TIME T

        new_aqconc1, new_aqconc2, aqconc_avg1, aqconc_avg2 = simultaneous_diffeq(a, b, e, f, aqconc1, aqconc2, delt)

        # Convert back to masses
        mn1 = new_aqconc1 / fw1[d] * daily_depth[d] * area_wb
        mn2 = new_aqconc2 / fw2 * benthic_depth * area_wb * porosity

        mavg1_store[d] = aqconc_avg1[d] / fw1[d] * daily_depth[d] * area_wb

        # Simple concentration calc - for comparsion
        conc_simple[d] = Total_Runoff_Mass[d] / q_adj_tot[d]

    # Adjust for change in units
    runoff_conc = np.nan_to_num((total_runoff_mass * 1000000.) / total_runoff)
    avgconc_adj *= 1000000.

    return total_flow, baseflow, avgconc_adj, runoff_conc


def pesticide_applications(plant_factor, stageflag, distribflag, appnumrec_global, appmass_global,
                           appmethod_init, app_windows, stagedays, cropstage, rain, soil_2cm, covmax, foliar_deg,
                           washoff_coeff):
    """
    Creates a time-series array of pesticide applications and transfers those to the soil
    """

    def days_since(a):
        """
        Creates a 1d array that counts the number of days since the last positive value
        """
        return np.arange(len(a)) - np.hstack(([0.0], np.maximum.accumulate(np.arange(len(a)) * (a > 0))[:-1]))

    def application_dates():
        """
        Predicts the dates of pesticide application based on crop stage
        """
        cumulative_maturity = np.int16((plant_factor == 1).cumsum())
        overcount = np.maximum.accumulate(cumulative_maturity * (plant_factor == 0))
        count_mature = np.concatenate(([0], (cumulative_maturity - overcount)[:-1]))
        n_plant = np.concatenate(([0], ((plant_factor == 0) & (count_mature > 0)).cumsum()))

        harvest_dates = np.where(n_plant[:-1] != n_plant[1:])[0]
        plant_dates = harvest_dates - (2 * count_mature[harvest_dates])
        emerg_dates = plant_dates + 7
        mature_dates = harvest_dates - count_mature[harvest_dates]

        return {1: plant_dates, 2: emerg_dates, 3: mature_dates, 4: harvest_dates}[cropstage]


    # Create local copies of application data arrays
    appmass = np.copy(appmass_global)
    appnumrec = np.copy(appnumrec_global)

    twindow1, pct1, twindow2, pct2 = app_windows
    appmass_init = appmass[0]
    appcount = 0

    for date in application_dates():

        # @@@ - if (cropstage, distribflag, stageflag) == (1, 2, 1): there are gaps (e.g. at 840)
        # @@@ with start and end date offsets (tp1 = 0 v tp1 = 1, etc)
        # MMF - gaps likely stem from the version of the Fortran code, which I don't think was completely updated for diff app types,
        # let's try to remove gaps here in possible combinations

        if (cropstage, distribflag) in [(1, 1), (2, 1), (3, 2), (3, 3), (4, 2), (4, 3)]:  # @@@
            appcount = 0

        if distribflag == 1:  # uniform application
            appmass[appcount:appcount + twindow1] = appmass_init / twindow1

        elif distribflag == 2:  # step application
            appmass[appcount:appcount + twindow1] = (pct1 / 100 * appmass_init) / twindow1
            appmass[appcount + twindow1:appcount + twindow1 + twindow2] = (pct2 / 100 * appmass_init) / twindow2

        elif distribflag == 3:  # triangular application
            thalf = twindow1 / 2
            leg = ((appmass_init / twindow1) / thalf) * ((np.arange(thalf) + 1) ** 2 - (np.arange(thalf) ** 2))
            appmass[appcount:appcount + twindow1] = np.hstack((leg, leg[::-1]))

        full_window = twindow1 + twindow2 if distribflag == 2 else twindow1
        multiplier = 1 if stageflag == 1 else -1
        appnumrec[appcount:appcount + full_window] = np.arange(full_window) + date + (multiplier * stagedays) - 1
        appcount += full_window

    applications = np.zeros((2, plant_factor.size))
    applications[0][appnumrec] = appmass
    applications[1][appnumrec] = appmethod_init  # @@@ - how does any other application method get used?
    # MMF - as of now we only offer the user to specify one app method for all app dates in a given simulation,
    # we could expand to allow users to enter diff app methods for each app date,
    # but then user input will begin to look like PWC interface (list of each app date and corresponding method),
    # may want to keep simple for now, and just have one app method/simulation

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



def solute_holding_capacity(koc):
    '''
    # Calculates Solute Holding capacities and mass transfer between water column and benthic regions
    '''

    # Input variables needed from NHD or other source
    area_wb = 50000  # placeholder, surface area of water body (m2) - will need to derive and/or use NHD if possible
    daily_depth = 2  # placeholder, daily depth of water body (m)
    depth_0 = 2  # placeholder, initial water body depth (m)

    # Water Column Parameters - USEPA OPP defaults
    dfac = 1.19  # photolysis parameter from VVWM
    sused = 30  # water column susp solid conc (mg/L)
    chloro = 0.005  # water column chlorophyll conc (mg/L)
    froc1 = 0.04  # water column organic carbon fraction on susp solids
    doc1 = 5  # water column dissolved organic carbon content (mg/L)
    plmas = 0.4  # water column biomass conc (mg/L)

    # Benthic Parameters - USEPA OPP defaults
    benthic_depth = 0.05  # benthic depth (m)
    porosity = 0.5  # benthic porosity
    bulk_density = 1.35  # bulk density, dry solid mass/total vol (g/cm3)
    froc2 = 0.04  # benthic organic carbon fraction
    doc2 = 5  # benthic dissolved organic carbon content (mg/L)
    bnmas = 0.006  # benthic biomass intensity (g/m2)
    d_over_dx = 1e-8  # mass transfer coefficient for exchange between benthic and water column (m/s)

    # Aqueous volumes in each region
    vol1a = daily_depth * area_wb  # total volume in water column, approximately equal to water volume alone
    vol2a = benthic_depth * area_wb  # total benthic volume
    vol2 = vol2a * porosity  # with EXAMS parameters   v2  = VOL2*BULKD*(1.-100./PCTWA)

    # Default EXAMS conditions for partitioning
    kow = koc / .35  # DEFAULT EXAMS CONDITION ON Kow  p.35
    kpdoc1 = kow * .074  # DEFAULT RELATION IN EXAMS (LITTORAL)
    kpdoc2 = koc  # DEFAULT RELATION IN EXAMS (BENTHIC) p.16 of EXAMS 2.98 (or is it Kow*.46 ?)
    xkpb = 0.436 * kow ** .907  # DEFAULT RELATION IN EXAMS

    # mass in littoral region
    vol1a = depth_0 * area_wb  # initial volume corresponding with susspended matter reference
    m_sed_1 = sused * vol1a * .001  # SEDIMENT MASS LITTORAL
    m_bio_1 = plmas * vol1a * .001  # BIOLOGICAL MASS LITTORAL
    m_doc_1 = doc1 * vol1a * .001  # DOC MASS LITTORAL

    # partitioning coefficients of individual media
    kd_sed_1 = koc * froc1 * .001  # Kd of sediment in littoral [m3/kg]
    kd_sed_2 = koc * froc2 * .001  # Kd of sediment in benthic
    kd_bio = xkpb / 1000.  # Kd of biological organisms
    kd_doc_1 = kpdoc1 / 1000.  # Kd of DOC in littoral region
    kd_doc_2 = kpdoc2 / 1000.  # Kd of DOC in benthic region

    # mass in benthic region
    m_sed_2 = bulk_density * vol2a * 1000.  # as defined by EXAMS parameters m_sed_2 = BULKD/PCTWA*VOL2*100000.
    m_bio_2 = bnmas * area_wb * .001
    m_doc_2 = doc2 * vol2 * .001

    # solute holding capacity in region 1
    capacity_1 = kd_sed_1 * m_sed_1 + kd_bio * m_bio_1 + kd_doc_1 * m_doc_1 + vol1

    # solute holding capacity in region 2
    capacity_2 = kd_sed_2 * m_sed_2 + kd_bio * m_bio_2 + kd_doc_2 * m_doc_2 + vol2

    fw1 = vol1 / capacity_1
    fw2 = vol2 / capacity_2

    theta = capacity_2 / capacity_1

    sed_conv_factor = vol2 / fw2 / m_sed_2  # converts pore water to Total Conc normalized to sed mass

    # Omega mass transfer
    omega = d_over_dx / benthic_depth  # (m3/hr)/(3600 s/hr)

    return capacity_1, capacity_2, fw1, fw2, theta, sed_conv_factor, omega


def simultaneous_diffeq(a, b, e, f, m1, m2, t_end):
    """
    # ANALYTICAL SOLUTION FOR THE TWO SIMULTANEOUS DIFFERNTIAL EQNS:
    #          dm1/dt = Am1 + Bm2
    #          dm2/dt = Em1 + Fm2
    #  WITH INITIAL VALUES m1 AND m2 FOR m1 AND m2
    #  mn1 IS OUTPUT VALUE FOR m1 AFTER TIME T
    #  mn2 IS OUTPUT VALUE FOR m2 AFTER TIME T
    #  mavg1 IS AVERAGE VALUE OF m1 OVER TIME T
    """

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
