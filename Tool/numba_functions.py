import numpy as np
from numba import jit, njit, int32, double
from numba.types import NPDatetime


@njit
def initialize_soil(delta_x, increments_1, increments_2, bd_5, fc_5, wp_5, bd_20, fc_20, wp_20):
    """ Initialize soil properties """
    soil_properties = np.zeros((5, increments_1 + increments_2))

    for i in range(increments_1):
        soil_properties[0, i] = bd_5 * 1000.  # kg/m3
        soil_properties[1, i], soil_properties[2, i], soil_properties[3, i] = \
            fc_5 * delta_x[i], fc_5 * delta_x[i], wp_5 * delta_x[i]
    for i in range(increments_1, increments_1 + increments_2):
        soil_properties[0, i] = bd_20 * 1000.  # kg/m3
        soil_properties[1, i], soil_properties[2, i], soil_properties[3, i] = \
            fc_20 * delta_x[i], fc_20 * delta_x[i], wp_20 * delta_x[i]

    # Generate a cumulative vector of depths
    cumsum = 0  # JCH - Not using cumsum here because I don't know how to put it in an array
    for i in range(delta_x.size):
        cumsum += delta_x[i]
        soil_properties[4, i] = cumsum

    return soil_properties


@njit
def rain_and_snow(precip, temp, sfac):
    """ Simplified for use with numba"""
    rain_and_melt = np.zeros((2, precip.size))
    snow_accumulation = 0.
    for i in range(precip.size):
        if temp[i] <= 0:
            snow_accumulation += precip[i]
        else:
            rain_and_melt[0, i] = precip[i]
            snow_melt = min(snow_accumulation, sfac * temp[i])
            snow_accumulation -= snow_melt
            rain_and_melt[1, i] = precip[i] + snow_melt
    return rain_and_melt


@njit
def interpolate_plant_stage(n_days, emergence, maturity, harvest, fallow_val, crop_val):
    plant_factor = np.ones(n_days) * fallow_val
    for i in range(emergence.size):
        period = maturity[i] - emergence[i]
        for j in range(period + 1):
            plant_factor[emergence[i] + j] = fallow_val + (j * (crop_val / period))
        plant_factor[maturity[i]:harvest[i] + 1] = crop_val
    return plant_factor


@njit
def find_node(n, depth, target_depth):
    n -= 1  # Zero indexing
    if target_depth >= depth[n]:
        node = n - 1
    elif target_depth <= depth[1]:
        node = 0
    else:
        for node in range(depth.size):
            if target_depth <= depth[node]:
                break
        if depth[node] - target_depth > target_depth - depth[node - 1]:  # select the node that's closer to the depth
            node -= 1
    return node


@njit
def surface_hydrology(field_m, wilt_m, plant_factor, cn, depth, soil_water_m,  # From other function output
                      irr_type, irr_depletion, anetd, root_max, leaching_factor, cintcp,  # From scenario
                      num_records, effective_rain, rain, potential_et,  # From metfile
                      number_soil_incr, delta_x):  # From parameters


    # Initialize arrays
    maxdays = plant_factor.size
    velocity_all = np.zeros((number_soil_incr, maxdays))
    soil_water_m_all = np.zeros((number_soil_incr, maxdays))
    runoff = np.zeros(maxdays)  # by initialization, runoff is zero when rain < 0.2S
    velocity = np.zeros(number_soil_incr)
    et_factor = np.zeros(number_soil_incr)
    available_water_m = np.zeros(number_soil_incr)
    et_depth = np.zeros(num_records)
    et_node = np.zeros(num_records, dtype=np.int32)

    fc_minus_wp = field_m - wilt_m
    if irr_type > 0:
        irrigation_node = find_node(number_soil_incr, depth, root_max)
        target_dryness = 0
        for i in range(irrigation_node):
            target_dryness += fc_minus_wp[i] * irr_depletion + wilt_m[i]
        total_fc = np.sum(field_m[:irrigation_node])

    evapo_node = find_node(number_soil_incr, depth, anetd)  # node only for evaporation
    et_node[:] = evapo_node  # initially set all to the minimum

    for i in range(plant_factor.size - 1):
        et_depth[i] = anetd
        if (plant_factor[i] > 0) and (plant_factor[i] * root_max) > anetd:
            et_depth[i] = plant_factor[i] * root_max
            if et_depth[i] > anetd:
                et_node[i] = find_node(number_soil_incr, depth, et_depth[i])

    canopy_holdup = 0

    for day in range(num_records):

        soil_layer_loss = np.zeros(number_soil_incr)

        # Process irrigation and modify
        overcanopy_irrigation = 0
        if irr_type > 0:
            current_dryness = np.sum(soil_water_m[:irrigation_node])
            daily_max_irrigation = 0.2 * ((2540. / cn[day]) - 25.4) / 100.
            if current_dryness < target_dryness and effective_rain[day] <= 0.:
                # water to be added to bring irrigation zone to field capacity
                irrig_required = (total_fc - current_dryness) * leaching_factor + 1.
            if irr_type == 3:
                overcanopy_irrigation = min(irrig_required, daily_max_irrigation)
                effective_rain[day] = overcanopy_irrigation
            elif irr_type == 4:  # undercanopy irrigatoin
                effective_rain[day] = min(irrig_required, daily_max_irrigation)

        # Determine daily runoff (Soil moisture curve number option)
        if effective_rain[day] > 0:
            s = 25.4 / cn[day] - .254  # Revised, CN soil moisture adjustment removed, mmf 9/2015

            if effective_rain[day] > (0.2 * s):  # Runoff by the Curve Number Method
                runoff[day] = max(0, (effective_rain[day] - (0.2 * s)) ** 2 / (effective_rain[day] + (0.8 * s)))

        # Leaching into top layer
        leaching = effective_rain[day] - runoff[day]
        if rain[day] > 0. or overcanopy_irrigation > 0:
            available_canopy_gain = (rain[day] + overcanopy_irrigation) * (1. - runoff[day] / effective_rain[day])
            delta_water = min(available_canopy_gain, cintcp * plant_factor[day] - canopy_holdup)
            canopy_holdup += delta_water
            leaching -= delta_water

        # update canopy holdup here
        et_from_canopy = canopy_holdup - potential_et[day]
        canopy_holdup = max(0., et_from_canopy)
        available_soil_et = max(0., -et_from_canopy)

        # Reduction factor below 0.6 available water
        check_moisture_et, target_moisture_et = 0., 0.
        for i in range(et_node[day] + 1):
            available_water_m[i] = soil_water_m[i] - wilt_m[i]
            check_moisture_et += soil_water_m[i] - wilt_m[i]
            target_moisture_et += 0.6 * fc_minus_wp[i]

        if check_moisture_et < target_moisture_et:
            available_soil_et *= check_moisture_et / target_moisture_et

        # Reductions by depth and available water
        for i in range(et_node[day] + 1):
            et_factor[i] = (depth[et_node[day]] - depth[i] + delta_x[i]) * available_water_m[i]

        # Normalize ET factor and set to zero if it's dry
        et_sum = np.sum(et_factor[:et_node[day] + 1])
        if et_sum > 0:
            for i in range(et_node[day] + 1):
                et_factor[i] /= et_sum
        else:
            et_factor[:] = 0.

        # Calculate soil layer loss (# JCH - somewhere in here is an array that doesn't need to be an array. increments zero?
        for i in range(et_node[day] + 1):
            soil_layer_loss[i] = available_soil_et * et_factor[i]  # potential loss

        # Leaching loop
        last_velocity = leaching
        for node in range(number_soil_incr):
            water_level = last_velocity - soil_layer_loss[node] + soil_water_m[node]
            if water_level > field_m[node]:
                velocity[node] = water_level - field_m[node]
                soil_water_m[node] = field_m[node]
            else:
                velocity[node] = 0.
                soil_water_m[node] = max(water_level, wilt_m[node])
            if velocity[node] <= 0. and node > et_node[day]:
                velocity[node:number_soil_incr] = 0.
                break
            last_velocity = velocity[node]

        # Here's the output for this loop
        for i in range(velocity.size):
            velocity_all[i, day] = velocity[i]
            soil_water_m_all[i, day] = soil_water_m[i]

    return runoff, soil_water_m_all, velocity_all


@njit
def process_erosion(num_records, slope, manning_n, runoff, rain, cn, usle_klscp, raintype):

    # Initialize output
    erosion_loss = np.zeros(num_records)

    # Time of concentration, TR-55
    l_sheet = min(100. * np.sqrt(slope) / manning_n, 100.)
    l_shallow = (100. * np.sqrt(slope) / manning_n) - l_sheet

    for i in range(runoff.size):
        if runoff[i] > 0.:
            t_conc = 0
            ia_over_p = 0
            if slope > 0:
                t_conc = 0.007 * (manning_n * l_sheet) ** 0.8 / \
                         np.sqrt(rain[i] / 0.0254) / \
                         slope ** 0.4 + l_shallow / 58084.2 / np.sqrt(slope)
            if rain[i] > 0:
                ia_over_p = .0254 * (200. / cn[i] - 2.) / rain[i]  # 0.2 * s, in inches

            # lower and upper limit according to TR-55
            if ia_over_p <= 0.1:
                c = raintype[0]
            elif ia_over_p >= 0.5:
                c = raintype[8]
            else:  # interpolation of intermediate. clunky because numba
                lower = (20. * (ia_over_p - 0.05)) - 1
                delta = raintype[int(lower) + 1] - raintype[int(lower)]
                interp = (lower % 1) * delta
                c = raintype[int(lower)] + interp

            # Erosion loss (kg) - MUSLE equation (Williams, 1975)
            peak_discharge = 10. ** (c[0] + c[1] * np.log10(t_conc) + c[2] * (np.log10(t_conc)) ** 2)
            qp = 1.54958679 * runoff[i] * peak_discharge
            erosion_loss[i] = 1.586 * (runoff[i] * 1000. * qp) ** .56 * usle_klscp[i] * 1000.  # kg/d

    return erosion_loss

