import math
import numpy as np
from collections import defaultdict


def classify_reaches(all_reaches, outlet_dict, waterbody_dict):

    # Get a dictionary of all lentic reaches
    lentic_reaches = {}
    for lake, outlet in outlet_dict.items():
        for reach in waterbody_dict.get(lake) - {outlet}:
            lentic_reaches[reach] = outlet

    # Get all lotic reaches
    lotic_reaches = set(all_reaches) - set(lentic_reaches.keys())

    # Get all reaches in the lake for each lake outlet
    upstream_lentics = defaultdict(set)
    for reach, outlet in lentic_reaches.items():
        upstream_lentics[outlet].add(reach)

    return lotic_reaches, upstream_lentics


def compute_concentration(output_time_series, local_baseflow):

    mass_time_series = output_time_series[0]
    runoff_time_series = output_time_series[1]
    baseflow = local_baseflow.sum(axis=0)
    total_flow = runoff_time_series + baseflow
    concentration = mass_time_series / total_flow
    runoff_conc = mass_time_series / runoff_time_series

    return mass_time_series, runoff_time_series, baseflow, total_flow, concentration, runoff_conc


def convolve_flowing(upstream_reaches, local_output, local_lookup, time_period, interval=1):

    # Sort upstream reaches into 'tanks' based on travel time upstream
    upstream_reaches[1] = np.round(upstream_reaches[1] / interval) * interval  # Round times to the nearest interval
    intervals = sorted(filter(lambda x: x > 0, np.unique(upstream_reaches[1])))  # Every interval containing reaches
    tanks = ((tank, list(upstream_reaches[0, upstream_reaches[1] == tank])) for tank in intervals)  # ((ivl, reach),...)

    # Loop through each tank and perform convolution
    # Trim SAM array to just the reaches in the tank, and sum all reaches together
    output_time_series = np.zeros((2, local_output.shape[2]))
    for interval, tank in tanks:
        row_numbers = list(filter(None, map(local_lookup.get, list(tank))))
        if row_numbers:  # Only proceed if tank reaches were found in the local output (they might all be missing)
            tank_output = local_output[:, row_numbers].sum(axis=1)
            if interval != 0:  # Only perform convolution if timestep is not 0
                irf = impulse_response_function(interval + 1, 1, time_period)
                tank_output[0] = np.convolve(tank_output[0], irf)[:time_period]
                tank_output[1] = np.convolve(tank_output[1], irf)[:time_period]
            output_time_series += tank_output  # Add the convolved tank time series to the total for the reach

    return output_time_series


def convolve_reservoirs(local_output, local_lookup, local_paths, local_lakes, residence_times, time_period):

    # Compile all of the reaches upstream of each lake upstream of the active reach
    upstream_reaches = defaultdict(set)
    all_lakes = np.trim_zeros(np.unique(local_lakes), 'f')  # An array with all the LAKE ids in local_lakes
    for lake in all_lakes:
        path_indices = np.nonzero(np.any(local_lakes == lake, axis=1))[0]  # Row nums of all paths containing the lake
        for path_index in path_indices:
            reaches = local_paths[path_index]
            lake_index = list(local_lakes[path_index]).index(lake)  # The location of the lake in the row (start point)
            upstream_reaches[lake] |= set(reaches[lake_index:])

    # Convolve each reach upstream of a lake against that lakes's residence time
    for lake, reaches in upstream_reaches.items():
        residence_time = residence_times.get(lake, 0.0)
        indices = map(local_lookup.get, reaches)
        if residence_time > 1.5:
            irf = impulse_response_function(1, residence_time, time_period)
            for i in range(2):  # mass and flow
                local_output[i, indices] = np.convolve(local_output[i, indices], irf)[:time_period]

    return local_output


def impulse_response_function(alpha, beta, length):
    def gamma_distribution(t, a, b):
        a, b = map(float, (a, b))
        tau = a * b
        return ((t ** (a - 1)) / (((tau / a) ** a) * math.gamma(a))) * math.exp(-(a / tau) * t)
    return np.array([gamma_distribution(i, alpha, beta) for i in range(length)])


def trim_to_reach(upstream, start_row, col, start_time, end_row):

    # Extract all paths containing the reach in question
    upstream = np.rollaxis(upstream, 2)  # [(paths, times, lakes), maxpaths, maxlength]
    local_upstream = upstream[:, start_row:end_row, col:]
    longest_path = np.argmin((local_upstream[0] > 0), axis=1).max()
    local_paths, local_times, local_lakes = local_upstream[:,:,:longest_path]
    local_times -= start_time

    # Get all unique reach/time combinations
    local_reaches = np.array(list(filter(lambda x: x[0] > 0, set(zip(local_paths.flat, local_times.flat))))).T

    return local_paths, local_times, local_lakes, local_reaches


def trim_to_upstream(sam_output, sam_lookup, all_upstream):

    # Check for missing output
    missing = all_upstream - set(sam_lookup.keys())
    if missing:
        #print("Missing SAM output for {} of {} upstream reaches".format(len(missing), len(all_upstream)))
        pass
    all_upstream = all_upstream & set(sam_lookup.keys())
    indices, reaches = zip(*sorted(((sam_lookup.get(reach), reach) for reach in all_upstream), key=lambda x: x[0]))
    local_lookup = dict(((reach, i) for i, reach in enumerate(reaches)))

    local_output = sam_output[[[0], [1]], indices]
    baseflow = sam_output[2, indices]

    return local_output, baseflow, local_lookup

if __name__ == "__main__":
    print("This is a library. Run travel_time.py")