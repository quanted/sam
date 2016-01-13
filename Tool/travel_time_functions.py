import math
import numpy as np
from collections import defaultdict


def classify_reaches(all_reaches, outlet_dict, waterbody_dict):

    all_reaches = set(all_reaches)

    # Get a dictionary of all lentic reaches
    lentic_reach_dict = {}
    for lake, outlet in outlet_dict.items():
        for reach in waterbody_dict.get(lake) - {outlet}:
            lentic_reach_dict[reach] = outlet

    lentic_reaches = set(lentic_reach_dict.keys())
    lotic_reaches = all_reaches - lentic_reaches

    # Get all reaches in the lake for each lake outlet
    upstream_lentics = defaultdict(set)
    for reach, outlet in lentic_reach_dict.items():
        upstream_lentics[outlet].add(reach)

    return lentic_reaches, lotic_reaches, upstream_lentics


def compute_concentration(mass_time_series, runoff_time_series, baseflow):

    total_flow = runoff_time_series + baseflow
    concentration = np.nan_to_num(mass_time_series / total_flow)
    runoff_conc = np.nan_to_num(mass_time_series / runoff_time_series)

    return mass_time_series, runoff_time_series, baseflow, total_flow, concentration, runoff_conc


def convolve_flowing(upstream_paths, upstream_times, upstream_output, upstream_lookup, interval=1):

    # Divide lotic reaches into tanks based on travel times
    reach_times = np.array(list(filter(lambda x: x[0] > 0, set(zip(upstream_paths.flat, upstream_times.flat))))).T
    reach_times[1] = np.round(reach_times[1] / interval) * interval  # Round times to the nearest interval
    intervals = sorted(np.unique(reach_times[1]))  # Every interval containing reaches
    tanks = ((tank, list(reach_times[0, reach_times[1] == tank])) for tank in intervals)  # ((ivl, reach),...)

     # Loop through each tank and perform convolution
    n_dates = upstream_output.shape[2]
    output_time_series = np.zeros((2, n_dates))
    for interval, tank in tanks:
        tank_indices = list(filter(lambda x: isinstance(x, int), map(upstream_lookup.get, tank)))
        if tank_indices:  # Only proceed if tank reaches were found in the local output (they might all be missing)
            tank_output = upstream_output[:, tank_indices].sum(axis=1)
            if interval > 0:  # Only perform convolution if timestep is not 0
                irf = impulse_response_function(interval, 1, n_dates)
                tank_output[0] = np.convolve(tank_output[0], irf)[:n_dates]
                tank_output[1] = np.convolve(tank_output[1], irf)[:n_dates]
            output_time_series += tank_output  # Add the convolved tank time series to the total for the reach

    return output_time_series


def convolve_reservoirs(local_output, local_lookup, local_paths, local_lakes, residence_times, n_dates):

    n_dates = local_output.shape[2]

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
            irf = impulse_response_function(1, residence_time, n_dates)
            for i in range(2):  # mass and flow
                local_output[i, indices] = np.convolve(local_output[i, indices], irf)[:n_dates]

    return local_output


def impulse_response_function(alpha, beta, length):
    def gamma_distribution(t, a, b):
        a, b = map(float, (a, b))
        tau = a * b
        return ((t ** (a - 1)) / (((tau / a) ** a) * math.gamma(a))) * math.exp(-(a / tau) * t)
    return np.array([gamma_distribution(i, alpha, beta) for i in range(length)])


def preconvolution_report(reach, dates, output_dir, output_format, upstream_output, baseflow):

    import os
    import write

    # Create output directory
    output_dir.replace("Convolved", "Aggregated")
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    runoff_mass, total_runoff = upstream_output.sum(axis=1)

    runoff_mass, total_runoff, baseflow, total_flow, total_conc, runoff_conc = \
        compute_concentration(runoff_mass, total_runoff, baseflow)

    write.daily(output_dir, output_format, str(reach) + "_ag", total_conc, runoff_conc, runoff_mass, dates,
                 total_flow, baseflow, total_runoff)


def trim_to_reach(upstream, start_row, col, start_time, end_row):

    # Extract all paths containing the reach in question
    upstream = np.rollaxis(upstream, 2)  # [(paths, times, lakes), maxpaths, maxlength]
    local_upstream = upstream[:, start_row:end_row, col:]
    longest_path = np.argmin((local_upstream[0] > 0), axis=1).max()
    local_paths, local_times, local_lakes = local_upstream[:,:,:longest_path]
    adjusted_times = local_times - start_time

    return local_paths, local_lakes, adjusted_times


def trim_to_upstream(sam_output, sam_lookup, upstream_paths):

    all_upstream = set(np.unique(upstream_paths))

    # Check for missing output
    missing = all_upstream - set(sam_lookup.keys())
    if missing:
        #print("Missing SAM output for {} of {} upstream reaches".format(len(missing), len(all_upstream)))
        pass
    all_upstream = all_upstream & set(sam_lookup.keys())
    indices, reaches = zip(*sorted(((sam_lookup.get(reach), reach) for reach in all_upstream), key=lambda x: x[0]))
    local_lookup = dict(((reach, i) for i, reach in enumerate(reaches)))
    local_output = np.copy(sam_output[[[0], [1]], indices])

    return local_output, local_lookup


if __name__ == "__main__":
    print("This is a library. Run travel_time.py")