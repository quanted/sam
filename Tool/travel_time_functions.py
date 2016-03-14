import math
import numpy as np
from collections import defaultdict


def bin_reaches(upstream_paths, upstream_times, interval=1, round_down=True):
    # Divide lotic reaches into tanks based on travel times
    round_func = np.int32 if round_down else np.round
    reach_times = np.array(list(filter(lambda x: x[0] > 0, set(zip(upstream_paths.flat, upstream_times.flat))))).T
    reach_times[1] = round_func(reach_times[1] / interval) * interval  # Round times to the nearest interval
    intervals = sorted(np.unique(reach_times[1]))  # Every interval containing reaches
    tanks = ((tank, list(reach_times[0, reach_times[1] == tank])) for tank in intervals)  # ((ivl, reach),...)

    return tanks

def classify_reaches(all_reaches, outlet_dict, waterbody_dict):
    """
    Classifies reaches as lentic (in a waterbody) or lotic (flowing) and provides a lookup dictionary containing all
    of the lentic reaches upstream of every reach
    :param all_reaches: Set containing the reach IDs of all reaches
    :param outlet_dict: Dictionary containing the outlet reach ID for every water body
    :param waterbody_dict: Dictionary containing the reaches contained within each waterbody
    :return: Set of lentic reach IDs, set of lotic reach IDs, dictionary with reach IDs as keys and a set of all
             upstream lentic reaches as values
    """
    all_reaches = set(all_reaches)

    # Get a dictionary of all lentic reaches
    lentic_reach_dict = {}
    for lake, outlet in outlet_dict.items():
        for reach in waterbody_dict.get(lake) - {outlet}:
            lentic_reach_dict[reach] = outlet

    # Identify lentic and lotic reaches
    lentic_reaches = set(lentic_reach_dict.keys())
    lotic_reaches = all_reaches - lentic_reaches

    # Get all reaches in the lake for each lake outlet
    upstream_lentics = defaultdict(set)
    for reach, outlet in lentic_reach_dict.items():
        upstream_lentics[outlet].add(reach)

    return lentic_reaches, lotic_reaches, upstream_lentics


def compute_concentration(mass_time_series, runoff_time_series, baseflow):
    """
    Divides mass by flow to compute runoff concentrations
    :param mass_time_series: 1d array containing a time series of pesticide mass
    :param runoff_time_series: 1d array containing a time series of runoff volume
    :param baseflow: 1d array containing a time series of baseflow
    :return: Total flow (sum of baseflow and runoff), pesticide concentration, and pesticide concentration in runoff
    """

    total_flow = runoff_time_series + baseflow
    concentration = np.nan_to_num(mass_time_series / total_flow)
    runoff_conc = np.nan_to_num(mass_time_series / runoff_time_series)

    return total_flow, concentration, runoff_conc


def process_flowing(upstream_paths, upstream_times, upstream_output, upstream_lookup, interval=1, convolve=True, diagnose=False):
    """
    Apply convolution to flowing reaches
    :param upstream_paths: 2d array containing all of the flow paths upstream of the active reach
    :param upstream_times: 2d array containing the reach times that correspond to upstream_paths
    :param upstream_output: 2d array containing time series of SAM-generated pesticide loads for each reach
    :param upstream_lookup: A dictionary which provides the row number for each reach in upstream_output
    :param interval:
    :return: A 2d array containing a time series of convolved mass and convolved runoff for each reach
    """

    # Loop through each tank and perform convolution
    n_dates = upstream_output.shape[2]
    output_time_series = np.zeros((2, n_dates))
    tanks = bin_reaches(upstream_paths, upstream_times, interval) # Divide lotic reaches into tanks based on times

    for interval, tank in tanks:  #  Interval: upstream travel time, tank: set of reaches
        if diagnose:
            print(0, interval, list(tank))
        tank_indices = list(filter(lambda x: isinstance(x, int), map(upstream_lookup.get, tank)))
        if tank_indices:  # Only proceed if tank reaches were found in the local output (they might all be missing)
            tank_output = upstream_output[:, tank_indices].sum(axis=1)  # Sum up all time series in the tank
            if interval > 0:
                if convolve:  # Only perform convolution if timestep is not 0
                    irf = impulse_response_function(interval + 1, 1, n_dates)  # Get the convolution function
                    tank_output[0] = np.convolve(tank_output[0], irf)[:n_dates]  # Convolve mass
                    tank_output[1] = np.convolve(tank_output[1], irf)[:n_dates]  # Convolve runoff
                else:  # If convolution is off, just offset the time series for the tank.
                    offset = np.array([0.0] * interval)
                    tank_output[0] = np.hstack((offset, tank_output[0]))[:n_dates]
                    tank_output[1] = np.hstack((offset, tank_output[1]))[:n_dates]

            output_time_series += tank_output  # Add the convolved tank time series to the total for the reach

    return output_time_series


def process_reservoirs(local_output, local_lookup, reach_to_waterbody, local_paths, residence_times):
    """
    Convolves the mass and runoff time series for reaches upstream of reservoirs
    :param local_output: SAM output array trimmed to reaches upstream of the active reach
    :param local_lookup: Lookup dictionary for matching a reach ID to a row number in the SAM output array
    :param local_paths:
    :param local_lakes:
    :param residence_times:
    :return:
    """

    n_dates = local_output.shape[2]

    # Get every lake upstream of the active reach, and all the reaches upstream of each of those lakes
    reaches_above_lake = defaultdict(set)
    upstream_lentics = set(local_lookup.keys()) & set(reach_to_waterbody.keys())  # All lentic reaches upstream of active reach
    for upstream_lentic in upstream_lentics:
        lake = reach_to_waterbody[upstream_lentic]  # The lake corresponding to the lentic reach
        path_indices = np.nonzero(np.any(local_paths == upstream_lentic, axis=1))[0]  # Row nums of all paths containing the reach
        for path_index in path_indices:
            full_path = local_paths[path_index]
            upstream_of_reach = full_path[np.argmax(full_path == upstream_lentic):]
            reaches_above_lake[lake] |= set(upstream_of_reach)

    # Convolve each reach upstream of a lake against that lakes's residence time
    for lake, reaches in reaches_above_lake.items():
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


def preconvolution_report(reach, dates, output_file, upstream_output, baseflow):

    from Tool import write

    # Create output directory
    output_file = output_file.replace("Convolved", "Aggregated")
    runoff_mass, total_runoff = upstream_output.sum(axis=1)
    total_flow, total_conc, runoff_conc = compute_concentration(runoff_mass, total_runoff, baseflow)
    write.daily(output_file.format("aggr", reach), total_conc, runoff_conc, runoff_mass, dates, total_flow, baseflow, total_runoff)


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

def report_progress(i, n_reaches, ivl=50):
    if not ((i + 1) % ivl):
        print("{}/{}".format(i + 1, n_reaches))

if __name__ == "__main__":
    print("This is a library. Run travel_time.py")