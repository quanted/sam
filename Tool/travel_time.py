import logging
import time
from functools import partial

from Tool import read
from Tool.parameters import paths
from Tool import travel_time_functions as functions

# To do:
# * benthic transport
# * migrate read functions to read.py?


def timeit(method):

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        print('Time start = %s, Time end = %s: %2.2f sec' %
              (ts, te, te-ts))
        return result

    return timed


def callback_mp(future):
    logging.info(future.exception())


def reach_calc(reaches):
    """
    Function to send to each process when using multiprocessing. The number total number of reaches processed is
    first sent through reach_ranges_mp() to generate ranges to divide up the work that is sent into this function.
    This function is a batch of reaches to be processed in a single process.
    :param reaches: list of Reach class instances
    """
    for reach in reaches:
        if not reach.run:
            try:
                reach.process()
            except Exception as e:
                print(str(e))


def time_of_travel(input_data):
    """
    Jon F. - Changed region.cascade() method to yield on Lake Bins, and moved Lake looping into this method
    This was done to all the multiprocess to easily separate each Lake Bin into separate sections to keep the
    sequence of calculations correct (so the Reservoir convolution would never occur before its upstream reaches
    had been convoluted, as the multiprocess processing could not guarantee that without breaking on Lake Bins.

    :param input_data:
    :return:
    """

    # Object containing input data parameters
    inputs = read.InputParams(input_data, mode="TimeOfTravel")

    # Output on/off. Included primarily for diagnostic purposes to limit output
    write_to_file = input_data['write']

    # Structure which generates and stores impulse response functions
    irf = functions.ImpulseResponse()

    # Initializes a time of travel run for an NHD region
    region = functions.Region(inputs, paths, irf, write_to_file)

    if input_data["multiprocess"]:
        print("multiprocess")
        from Tool import mp
        from concurrent.futures import wait, ALL_COMPLETED

        mp = mp.Multiprocessing(input_data['no_of_processes'])
        pool = mp.setup()

        # lake_bin_counter = 1
        for lake_bin, lakes, no_of_lakes in region.cascade():
            # print("Processing Lake Bin #%s, w/ No of lakes = %s" % (lake_bin_counter, no_of_lakes))

            futures = []  # Container to store each future (each job submitted to executor)

            if lakes is not None:

                for lake in lakes:
                    # Loop over each lake within the Lake Bin, sending each lake to its own process
                    # If the # of upstream reaches in a lake is large, batch them into multiple processes

                    # reaches = filter(None, map(region.reaches.get, lake.upstream_reaches))  # Generator version

                    # Changed "reaches" to a list so it can be indexed when chunking it
                    reaches = [region.reaches.get(x) for x in lake.upstream_reaches]  # list version

                    no_upstream_reaches = lake.upstream_reaches.size

                    if no_upstream_reaches < 64:
                        # When the # of upstream_reaches is small, send all of them to a single process
                        futures.append(
                            pool.submit(reach_calc(reaches))
                        )

                        # TODO: Old approach, where each reach is sent to its own process == no bueno
                        # for reach in reaches:
                        #     # Process reaches upstream of the reservoir (lake)
                        #     # reach.run = True  # TODO: Fake running the reaches
                        #     if not reach.run:
                        #         # futures.append(pool.submit(reach.process()))  # Submit job to executor
                        #         pool.submit(reach.process())

                    else:
                        # When the # of upstream_reaches is large, chunk them and send them to multiple processes
                        ranges = mp.chunk(int(no_upstream_reaches))

                        for rng in ranges:
                            # Submit batches of reaches to a process based on 'chunk_size'
                            futures.append(
                                pool.submit(
                                    reach_calc(
                                        reaches[rng[0]:rng[1]]
                                    )
                                )
                            )

                    # Process the reservoir (lake)
                    # Note: Last batch of reaches will not have an associated reservoir (lake)
                    futures.append(pool.submit(lake.process()))  # Submit job to executor

            else:
                # Loop through reaches that have not yet been run (ostensibly those not upstream of any reservoir (lake)

                print("No Lakes")
                t_lakes = time.time()
                print(t_lakes)  # Print time it takes to get to last lake bin for Ohio River Valley (Region 05)
                # remaining_reaches = filter(lambda x: not x.run, region.reaches.values())  # Generator version
                remaining_reaches = [x for x in region.reaches.values() if not x.run]  # list version

                # TODO: Sequential

                # no_of_remaining_reaches = 1
                # for reach in remaining_reaches:
                #     no_of_remaining_reaches += 1
                #     if not reach.run:  # TODO: This is unnecessary
                #         reach.process()
                # print("Number of reaches remaining: %s" % no_of_remaining_reaches)

                # TODO: Parallel
                print("Number of 'remaining_reaches': %s" % len(remaining_reaches))
                ranges = mp.chunk(len(remaining_reaches))

                for rng in ranges:
                    # Submit batches of reaches to a process based on 'chunk_size'
                    futures.append(
                        pool.submit(
                            reach_calc(
                                remaining_reaches[rng[0]:rng[1]]
                            )
                        )
                    )

            # Wait until all futures (jobs submitted to executor) are finished before continuing to next Lake Bin
            wait(futures, return_when=ALL_COMPLETED)
            # mp.wait_to_finish(futures)

            # lake_bin_counter += 1
        # print("Total number of Lake Bins: %s" % lake_bin_counter)
        pool.shutdown()

    else:  # Sequential
        print("Sequential")
        lake_bin_counter = 1
        for lake_bin, lakes, no_of_lakes in region.cascade():

            # print("No of lakes = %s" % no_of_lakes)
            if lakes is not None:
                lake_counter = 0
                for lake in lakes:  # Loop over each lake with in the Lake Bins
                    lake_counter += 1
                    # print("%s, " % lake_bin_counter, end="")  # Print the Lake Bin #
                    # print("%s, " % lake_counter, end="")  # Print lake #
                    reaches = filter(None, map(region.reaches.get,
                                     lake.upstream_reaches))  # filter(None, map(self.reaches.get, lake.upstream_reaches))
                    reach_counter = 0
                    for reach in reaches:
                        # Process reaches upstream of the reservoir (lake)
                        if not reach.run:
                            reach_counter += 1
                            # reach.run = True
                            reach.process()
                    # print("%s, " % reach_counter)  # Print the number of upstream reaches for the lake

                    # Process the reservoir (lake)
                    # Note: Last batch of reaches will not have an associated reservoir (lake)
                    lake.process()

            else:
                t_lakes = time.time()
                print(t_lakes)
                # print("%s, " % lake_bin_counter, end="")  # Print the Lake Bin #
                # Loop through reaches that have not yet been run (ostensibly those not upstream of any reservoir (lake)
                # print("0, ", end="")  # No lakes

                # print("No Lakes")
                # print(time.time())  # Print time it takes to get to last lake bin for Ohio River Valley (Region 05)
                remaining_reaches = filter(lambda x: not x.run, region.reaches.values())
                for reach in remaining_reaches:
                    if not reach.run:
                        # reach.run = True
                        reach.process()

            lake_bin_counter += 1

    dayshed_counter, convolution_counter = 0, 0
    for reach in region.reaches.values():
        dayshed_counter += reach.daysheds
        convolution_counter += reach.times_convolved

    return t_lakes, dayshed_counter, convolution_counter


# @timeit
def main(input_data=None, write=False, mp=False, nproc=16, log=False):

    if input_data is None:
        input_data = {"inputs":
                          {"mode": "convolved",
                           "region": "05",
                           "sam_output_id": "dummy_region_05",
                           },
                      "run_type": "single",
                      "write": write,
                      "multiprocess": mp,
                      "no_of_processes": nproc  # Set to False if you want to use the no. of CPU cores on machine
                      }
    print("Write = %s, Multiprocessing = %s, Logging = %s" % (write, mp, log))
    t_start = time.time()
    print(t_start)
    output, no_daysheds, times_convolved = time_of_travel(input_data)
    t_end = time.time()
    print(t_end)

    if log:
        ts = time.time()
        f = open("tot_mp_%s_%4d.txt" % (nproc, ts), 'w')
        f.write("%2.2f, %2.2f, %2.2f, Total Daysheds:, %s, Times Convolved:, %s" % (t_start, output, t_end, no_daysheds, times_convolved))
        f.close()

if __name__ == "__main__":
    main(write=False, log=True, mp=False, nproc=16)
