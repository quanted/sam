from Tool import read
from Tool.parameters import paths
from Tool import travel_time_functions as functions


# To do:
# * benthic transport
# * migrate read functions to read.py?

def time_of_travel(input_data):

    # Object containing input data parameters
    inputs = read.InputParams(input_data, mode="TimeOfTravel")

    # Output on/off. Included primarily for diagnostic purposes to limit output
    write_to_file = False

    # Structure which generates and stores impulse response functions
    irf = functions.ImpulseResponse()

    # Initializes a time of travel run for an NHD region
    region = functions.Region(inputs, paths, irf, write_to_file)

    # Loop through all reservoirs in a downstream cascade
    for reaches, lake in region.cascade():

        # Process reaches upstream of the reservoir
        for reach in reaches:
            reach.process()

        # Process the reservoir
        if lake:  # Last batch of reaches will not have an associated reservoir
            lake.process()


def main(input_data=None):
    if input_data is None:
        input_data = {"inputs":
                          {"mode": "convolved",
                           "region": "07",
                           "sam_output_id": "dummy_region_07",
                           },
                      "run_type": "single"}

    time_of_travel(input_data)


def batch_compare():
    for mode in ("unconvolved", "aggregated"):
        try:
            main({"inputs": {"mode": mode, "region": "07", "sam_output_id": "dummy_region_07"}})
        except Exception as e:
            print(e)


if __name__ == "__main__":
    time_it = True
    batch_compare_on = False
    if time_it:
        import cProfile
        cProfile.run('main()')
    elif batch_compare_on:
        batch_compare()
    else:
        main()
