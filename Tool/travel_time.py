from Tool import read
from Tool.parameters import paths as p
from Tool import travel_time_functions as functions

# To do:
# * stream_calc = 0
# * benthic transport
# * migrate read functions to read.py?


def time_of_travel(input_data):

    inputs = read.InputParams(input_data, mode="TimeOfTravel")

    region = functions.Region(inputs.mode, inputs.region, inputs.sam_output_id,
                              p.lakefile_path, p.lentics_path, p.sam_output_path, p.upstream_path)

    irf = functions.ImpulseResponse()

    write_to_file = True  # Included primarily for diagnostic purposes to limit output

    # Loop through all reservoirs in cascading order
    for waterbody in region.cascade():

        if type(waterbody) == functions.Reach:
            waterbody.process(irf, p.tot_output_path, inputs.sam_output_id, write_to_file)

        elif type(waterbody) == functions.Reservoir:
            waterbody.process(irf)


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
    time_it = False
    batch_compare_on = True
    if time_it:
        import cProfile
        cProfile.run('main()')
    elif batch_compare_on:
        batch_compare()
    else:
        main()
