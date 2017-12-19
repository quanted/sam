from .parameters import paths as p
from .functions import InputParams, Hydroregion, Scenarios, Recipes, Outputs, initialize


def pesticide_calculator(input_data):

    # Initialize file structure
    initialize()

    # Initialize parameters from front end
    inputs = InputParams(input_data)

    # Loop through all NHD regions included in selected runs
    for region_id in inputs.active_regions:

        # Load watershed topology maps and account for necessary files
        print("Processing hydroregion {}...".format(region_id))
        region = Hydroregion(region_id, inputs.sim_type, p.map_path, p.flow_dir, p.upstream_path, p.lakefile_path, p.geometry_path)

        # Simulate application of pesticide to all input scenarios
        print("Processing scenarios...")
        scenarios = Scenarios(inputs, region_id, p.input_scenario_path, region.active_reaches,
                              retain=r"..\bin\deux")  # retain=r"..\bin\une"

        # Initialize output object
        outputs = Outputs(inputs, scenarios.names, p.output_path, region.geometry, region.feature_type, demo_mode=True)

        # Cascade downstream processing watershed recipes and performing travel time analysis
        for year in [2011]:  # manual years

            print("Processing recipes for {}...".format(year))
            recipes = Recipes(inputs, outputs, year, region, scenarios, p.output_path, region.active_reaches)

            # Iterate through batches of reaches upstream of a reservoir
            for reaches, lake in region.cascade():

                # Process all the recipes in the batch
                recipes.process_recipes(reaches)

                # Modify all stored recipe data in the batch to simulate passage through reservoir
                recipes.burn_reservoir(lake, reaches)

        # Write output
        print("Writing output...")
        outputs.write_output()


def main(input_data=None):
    if input_data is None:
        from .chemicals import atrazine_demo as input_data
    pesticide_calculator(input_data)


if __name__ == "__main__":
    from .chemicals import atrazine_demo
    main(atrazine_demo)
