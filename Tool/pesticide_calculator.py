from Tool.parameters import paths as p
from Tool.functions import InputParams, Hydroregion, Scenarios, Recipes, Outputs, initialize


def pesticide_calculator(input_data):

    # Initialize file structure
    initialize()

    # Initialize parameters from front end
    inputs = InputParams(input_data)

    # Loop through all NHD regions included in selected runs
    for region_id in inputs.active_regions:

        # Simulate application of pesticide to all input scenarios
        print("Processing scenarios...")
        scenarios = Scenarios(inputs, p.input_scenario_path)  # retain="onze"

        # Load watershed topology maps and account for necessary files
        print("Processing hydroregion {}...".format(region_id))
        region = Hydroregion(region_id, p.map_path, p.flow_dir, p.upstream_path, p.lakefile_path)

        # Initialize output object
        outputs = Outputs(inputs, region.active_reaches, scenarios.names, p.output_path)

        # Cascade downstream processing watershed recipes and performing travel time analysis
        for year in [2010]:  # manual years

            print("Processing recipes for {}...".format(year))
            recipes = Recipes(inputs, outputs, year, region, scenarios, p.output_path, inputs.active_reaches)

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
        from Tool.chemicals import atrazine_demo as input_data
    pesticide_calculator(input_data)


if __name__ == "__main__":
    from Tool.chemicals import atrazine_demo
    main(atrazine_demo)
