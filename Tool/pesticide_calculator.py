from .parameters import nhd_regions, write_list, paths as p
from .functions import InputFile, Hydroregion, Scenarios, Recipes, confine_regions


def pesticide_calculator(input_data):

    # Initialize parameters from front end
    inputs = InputFile(input_data)

    # Simulate application of pesticide to all input scenarios
    print("Processing scenarios...")
    scenarios = Scenarios(inputs, p.input_scenario_path, retain="six")  # retain="six"

    # Loop through all NHD regions included in selected runs
    for region in confine_regions(nhd_regions, p.map_path, write_list):

        # Load watershed topology maps and account for necessary files
        print("Processing hydroregion {}...".format(region))
        region = Hydroregion(region, p.map_path, p.flow_dir, p.upstream_path, p.lakefile_path)

        # Cascade downstream processing watershed recipes and performing travel time analysis
        for year in [2010]:  # manual years

            print("Processing recipes for {}...".format(year))
            recipes = Recipes(inputs, year, region, scenarios, p.output_path, write_list=write_list)

            # Iterate through batches of reaches upstream of a reservoir
            for reaches, lake in region.cascade():

                # Process all the recipes in the batch
                recipes.process_recipes(reaches)

                # Modify all stored recipe data in the batch to simulate passage through reservoir
                recipes.burn_reservoir(lake, reaches)


        # Write output
        recipes.write_output()


def main(input_data=None):
    if input_data is None:
        from .chemicals import atrazine as input_data
    pesticide_calculator(input_data)


if __name__ == "__main__":
    from .chemicals import atrazine
    main(atrazine)
