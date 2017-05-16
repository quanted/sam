from parameters import write_list, paths as p
from functions import InputFile, Hydroregion, ScenarioMatrices, RecipeMatrices

# Existential questions:
# Better inspections: check for (and report on) dates, missing data (scenarios, etc), other failures
# Benthic partitioning


def pesticide_calculator(input_data):

    # Initialize parameters from front end
    inputs = InputFile(input_data)

    # Load watershed topology maps and account for necessary files
    print("Initializing watershed...")
    region = Hydroregion(inputs, p.map_path, p.flow_dir, p.upstream_path, p.lakefile_path)

    # Simulate application of pesticide to all input scenarios
    print("Processing scenarios...")
    scenarios = ScenarioMatrices(inputs, p.input_scenario_path, retain="mtb_" + inputs.chemical_name)

    # Cascade downstream processing watershed recipes and performing travel time analysis
    for year in inputs.manual_years:
        print("Processing recipes for {}...".format(year))
        recipes = RecipeMatrices(inputs, year, region, scenarios, p.output_path, write_list=write_list)
        for reaches, lake in region.cascade():
            recipes.process_recipes(reaches)
            recipes.burn_reservoir(lake, reaches)

        recipes.write_time_series(4867727)

def main(input_data=None):
    if input_data is None:
        import chemicals
        input_data = chemicals.atrazine
    pesticide_calculator(input_data)


if __name__ == "__main__":
    from chemicals import chlorpyrifos, atrazine

    time_it = True
    if time_it:
        import cProfile
        for chemical in (atrazine, ):
            cProfile.run('main({})'.format(chemical))
    else:
        for chemical in (atrazine,):
            main(chemical)
