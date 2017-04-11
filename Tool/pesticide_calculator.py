import functions as func
from parameters import write_list
from parameters import paths as p

# Existential questions:
# Recipe map takes a while to load
# Scenarios loading from memmap, but there's going to be a structural change (region to metfile?)
# Run area
# Output types
# Checks for completeness


def pesticide_calculator(input_data):

    # Initialize parameters from front end
    inputs = func.InputFile(input_data)

    # Load watershed topology maps and account for necessary files
    print("Initializing watershed...")
    region = func.Hydroregion(inputs, p.map_path, p.flow_dir, p.upstream_path, p.lakefile_path, p.input_scenario_path)

    print("Processing scenarios...")
    scenario_matrix = \
        func.ScenarioMatrix(inputs, region.scenario_memmap,
                            stored=r'..\bin\Preprocessed\mtb_{}.dat'.format(inputs.chemical_name), overwrite=False)

    region.years = [2010]  # JCH - temporary
    for year in region.years:
        print("Processing recipes for {}...".format(year))
        recipe_matrix = func.RecipeMatrix(inputs, year, region, scenario_matrix, p.output_path, write_list=write_list)
        for reaches, lake in region.cascade():
            recipe_matrix.process_recipes(reaches)
            recipe_matrix.time_of_travel(reaches, lake)


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
        for chemical in (chlorpyrifos,):
            cProfile.run('main({})'.format(chemical))
    else:
        for chemical in (chlorpyrifos,):
            main(chemical)
