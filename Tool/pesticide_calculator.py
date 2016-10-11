from RestructuredTool import functions as func
from RestructuredTool.parameters import paths as p


def pesticide_calculator(input_data):
    # Initialize parameters from front end
    inputs = func.InputFile(input_data)

    print("Initializing watershed...")
    region = func.Hydroregion(
        inputs, p.recipe_path, p.map_path, p.scenario_dir, p.flow_dir, p.upstream_path, p.lakefile_path, p.lentics_path)

    print("Processing scenarios...")
    scenario_matrix = func.ScenarioMatrix(inputs, region, p.scenario_dir)
    inputs.years = [2010]
    for year in inputs.years:

        print("Processing recipes for {}...".format(year))
        recipe_matrix = func.RecipeMatrix(inputs, region, scenario_matrix, p.output_path)

        for reaches, lake in region.cascade():

            recipe_matrix.process_recipes(reaches, year)

            recipe_matrix.time_of_travel(reaches, lake, inputs.convolution_mode)


def main(input_data=None):
    if input_data is None:
        import chemicals
        input_data = chemicals.chlorpyrifos
    pesticide_calculator(input_data)


from chemicals import chlorpyrifos, atrazine

if __name__ == "__main__":
    time_it = True
    if time_it:
        import cProfile
        cProfile.run('main()')
    else:
        main(chlorpyrifos)
