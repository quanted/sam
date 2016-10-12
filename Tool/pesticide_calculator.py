import Tool.functions as func
from Tool.parameters import paths as p


def pesticide_calculator(input_data):
    # Initialize parameters from front end
    inputs = func.InputFile(input_data)
    import os
    p.output_path.dir = os.path.join(p.output_path.dir, inputs.chemical_name)  # JCH - move
    print("Initializing watershed...")
    region = func.Hydroregion(
        inputs, p.recipe_path, p.map_path, p.scenario_dir, p.flow_dir, p.upstream_path, p.lakefile_path, p.lentics_path)

    print("Processing scenarios...")
    scenario_matrix = func.ScenarioMatrix(inputs, region, p.scenario_dir)
    for year in inputs.years:

        print("Processing recipes for {}...".format(year))
        recipe_matrix = func.RecipeMatrix(inputs, year, region, scenario_matrix, p.output_path)

        for reaches, lake in region.cascade():

            recipe_matrix.process_recipes(reaches, p.output_path)

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
        for chemical in (chlorpyrifos, atrazine):
            cProfile.run('main({})'.format(chemical))
    else:
        for chemical in (chlorpyrifos, atrazine):
            main(chemical)
