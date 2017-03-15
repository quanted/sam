import os
import functions as func
from parameters import paths as p

# Existential questions:
# Recipe map takes a while to load
# Scenarios loading from memmap, but there's going to be a structural change (region to metfile?)

def pesticide_calculator(input_data):
    # Initialize parameters from front end
    inputs = func.InputFile(input_data)

    p.output_path.dir = os.path.join(p.output_path.dir, inputs.chemical_name)  # JCH - move

    print("Initializing watershed...")
    region = func.Hydroregion(inputs, p.map_path, p.flow_dir, p.upstream_path, p.lakefile_path, p.lentics_path,
                              scenario_memmap=r"..\..\..\Preprocessed\Scenarios\mtb_test")

    print("Processing scenarios...")
    scenario_matrix = func.ScenarioMatrix(inputs, region, region.scenario_memmap,
                                          stored=r'..\..\..\Preprocessed\mtb_new.dat', overwrite_stored=False)

    region.years = [2010]
    for year in region.years:
        print("Processing recipes for {}...".format(year))
        recipe_matrix = func.RecipeMatrix(inputs, year, region, scenario_matrix, p.output_path, filter=[5042400])
        for reaches, lake in region.cascade():
            recipe_matrix.process_recipes(reaches, p.output_path)
            if any(inputs.convolution_mode):
                recipe_matrix.time_of_travel(reaches, lake, inputs.convolution_mode, p.output_path)


def main(input_data=None):
    if input_data is None:
        import chemicals

        input_data = chemicals.chlorpyrifos
    pesticide_calculator(input_data)


if __name__ == "__main__":
    from chemicals import chlorpyrifos, atrazine

    time_it = False
    if time_it:
        import cProfile

        for chemical in (atrazine,):
            cProfile.run('main({})'.format(chemical))
    else:
        for chemical in (atrazine,):
            main(chemical)
