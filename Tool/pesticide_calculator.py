import os
import numpy as np

import read
import pesticide_functions as functions
import write


def pesticide_calculator(input_file, flow_file, scenario_dir, recipe_path, hydro_path, output_path, input_years):

    # Read SAM input file
    input = read.input_file(input_file)

    # Find and assemble recipes
    recipe_files = read.recipes(recipe_path, input_years, scenario_dir, input.cropdesired)

    # Loop through recipes and corresponding flows listed in flow file
    for recipe_id, q, v, xc in read.flows(flow_file, input.dates):

        print(recipe_id)

        total_runoff_by_year = read.hydro(hydro_path, recipe_id, input_years, input.start_count)

        for year in input_years:  # recipe_files[recipe_id]:

            scenarios = recipe_files[recipe_id][year]

            total_runoff = total_runoff_by_year[year]

            total_runoff_mass = np.zeros_like(total_runoff)  # Initializes an array to hold daily total runoff mass

            for scenario_file, area in scenarios:

                # Read scenario
                scenario = read.scenario(scenario_file, input.start_count)

                # Compute pesticide applications
                pesticide_mass_soil = \
                    functions.applications(input, scenario)

                # Determine the loading of pesticide into runoff
                runoff_mass = functions.transport(pesticide_mass_soil, scenario)

                # Update total runoff
                total_runoff_mass += runoff_mass * area

            # Compute concentration in water
            q_tot, baseflow, total_conc, runoff_conc = \
                functions.waterbody_concentration(q, xc, total_runoff, total_runoff_mass)

            # Write daily output
            write.daily(output_path.format(recipe_id, year), total_conc, runoff_conc, total_runoff_mass,
                        input.dates, q_tot, baseflow, total_runoff, year)


def main():

    input_file = r"T:\SAM\FortranToPy\Inputs\SAM.inp"
    flow_file = r"T:\SAM\FortranToPy\MarkTwain\MO_flows.csv"

    scenario_dir = r"T:\SAM\FortranToPy\MarkTwain\Scenarios\Pickled"
    recipe_dir = r"T:\SAM\FortranToPy\MarkTwain\Recipes"
    hydro_dir = r"T:\SAM\FortranToPy\MarkTwain\Hydro"
    output_dir = r"T:\SAM\Outputs\Python"

    recipe_format = "nhd_recipe_(\d+?)_(\d{4}).txt"
    hydro_format = "{}_hydro.txt"
    output_format = "Eco_{}_{}_daily.out"

    recipe_path = os.path.join(recipe_dir, recipe_format)
    hydro_path = os.path.join(hydro_dir, hydro_format)
    output_path = os.path.join(output_dir, output_format)

    input_years = [2010, 2011, 2012, 2013]

    pesticide_calculator(input_file, flow_file, scenario_dir, recipe_path, hydro_path, output_path, input_years)

if __name__ == "__main__":
    main()