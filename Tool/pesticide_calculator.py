import os
import numpy as np

import read
import pesticide_functions as functions
import write


def pesticide_calculator(input_file, flow_file, scenario_dir, recipe_path, hydro_path, output_file, input_years,
                         process_benthic=True, process_erosion=True):

    # Read SAM input file
    input = read.input_file(input_file)

    # Find and assemble recipes
    recipe_files = read.recipes(recipe_path, input_years, scenario_dir, input.cropdesired)

    # Loop through recipes and corresponding flows listed in flow file
    for recipe_id, q, v, xc in read.flows(flow_file, input.dates):

        print(recipe_id)

        total_runoff_by_year, total_erosion_by_year = read.hydro(hydro_path, recipe_id, input_years, input.start_count)

        for year in input_years:

            # Initialize arrays for runoff and erosion totals
            # JCH - What are we doing with total_erosion_mass? Is this going into the daily output?
            total_runoff, total_erosion = total_runoff_by_year[year], total_erosion_by_year[year]
            total_runoff_mass, total_erosion_mass = np.zeros_like(total_runoff), np.zeros_like(total_erosion)

            # Loop through scenarios contained in the recipe
            scenarios = recipe_files[recipe_id][year]
            for scenario_file, area in scenarios:

                # Read scenario
                scenario = read.scenario(scenario_file, input)

                # Compute pesticide applications
                pesticide_mass_soil = functions.applications(input, scenario)

                # Determine the loading of pesticide into runoff and erosion - MMF added erosion
                runoff_mass, erosion_mass = functions.transport(pesticide_mass_soil, scenario)

                # Update runoff and erosion totals
                total_runoff_mass += runoff_mass * area
                total_erosion_mass += erosion_mass * area

            # Compute concentration in water
            total_flow, baseflow, total_conc, runoff_conc, daily_depth, aqconc_avg1, aqconc_avg2, aq1_store = \
                functions.waterbody_concentration(q, xc, total_runoff, runoff_mass, erosion_mass,
                                                  process_benthic, area_wb, daily_depth,
                                                  input.degradation_aqueous, input.koc)

            # Write daily output
            write.daily(output_file, input.dates, total_flow, baseflow, total_runoff, total_conc, runoff_conc,
                        total_runoff_mass, aqconc_avg1, aqconc_avg2, aq1_store)


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

    recipe_path = read.FilePath(recipe_dir, recipe_format)
    hydro_path = read.FilePath(hydro_dir, hydro_format)
    output_path = read.FilePath(output_dir, output_format)

    input_years = [2010, 2011, 2012, 2013]

    process_benthic = True
    process_erosion = True

    pesticide_calculator(input_file, flow_file, scenario_dir, recipe_path, hydro_path, output_path, input_years,
                         process_benthic, process_erosion)

if __name__ == "__main__":
    main()