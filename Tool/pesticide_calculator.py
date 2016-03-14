import numpy as np

from Tool import read
from Tool import write
from Tool import pesticide_functions as functions


def pesticide_calculator(input_file, flow_file, scenario_dir, recipe_path, hydro_path, output_file, input_years,
                         process_benthic, process_erosion, write_daily_files, convolution):

    # Read SAM input file
    input = read.input_file(input_file)

    # Find and assemble recipes
    recipe_files = read.recipes(recipe_path, input_years, scenario_dir, input.cropdesired)

    # Loop through recipes and corresponding flows listed in flow file
    for recipe_id, q, v, l in read.flows(flow_file, input.dates, filter=recipe_files.keys()):

        if recipe_id in recipe_files:

            print(recipe_id)

            total_runoff_by_year, total_erosion_by_year = \
                read.hydro(hydro_path, recipe_id, input_years, input.start_count, process_erosion)

            for year in input_years:

                # Initialize arrays for runoff and erosion totals
                # JCH - What are we doing with total_erosion_mass? Is this going into the daily output?
                total_runoff = total_runoff_by_year[2010]  # JCH - Outputs match Fortran if we keep fixed at 2010
                total_runoff_mass = np.zeros_like(total_runoff)
                total_erosion = total_erosion_by_year[year] if process_erosion else None
                total_erosion_mass = np.zeros_like(total_erosion) if process_erosion else None

                # Loop through scenarios contained in the recipe
                scenarios = recipe_files[recipe_id][year]

                for scenario_file, area in scenarios:

                    # Read scenario
                    scenario = read.scenario(scenario_file, input, process_erosion)

                    # Compute pesticide applications
                    pesticide_mass_soil = functions.applications(input, scenario)

                    # Determine the loading of pesticide into runoff and erosion - MMF added erosion
                    runoff_mass, erosion_mass = functions.transport(pesticide_mass_soil, scenario, input, process_erosion)

                    # Update runoff and erosion totals
                    total_runoff_mass += runoff_mass * area
                    if process_erosion:
                        total_erosion_mass += erosion_mass * area

                # Compute concentration in water
                total_flow, baseflow, total_conc, runoff_conc, aqconc_avg1, aqconc_avg2, aq1_store = \
                    functions.waterbody_concentration(q, v, l, total_runoff, total_runoff_mass, total_erosion_mass,
                                                      process_benthic, input.degradation_aqueous, input.koc)

                # Write daily output
                if write_daily_files:
                    write.daily(output_file, recipe_id, year, input.dates, total_flow, baseflow, total_runoff,
                                total_conc, runoff_conc, total_runoff_mass, aqconc_avg1, aqconc_avg2, aq1_store)


def main():

    input_file = r"..\bin\MarkTwain\Inputs\SAM.inp"
    flow_file = r"..\bin\MarkTwain\Flows\region_07.csv"
    scenario_dir = r"..\bin\MarkTwain\Scenarios\Pickled"
    recipe_dir = r"..\bin\MarkTwain\Recipes"
    hydro_dir = r"..\bin\MarkTwain\Hydro"
    output_dir = r"..\bin\Outputs\Python"

    recipe_format = "nhd_recipe_(\d+?)_(\d{4}).txt"
    hydro_format = "{}_hydro.txt"
    output_format = "Eco_{}_{}_daily.out"

    recipe_path = read.FilePath(recipe_dir, recipe_format)
    hydro_path = read.FilePath(hydro_dir, hydro_format)
    output_path = read.FilePath(output_dir, output_format)

    input_years = [2010, 2011, 2012, 2013]

    process_benthic = False
    process_erosion = False

    write_daily_files = True  # JCH - Pending a more detailed development of output type specifications
    convolution = False

    pesticide_calculator(input_file, flow_file, scenario_dir, recipe_path, hydro_path, output_path, input_years,
                         process_benthic, process_erosion, write_daily_files, convolution)

if __name__ == "__main__":
    main()