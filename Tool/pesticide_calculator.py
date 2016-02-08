import os
import numpy as np
import read
import pesticide_functions as functions
import write     


<<<<<<< HEAD
def pesticide_calculator(input_file, scenario_dir, flow_file, recipe_path, hydro_path, output_path, input_years):

    # Read in hardwired parameters (set in read.py)
    delta_x, foliar_deg, washoff, soil_2cm, runoff_effic = read.pesticide_parameters()
=======
def pesticide_calculator(input_file, flow_file, scenario_dir, recipe_path, hydro_path, output_path, input_years):
>>>>>>> master

    # Read SAM input file
    input = read.input_file(input_file)

    # Find and assemble recipes
<<<<<<< HEAD
    recipe_files = read.recipes(recipe_path, input_years, scenario_dir, cropdesired)
=======
    recipe_files = read.recipes(recipe_path, input_years, scenario_dir, input.cropdesired)
>>>>>>> master

    # Loop through recipes and corresponding flows listed in flow file
    for recipe_id, q, v, xc in read.flows(flow_file, input.dates):

        print(recipe_id)

        total_runoff_by_year = read.hydro(hydro_path, recipe_id, input_years, input.start_count)

        for year in input_years:  # recipe_files[recipe_id]:

            scenarios = recipe_files[recipe_id][year]

<<<<<<< HEAD
            total_runoff = total_runoff_by_year[year]  # @@@ - total_runoff_by_year[year]
=======
            total_runoff = total_runoff_by_year[year]
>>>>>>> master

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

# MMF - Perform Water Body Calcs on both sets of runoff masses: (1) NHDbasin outlets, and (2) larger drainages...

# MMF - Addition of solute holding capacity and mass transfer functions needed for Water Body Calcs
            capacity1, capacity2, fw1, fw2, theta, sed_conv_factor, omega = \
                functions.soluteholdingcap(koc)

            # Compute concentration in water
            q_tot, baseflow, total_conc, runoff_conc = \
                functions.waterbody_concentration(q, xc, total_runoff, total_runoff_mass, degradation_aqueous, omega, theta)

            # Write daily output
<<<<<<< HEAD
            output.daily(output_path, recipe_id, total_conc, runoff_conc, total_runoff_mass, dates, q_tot, baseflow,
                         total_runoff, year)
=======
            write.daily(output_path.format(recipe_id, year), total_conc, runoff_conc, total_runoff_mass,
                        input.dates, q_tot, baseflow, total_runoff, year)
>>>>>>> master


def main():

# MMF - relative paths

    input_file = r"..\..\MarkTwain\Inputs\SAM.inp"
    flow_file =  r"..\..\MarkTwain\MO_flows.csv"

    scenario_dir = r"..\..\MarkTwain\Scenarios\Pickled"
    recipe_dir = r"..\..\MarkTwain\Recipes"
    hydro_dir = r"..\..\MarkTwain\Hydro"
    output_dir = r"..\..\Outputs\Python"

    recipe_format = "nhd_recipe_(\d+?)_(\d{4}).txt"
    hydro_format = "{}_hydro.txt"
    output_format = "Eco_{}_{}_daily.out"

    recipe_path = os.path.join(recipe_dir, recipe_format)
    hydro_path = os.path.join(hydro_dir, hydro_format)
    output_path = os.path.join(output_dir, output_format)

    input_years = [2010, 2011, 2012, 2013]

<<<<<<< HEAD
    pesticide_calculator(input_file, scenario_dir, flow_file, recipe_path, hydro_path, output_path, input_years)
=======
    pesticide_calculator(input_file, flow_file, scenario_dir, recipe_path, hydro_path, output_path, input_years)
>>>>>>> master

if __name__ == "__main__":
    main()