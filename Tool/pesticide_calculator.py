import os
import numpy as np

import read
import pesticide_functions as functions
import write     



def pesticide_calculator(input_file, scenario_dir, flow_file, recipe_path, hydro_path, output_path, input_years):

    # Read in hardwired parameters (set in read.py)
    delta_x, foliar_deg, washoff, soil_2cm, runoff_effic = read.pesticide_parameters()

    # Read SAM input file
    start_count, dates, ndates, cropdesired, koc, kflag, appflag, distribflag, cropstage, stagedays, stageflag, \
    app_windows, appnumrec, appmass, appmethod_init, degradation_aqueous = read.input_file(input_file)

    # Find and assemble recipes
    recipe_files = read.recipes(recipe_path, input_years, scenario_dir, cropdesired)

    # Loop through recipes and corresponding flows listed in flow file
    for recipe_id, q, _, xc in read.flows(flow_file, dates):

        print(recipe_id)

        total_runoff_by_year = read.hydro(hydro_path, recipe_id, input_years, start_count)

        for year in input_years:  # recipe_files[recipe_id]:

            scenarios = recipe_files[recipe_id][year]

# MMF  Sub in [year]? Already read in all 4 columns in hydro file (in read.hydro)
# Should we try to construct new runoff array based on using 2010-2013 hydro backwards in time? 
# so 2013 in 2013,...2010 in 2010,...2013 in 2009, 2012 in 2008, etc 
            total_runoff = total_runoff_by_year[2010]  # @@@ - total_runoff_by_year[year] (not using 2011-2013)

            total_runoff_mass = np.zeros_like(total_runoff)  # Initializes an array to hold daily total runoff mass

# MMF reading pickle scenario here?
            for scenario_file, area in scenarios:

                # Read scenario
                runoff, leaching, rain, plant_factor, soil_water, covmax, org_carbon, bulk_density = \
                    read.scenario(scenario_file, start_count)

                # Compute pesticide applications
                pesticide_mass_soil = \
                    functions.pesticide_applications(plant_factor, stageflag, distribflag, appnumrec, appmass,
                                                     appmethod_init, app_windows, stagedays, cropstage, rain, soil_2cm,
                                                     covmax, foliar_deg, washoff)

                # Determine the loading of pesticide into runoff
                runoff_mass = functions.transport(koc, org_carbon, bulk_density, degradation_aqueous, soil_water, delta_x,
                                                  kflag, runoff, leaching, runoff_effic, pesticide_mass_soil)

                # Update total runoff
                total_runoff_mass += runoff_mass * area

# MMF - INSERT TIME OF TRAVEL calc - find total runoff mass for larger drainages...
# then perform water body calcs on both outputs of NHDbasin water bodies
# and also 'cumulated larger drainage water bodies

# MMF - change Diagnostics to True?
###     time_of_travel(lake_file, upstream_file, sam_output_file, output_dir, output_format, diagnostics=False)


# MMF - add solute holding capacity, mass transfer functions
            capacity1, capacity2, fw1, fw2, theta, sed_conv_factor = \
                functions.soluteholdingcap()


            # Compute concentration in water
            q_tot, baseflow, total_conc, runoff_conc = \
                functions.waterbody_concentration(q, xc, total_runoff, total_runoff_mass)

            # Write daily output
            output.daily(output_path, recipe_id, total_conc, runoff_conc, total_runoff_mass, dates, q_tot, baseflow,
                         total_runoff, year)


def main():

# MMF - relative paths

    input_file = r"..\..\MarkTwain\Inputs\SAM.inp"   #r"T:\SAM\FortranToPy\Inputs\SAM.inp"
    flow_file =  r"..\..\MarkTwain\MO_flows.csv"     #r"T:\SAM\FortranToPy\MarkTwain\MO_flows.csv"

    scenario_dir = r"..\..\MarkTwain\Scenarios\Pickled"  #r"T:\SAM\FortranToPy\MarkTwain\Scenarios\Pickled"
    recipe_dir = r"..\..\MarkTwain\Recipes" #r"T:\SAM\FortranToPy\MarkTwain\Recipes"
    hydro_dir = r"..\..\MarkTwain\Hydro"    #r"T:\SAM\FortranToPy\MarkTwain\Hydro"
    output_dir = r"..\..\Outputs\Python"    #r"T:\SAM\Outputs\Python"

    recipe_format = "nhd_recipe_(\d+?)_(\d{4}).txt"
    hydro_format = "{}_hydro.txt"
    output_format = "Eco_{}_{}_daily.out"

    recipe_path = os.path.join(recipe_dir, recipe_format)
    hydro_path = os.path.join(hydro_dir, hydro_format)
    output_path = os.path.join(output_dir, output_format)

    input_years = [2010, 2011, 2012, 2013]

    pesticide_calculator(input_file, scenario_dir, flow_file, recipe_path, hydro_path, output_path, input_years)

if __name__ == "__main__":
    main()