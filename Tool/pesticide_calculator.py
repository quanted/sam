import numpy as np
from Tool import read
from Tool import write
from Tool import pesticide_functions as functions


def pesticide_calculator(input_data):
    from Tool.parameters import data as d

    # Initialize parameter set from json input
    i = read.input_validate(input_data)

    # Find and assemble recipes
    recipe_files = read.recipes(d.recipe_path, i.input_years, d.scenario_dir, i.cropdesired)

    # Loop through recipes and corresponding flows listed in flow file
    for recipe_id, q, v, l in read.flows(d.flow_file, i.dates, filter=recipe_files.keys()):

        if recipe_id in recipe_files:

            print(recipe_id)

            total_runoff_by_year, total_erosion_by_year = \
                read.hydro(d.hydro_path, recipe_id, i.input_years, i.start_count, i.process_erosion)

            for year in i.input_years:

                # Initialize arrays for runoff and erosion totals
                # JCH - What are we doing with total_erosion_mass? Is this going into the daily output?
                total_runoff = total_runoff_by_year[2010]  # JCH - Outputs match Fortran if we keep fixed at 2010
                total_runoff_mass = np.zeros_like(total_runoff)
                total_erosion = total_erosion_by_year[year] if i.process_erosion else None
                total_erosion_mass = np.zeros_like(total_erosion) if i.process_erosion else None

                # Loop through scenarios contained in the recipe
                scenarios = recipe_files[recipe_id][year]

                for scenario_file, area in scenarios:

                    # Read scenario
                    scenario = read.scenario(scenario_file, input, i.process_erosion)

                    # Compute pesticide applications
                    pesticide_mass_soil = functions.applications(input, scenario)

                    # Determine the loading of pesticide into runoff and erosion - MMF added erosion
                    runoff_mass, erosion_mass = functions.transport(pesticide_mass_soil, scenario, input, i.process_erosion)

                    # Update runoff and erosion totals
                    total_runoff_mass += runoff_mass * area
                    if i.process_erosion:
                        total_erosion_mass += erosion_mass * area

                # Compute concentration in water
                total_flow, baseflow, avgconc_adj, runoff_conc, aqconc_avg1, aqconc_avg2, aq_conc1 = \
                    functions.waterbody_concentration(q, v, l, total_runoff, total_runoff_mass, total_erosion_mass,
                                                      i.process_benthic, i.degradation_aqueous, i.koc)

                # Write daily output
                if i.write_daily_files:
                    write.daily(i.output_file, recipe_id, year, i.dates, total_flow, baseflow, total_runoff,
                                avgconc_adj, runoff_conc, total_runoff_mass, aqconc_avg1, aqconc_avg2)


def main(input_data=None):

    if input_data is None:
        input_data = {
            "inputs": {
                "ndates": "5479",
                "cropdesired": "10, 40, 15, 18",
                "appmass_init": "1.32",
                "pct1": "0.0",
                "outtype": "1",
                "stageflag": "2",
                "firstmon": "1",
                "appmethod_init": "1",
                "chem": "atrazine",
                "appnumrec_init": "0",
                "outformat": "1",
                "number_crop_ids": "4",
                "lastyear": "2014",
                "avgpd": "4",
                "start_count": "14245",
                "distribflag": "1",
                "appflag": "1",
                "startjul": "36525",
                "firstday": "1",
                "firstyear": "2000",
                "twindow2": "0",
                "napps": "1",
                "eco_or_dw": "eco",
                "threshold": "0",
                "stagedays": "14",
                "appdate_init": "0",
                "endjul": "42003",
                "outputtype": "0",
                "run_type": "single",
                "input_years": "2010 2011 2012 2013",
                "process_benthic": "False",
                "process_erosion": "False",
                "write_daily_files": "True",
                "convolution": "False"
            }
        }

    pesticide_calculator(input_data)

if __name__ == "__main__":
    main()
