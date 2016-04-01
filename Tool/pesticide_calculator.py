import numpy as np
from Tool import pesticide_functions as functions
from Tool import read
from Tool import write
from Tool.parameters import paths as p


def pesticide_calculator(input_data):

    # Initialize parameter set from json input
    inputs = read.InputParams(input_data)

    # Locate and index recipe files by reach ID and year
    recipe_map = read.map_recipes(p.recipe_path, inputs.years)

    # Loop through recipes and corresponding flows listed in flow file
    count = 0
    for recipe_id, flow in read.flows(p.flow_file, inputs.dates, filter=recipe_map.keys()):
        count += 1
        if count > 10:
            abra

        print(recipe_id)

        recipe = read.Recipe(recipe_map, recipe_id, p.scenario_dir, inputs)

        for year, runoff_and_erosion in \
                read.hydro(p.hydro_path, recipe.id, inputs.date_offset, inputs.years, inputs.process_erosion):

            transported_mass = np.zeros_like(runoff_and_erosion)
            for scenario in recipe.scenarios[year]:

                # Compute pesticide applications
                pesticide_mass_soil = functions.applications(inputs, scenario)

                # Determine the loading of pesticide into runoff and erosion - MMF added erosion
                transported_mass += functions.transport(pesticide_mass_soil, scenario, inputs) * scenario.area

            # Compute concentration in water
            total_flow, baseflow, runoff_conc, aqconc_avg_wb, aqconc_avg1, aqconc_avg2, aq_peak1 = \
                functions.waterbody_concentration(flow, runoff_and_erosion, transported_mass,
                                                  inputs.process_benthic, inputs.degradation_aqueous, inputs.koc)

            # Write daily output
            if inputs.write_daily_files:
                write.daily(p.output_path, recipe_id, year, inputs.dates, total_flow, baseflow, runoff_and_erosion[0],
                            aqconc_avg_wb, runoff_conc, transported_mass[0], aqconc_avg1, aqconc_avg2, aq_peak1)

def main(input_data=None):

    if input_data is None:

        input_data = {"inputs": 
                       {"scenario_selection": "0",
                        "crop": "10 40 15 18",
                        "refine": "uniform_step",
                        "output_time_avg_conc": "1",
                        "application_rate": "1.3",
                        "crop_list_no": "10,40,15,18",
                        "output_avg_days": "4",
                        "workers": "16",
                        "crop_number": "4",
                        "chemical_name": "Custom",
                        "soil_metabolism_hl": "123",
                        "refine_time_window2": "0",
                        "refine_time_window1": "50",
                        "coefficient": "1",
                        "sim_date_1stapp": "04/20/1984",
                        "output_tox_value": "4",
                        "output_format": "3",
                        "sim_date_start": "01/01/2000",
                        "sim_type": "eco",
                        "output_time_avg_option": "2",
                        "output_tox_thres_exceed": "1",
                        "processes": "1",
                        "sim_date_end": "12/31/2014",
                        "application_method": "1",
                        "region": "Ohio Valley",
                        "apps_per_year": "1",
                        "output_type": "2",
                        "refine_percent_applied2": "50",
                        "koc": "100",
                        "refine_percent_applied1": "50"},
                        "run_type": "single"}

    pesticide_calculator(input_data)

if __name__ == "__main__":
    time_it = False
    if time_it:
        import cProfile
        cProfile.run('main()')
    else:
        main()