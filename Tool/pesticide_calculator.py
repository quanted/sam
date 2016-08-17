import numpy as np

from Tool import read, write
from Tool import pesticide_functions as functions
from Tool.parameters import paths as p


def pesticide_calculator(input_data):

    # Initialize parameters from front end
    inputs = read.InputParams(input_data)
    inputs.write_daily_files = True
    inputs.convolution = False

    # Locate and index recipe files by reach ID and year
    recipe_map = read.map_recipes(p.recipe_path, inputs.years)

    n_recipes = len(recipe_map.keys())

    # Create an array to transport output to convolution routine
    output_arrays = {year: functions.ConvolutionArray() for year in inputs.years}

    count = 0
    # Loop through recipes and corresponding flows listed in flow fil
    for recipe_id, flow in read.flows(p.flow_file, inputs.dates, filter=recipe_map.keys()):

        count +=1
        if not count % 100:
            print(count)

        # Read the recipe file that corresponds to the Recipe ID and construct a list of scenarios
        recipe = read.Recipe(recipe_map, recipe_id, p.scenario_dir, inputs)

        # Loop through each of the years of surface hydrology (runoff and erosion) in the scenario
        for year, runoff_and_erosion in read.hydro(p.hydro_path, recipe.id, inputs.hydro_date_offset, inputs.years,
                                                   process_erosion=True):

            # Initialize array for cumulative values of runoff (first row) and erosion (second row)
            transported_mass = np.zeros(runoff_and_erosion.shape)

            # Loop through each scenario for the catchment-year
            for scenario in recipe.scenarios[year]:
                # Compute pesticide application that winds up in soil
                pesticide_mass_soil = functions.applications(inputs, scenario)

                # Determine the loading of pesticide into runoff and erosion
                transported_mass += functions.transport(pesticide_mass_soil, scenario, inputs) * scenario.area

            # Compute concentration in water
            total_flow, baseflow, runoff_conc, aqconc_avg_wb, aqconc_avg, aq_peak = \
                functions.waterbody_concentration(flow, runoff_and_erosion, transported_mass,
                                                  inputs.process_benthic, inputs.degradation_aqueous, inputs.koc)

            # Write daily output
            if inputs.write_daily_files:
                write.daily(p.output_path, recipe_id, year, inputs.dates_str, total_flow, baseflow, runoff_and_erosion,
                            aqconc_avg_wb, runoff_conc, transported_mass, aqconc_avg, aq_peak, report=False)

            if inputs.convolution:

                if not output_arrays[year].initialized:
                    output_arrays[year].initialize(n_recipes, total_flow.size)

                output_arrays[year].update(recipe_id, transported_mass[0], runoff_and_erosion[0], baseflow)

    if inputs.convolution:
        for year, array in output_arrays.items():
            array.write_to_file(p.convolution_path.format(year))
            # time_of_travel(output_array)


def main(input_data=None):
    if input_data is None:
        input_data = {"inputs":
                          {"scenario_selection": "0",
                           "crop": "10 40 15 18",
                           "refine": "uniform_step",
                           "output_time_avg_conc": "1",
                           "application_rate": "1.3",
                           "crop_list_no": "10,40,15,18,140",
                           "output_avg_days": "4",
                           "workers": "16",
                           "crop_number": "4",
                           "chemical_name": "Custom",
                           "soil_hl": "417",             # Soil half life at 25 degC (atrazine)
                           "wc_metabolism_hl": "277",    # Water column metabolism half life at 25 degC (atrazine)
                           "ben_metabolism_hl": "588",   # Benthic metabolism half life at 25 degC (atrazine)
                           "aq_photolysis_hl": "168",    # Aqueous photolysis half life at 40 deg latitude (atrazine)
                           "hydrolysis_hl": "0",         # Hydrolysis half life (atrazine)
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
                           "region": "Mark Twain Basin",
                           "apps_per_year": "1",
                           "output_type": "2",
                           "refine_percent_applied2": "0",
                           "koc": "100",
                           "refine_percent_applied1": "100"},
                      "run_type": "single"}

    pesticide_calculator(input_data)


if __name__ == "__main__":
    time_it = False
    if time_it:
        import cProfile

        cProfile.run('main()')
    else:
        main()
