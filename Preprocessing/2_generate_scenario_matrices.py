import os
import numpy as np
import pandas as pd


# ISSUES
# Most weather/crop ids in matrix don't have a match in crop_dates (IA - 302903/800888)
# No emergence dates in crop dates

class ScenarioMatrix(object):
    def __init__(self, state, years, combo_path, soil_path, crop_params, crop_dates, aggregation, met_data,
                 output_path, aggregate_soils):
        """ Reads tables, joins them together, and formats a scenario matrix for a state """
        from Tool.fields import combo_fields

        # Initialize variables
        self.state = state
        self.aggregate_soils = aggregate_soils
        self.crop_params = crop_params
        self.crop_dates = crop_dates[crop_dates.state == state]
        self.aggregation = aggregation
        self.met_data = met_data

        # Set local table paths
        soil_file_format = '{}_soil_aggregation.txt' if aggregate_soils else '{}_soil_gSSURGO_2016.txt'
        self.soil_path = os.path.join(soil_path, soil_file_format.format(state))
        self.combo_paths = [os.path.join(combo_path, "combo_{}_genclass{}.txt".format(state, year)) for year in years]
        self.output_path = os.path.join(output_path, "{}_scenarios_agg.txt".format(state))

        # Initialize matrix by reading and merging the combo files for each year for the state
        print("\tReading combo tables...")
        self.matrix = pd.concat([pd.read_csv(path, dtype='str')[combo_fields] for path in self.combo_paths], axis=0)

        # Join crop, soil, and weather tables
        print("\tProcessing crops...")
        self.process_crops()

        print("\tProcessing soils...")
        self.process_soils()

        print("\tProcessing weather...")
        self.process_weather()

        # Clean up matrix and write to file
        print("\tCleaning up and writing to file...")
        self.write_to_file()

    def process_soils(self):
        """ Join combos with soil data and assign curve numbers """
        from Tool.params import num_to_hsg

        # Join soil data with combos
        soil = pd.read_csv(self.soil_path, dtype='str')
        soil = soil.merge(self.aggregation, left_on='mukey', right_on='soilAggStr', how='left')
        self.matrix = self.matrix.merge(soil, left_on='soilagg', right_on='soilAggCode', how='left')

        # Process curve number
        num_to_hsg.update({2: "A", 4: "B", 6: "C"})  # A/D -> A, B/D -> B, C/D -> C
        self.matrix.hsg = self.matrix.hsg.astype(np.float32).fillna(-1).astype(np.int32).astype('str')  # JCH - temporary
        for num, hsg in num_to_hsg.items():
            self.matrix.loc[self.matrix.hsg == str(num), 'cn_ag'] = self.matrix['cn_ag_' + hsg]
            self.matrix.loc[self.matrix.hsg == str(num), 'cn_fallow'] = self.matrix['cn_fallow_' + hsg]

        # /D soils are evenly numbered, selected with hsg % 2 = 0
        non_cultivated_slash_d = ((self.matrix.cultivated == '0') & (np.float32(self.matrix.hsg) % 2 == 0))
        self.matrix.loc[non_cultivated_slash_d, ['cn_ag', 'cn_fallow']] = \
            self.matrix[['cn_ag_D', 'cn_fallow_D']].loc[non_cultivated_slash_d]

        # Deal with maximum rooting depth
        self.matrix.loc[self.matrix.root_zone_max < self.matrix.amxdr, 'amxdr'] = self.matrix.root_zone_max

    def process_crops(self):
        """ Join CDL-class-specific parameters to the table and add new rows for double cropped classes """
        from Tool.params import double_crops

        # Process double crops
        self.matrix['orig_cdl'] = self.matrix['cdl']
        self.matrix['overlay'] = 0  # Overlay crops aren't used to generate runoff in pesticide calculator
        all_new = []
        for old_crop, new_crops in double_crops.items():
            for i, new_crop in enumerate(new_crops):
                new_rows = self.matrix[self.matrix.orig_cdl == str(old_crop)].copy()
                new_rows['cdl'] = str(new_crop)
                new_rows['overlay'] = i
                all_new.append(new_rows)
        new_data = pd.concat(all_new, axis=0)
        self.matrix = pd.concat([self.matrix, new_data], axis=0).astype('str')

        # Join crop params to master table
        self.matrix = self.matrix.merge(self.crop_params, on='cdl', how='left')

        # Create a weather/crop ID and join crop data to master table
        self.matrix['weather-crop'] = self.matrix.weather + '-' + self.matrix['cdl']
        self.matrix = pd.merge(self.matrix, self.crop_dates, on='weather-crop', how='left')

    def process_weather(self):
        """ Read metfile and assign snowmelt factor (sfac) """
        self.matrix = self.matrix.merge(self.met_data, left_on='weather', right_on='stationID', how='left')
        self.matrix['sfac'] = 0.36
        self.matrix.loc[self.matrix.cdl.isin((60, 70, 140, 190)), 'sfac'] = .16

    def write_to_file(self):
        """ Apply finishing touches to matrix """
        from Tool.fields import crop_event_fields, scenario_matrix_fields

        for field in self.matrix.columns:
            print(field, self.matrix[field].isnull().sum())

        # Create a unique ID field
        soil_id_field = 'soilagg' if not self.aggregate_soils else 'soilAggStr'
        self.matrix['scenario'] = \
            self.state + self.matrix[soil_id_field] + 'st' + self.matrix.weather + 'cdl' + self.matrix.cdl

        # Flag bad rows (rows that Shelly would have deleted)
        # JCH - update this to be good
        self.matrix['bad'] = 0
        self.matrix.loc[pd.isnull(self.matrix.cokey), 'bad'] = 1
        self.matrix.loc[self.matrix.mukey == '0', 'bad'] = 2
        self.matrix.loc[self.matrix.cdl == '180', 'bad'] = 3
        self.matrix.loc[self.matrix.soilagg == '0', 'bad'] = 4
        self.matrix.loc[(1 > self.matrix.rainfall.astype(np.float32)) |
                        (4 < self.matrix.rainfall.astype(np.float32)), 'bad'] = 5
        self.matrix.loc[self.matrix.scenario == '', 'bad'] = 6
        self.matrix.loc[self.matrix.scenario.duplicated(), 'bad'] = 7
        self.matrix.loc[self.matrix.kwfact == '0', 'bad'] = 8

        # Fill incomplete data
        self.matrix.loc[:, crop_event_fields.new].fillna(1, inplace=True)
        self.matrix.loc[:, 'irr_type'].fillna(0, inplace=True)
        self.matrix.loc[:, ['cn_ag', 'cn_fallow']].fillna(0, inplace=True)

        # Trim to fields needed for output
        self.matrix = self.matrix[['bad'] + scenario_matrix_fields].reset_index(drop=True)

        # Write to file
        self.matrix.to_csv(self.output_path, index=False, float_format='%3.2f')


def read_tables(crop_params_path, crop_dates_path, aggregation_path, metfile_path):
    from Tool.fields import kurt_fields, crop_event_fields, crop_params_fields

    crop_dates_fields = crop_event_fields + kurt_fields

    # Read tables
    crop_params = pd.read_csv(crop_params_path, dtype='str').rename(columns=crop_params_fields.convert)
    crop_dates = \
        pd.read_csv(crop_dates_path, dtype='str').rename(columns=crop_dates_fields.convert)[crop_dates_fields.new]
    aggregation = pd.read_csv(aggregation_path, dtype='str')
    met_data = pd.read_csv(metfile_path, dtype='str')

    return crop_params, crop_dates, aggregation, met_data


def main():
    from Tool.params import states, nhd_states

    combo_path = os.path.join("..", "Input", "combos_CDLagg_121317")
    soil_path = os.path.join("..", "Output", "Soils")
    crop_params_path = os.path.join("..", "Input", "cdl_class_params.csv")
    crop_dates_path = os.path.join("..", "Input", "crop_genclass_120117.csv")
    aggregation_path = os.path.join("..", "Input", 'MASTER_SOIL_LUT_051917.csv')
    metfile_path = os.path.join("..", "Input", "met_data.csv")


    output_path = os.path.join("..", "Output", "ScenarioMatrices")

    aggregate_soils = True
    years = range(2012, 2016)

    # Read tables
    crop_params, crop_dates, aggregation, met_data = \
        read_tables(crop_params_path, crop_dates_path, aggregation_path, metfile_path)

    # Iterate states
    for state in sorted(nhd_states['07']):
        print("Generating scenario matrix for {}...".format(state))
        state_soil = \
            ScenarioMatrix(state, years, combo_path, soil_path, crop_params, crop_dates, aggregation, met_data,
                           output_path, aggregate_soils)


main()