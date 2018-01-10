import os
import numpy as np
import pandas as pd

from scipy import stats


class StateSoil(object):
    def __init__(self, state, ssurgo, out_path):
        from Tool.params import depth_bins

        # Initialize variables
        self.state = state
        self.ssurgo = ssurgo
        self.depth_bins = depth_bins
        self.out_path = out_path

        # Initialize paths
        self.unaggregated_outfile = os.path.join(self.out_path, self.state + '_soil_gSSURGO_2016.txt')
        self.aggregated_outfile = os.path.join(self.out_path, self.state + '_soil_aggregation.txt')
        self.aggregation_map_file = os.path.join(self.out_path, self.state + '_aggregation_map.txt')

        # Read data from SSURGO
        self.soil_table, self.soil_params = self.read_tables()

        print("\tSelecting components...")
        components = self.identify_components()

        print("\tDepth weighting..")
        self.depth_weighting(components)

        print("\tAdding soil attributes...")
        self.soil_attribution()

        # Aggregate soil
        print("\tAggregating soil...")
        self.aggregation_map, self.aggregated_soils = self.aggregate_soils()

        # Write to file
        print("\tWriting to file...")
        self.write_to_file()

    def aggregate_soils(self):
        from Tool.params import uslep_values, bins

        # Sort data into bins

        out_data = [self.soil_table.hsg_letter]
        for field, field_bins in bins.items():
            labels = [field[:2 if field == "slope" else 1] + str(i) for i in range(1, len(field_bins))]
            sliced = pd.cut(self.soil_table[field], field_bins, labels=labels, right=False).astype("str")
            out_data.append(sliced)
        soil_agg = pd.concat(out_data, axis=1)

        soil_agg['aggregation_key'] = \
            soil_agg['hsg_letter'] + soil_agg['slope'] + soil_agg['orgC_5'] + soil_agg['sand_5'] + soil_agg['clay_5']

        # Add USLEP and its adjustment (Nelson email 10/3/17)
        for i, val in enumerate(uslep_values):
            soil_agg.loc[soil_agg.slope == ('sl' + str(i + 1)), 'uslep'] = val

        # Add aggregation ID and USLEP to soil table
        both = pd.concat([self.soil_table, soil_agg[['aggregation_key', 'uslep']]], axis=1)
        aggregation_map = both[['mukey', 'aggregation_key']]

        # Group by aggregation ID and take the mean of all properties except HSG, which will use mode
        grouped = both.groupby('aggregation_key')
        averaged = grouped.mean().reset_index()
        hsg = grouped['hsg'].agg(lambda x: stats.mode(x)[0][0]).to_frame().reset_index()
        del averaged['hsg']
        aggregated_soils = averaged.merge(hsg, on='aggregation_key')

        return aggregation_map, aggregated_soils

    def depth_weighting(self, all_components):
        from Tool.fields import chorizon_fields, depth_fields

        # Loop through all components
        all_data = []
        depth_labels = range(1, self.depth_bins.size + 1)
        for i, component in enumerate(all_components):
            if i and not i % 1000:
                print("\t\tProcessed {} of {} components".format(i, len(all_components)))

            # Limit to current component
            component_table = self.soil_table[self.soil_table.cokey == component]
            if component_table.shape[0] > 1:

                # Get a valid kwfact from one of the top two layers
                top_row, second_row = component_table.iloc[0], component_table.iloc[1]
                kwfact = top_row.kwfact if top_row.kwfact != 9999 else second_row.kwfact
                if any((np.abs(kwfact) == 9999, top_row.hzdept_r not in (0, 1), top_row.desgnmaster == 'R')):
                    kwfact = np.nan  # bad component

                # Trim to horizons in top 100 cm and interpolate
                component_table = component_table[component_table['hzdept_r'] < 100]
                thickness = np.int16((component_table.hzdepb_r - component_table.hzdept_r).values)
                soil = np.repeat(component_table[chorizon_fields.new].as_matrix(), thickness, axis=0)[:100]
                soil = pd.DataFrame(soil, columns=chorizon_fields.new)

                # Group by horizon and compute averages
                soil['horizon'] = \
                    pd.cut(np.arange(soil.shape[0]), np.array([0] + list(self.depth_bins)) - 1, labels=depth_labels)
                averages = soil.groupby('horizon').mean()
                local = np.concatenate((np.array([component, kwfact]), averages.T.values.flatten()))
                all_data.append(local)
        self.soil_table = pd.DataFrame(data=np.array(all_data), columns=['cokey', 'kwfact'] + depth_fields)

    def soil_attribution(self):

        def calculate_uslels(df):
            from Tool.params import uslels_matrix
            row = (uslels_matrix.index.values < df.slope).sum()
            col = (uslels_matrix.columns.values < df.slope_length).sum()
            try:
                return uslels_matrix.iloc[row, col]
            except IndexError:
                return np.nan

        """ Merge soil table with params and get hydrologic soil group (hsg) """
        from Tool.params import hsg_to_num, num_to_hsg

        # Join soil table to params table
        self.soil_table = self.soil_table.merge(self.soil_params, on='cokey')

        # New HSG code - take 'max' of two versions of hsg
        self.soil_table['hsg'] = \
            self.soil_table[['hydgrpdcd', 'hydgrp']].applymap(
                lambda x: hsg_to_num.get(x)).max(axis=1).fillna(-1).astype(np.int32)
        self.soil_table['hsg_letter'] = self.soil_table['hsg'].map(num_to_hsg)

        # Calculate USLE LS
        self.soil_table['uslels'] = self.soil_table.apply(calculate_uslels, axis=1)

    def identify_components(self):
        """  Identify component to be used for each map unit """
        selection_fields = ['mukey', 'cokey', 'majcompflag', 'comppct_r']

        # Select major compoments (majcompflag)
        major_components = self.soil_params[selection_fields][self.soil_params.majcompflag == 'Yes']

        # Select major component with largest area (comppct)
        major_components.sort_values('comppct_r', ascending=False, inplace=True)
        selected_components = major_components[~major_components['mukey'].duplicated()]
        return selected_components.cokey.unique()

    def read_tables(self):
        from Tool.fields import chorizon_fields, component_fields

        # Read horizon table
        chorizon_table = self.ssurgo.fetch(self.state, 'chorizon').sort_values(['cokey', 'hzdept_r'])

        # Load other tables and join
        component = self.ssurgo.fetch(self.state, 'component')
        muaggatt = self.ssurgo.fetch(self.state, 'muaggatt')
        valu1 = self.ssurgo.fetch(self.state, 'valu1')
        soil_params = component.merge(muaggatt, on='mukey').merge(valu1, on='mukey')

        # Adjust data
        chorizon_table.loc[:, 'om_r'] /= 1.724
        chorizon_table.loc[:, ['wthirdbar_r', 'wfifteenbar_r']] /= 100.

        # Change field names
        chorizon_table.rename(columns=chorizon_fields.convert, inplace=True)
        soil_params.rename(columns=component_fields.convert, inplace=True)

        return chorizon_table, soil_params

    def write_to_file(self):
        from Tool.fields import soil_table_fields
        # Confine to output fields
        self.soil_table = self.soil_table[['mukey'] + soil_table_fields]
        self.aggregated_soils = self.aggregated_soils[['aggregation_key'] + soil_table_fields]

        # Write aggregation map
        self.aggregation_map.to_csv(self.aggregation_map_file, index=False)

        # Write aggregated
        self.soil_table.reset_index(drop=True).to_csv(self.unaggregated_outfile, index=False, float_format='%3.2f')

        # Write averaged
        self.aggregated_soils.reset_index(drop=True).to_csv(self.aggregated_outfile, index=False, float_format='%3.2f')


class SSURGOReader(object):
    def __init__(self, ssurgo_dir):
        self.path = ssurgo_dir
        self.valu1_path = os.path.join(ssurgo_dir, "valu1.txt")

    def fetch(self, state, table_name, index=None):
        if table_name != 'valu1':
            table_path = os.path.join(self.path, state, table_name + '.csv')
        else:
            table_path = os.path.join(self.path, 'valu1.txt')
        return pd.read_csv(table_path, index_col=index)


def main():
    from Tool.params import states, nhd_states

    # Set paths
    ssurgo_path = os.path.join("..", "..", "NationalData", "CustomSSURGO")
    out_path = os.path.join("..", "Output", "Soils")

    # Initialize SSURGO reader
    ssurgo = SSURGOReader(ssurgo_path)

    # Iterate through states
    for state in nhd_states['07']:
    #for state in ['RI']:
        print("Processing {} soil...".format(state))

        # Process state
        StateSoil(state, ssurgo, out_path)


main()
