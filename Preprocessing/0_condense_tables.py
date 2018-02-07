import numpy as np
import os
from utilities import read_dbf, read_gdb


class CustomSSURGO(object):
    def __init__(self, ssurgo_path, output_path, overwrite):
        self.in_folder = ssurgo_path
        self.out_folder = output_path
        self.overwrite = overwrite

        self._value_table = None

        # Initialize output folder if it doesn't exist
        if not os.path.exists(self.out_folder):
            os.makedirs(self.out_folder)

        # Extract each state
        for state in states:
            state_gdb = os.path.join(ssurgo_path, "gSSURGO_{}.gdb".format(state).lower())
            if os.path.exists(state_gdb):
                self.extract_tables(state_gdb, state)
            else:
                print("No SSURGO data found for {}".format(state))

    @property
    def value_table(self):
        from fields import valu1_fields
        if self._value_table is None:
            value_table = os.path.join(self.in_folder, "valu_fy2016.gdb")
            valu1 = read_gdb(value_table, "valu1", ['mukey'] + valu1_fields.old)
            valu1 = valu1.rename(columns=valu1_fields.convert)
            self._value_table = valu1
        return self._value_table

    def extract_tables(self, gdb, state):
        from fields import chorizon_fields, muaggatt_fields, component_fields
        out_file = os.path.join(self.out_folder, state)
        if self.overwrite or not all(map(os.path.exists, (out_file + '_horizon.csv', out_file + '_params.csv'))):
            print("Generating {}...".format(out_file))

            # Read horizon table
            chorizon_table = read_gdb(gdb, 'chorizon', ['cokey'] + chorizon_fields.old + ['desgnmaster'])
            chorizon_table = chorizon_table.sort_values(['cokey', 'hzdept_r'])

            # Load other tables and join
            component_table = read_gdb(gdb, 'component', ['mukey'] + component_fields.old)
            muaggatt_table = read_gdb(gdb, 'muaggatt', ['mukey'] + muaggatt_fields.old)
            state_data = component_table.merge(muaggatt_table, on='mukey').merge(self.value_table, on='mukey')

            # Rename fields
            chorizon_table = chorizon_table.rename(columns=chorizon_fields.convert)
            state_data = state_data.rename(columns=(muaggatt_fields + component_fields).convert)

            chorizon_table.to_csv(out_file + '_horizon.csv')
            state_data.to_csv(out_file + '_params.csv')


class CustomNHDPlus(object):
    def __init__(self, nhd_path, out_path, overwrite):
        self.nhd_path = nhd_path
        self.out_path = out_path
        self.overwrite = overwrite

        for region in nhd_regions:
            out_folder = out_path.format(region)
            if not os.path.exists(out_folder):
                os.makedirs(out_folder)
            self.extract_tables(region)

    def extract_tables(self, region):
        from fields import plus_flow_fields, gridcode_fields, flowline_fields, vaa_fields
        outfile = os.path.join(self.out_path, "region_{}.npz".format(region))
        if self.overwrite or not os.path.exists(outfile):
            print(outfile)
            region_path = self.nhd_path.format(region)
            table_paths = [os.path.join("NHDPlusAttributes", "PlusFlow.dbf"),
                           os.path.join("NHDSnapshot", "Hydrography", "NHDFlowline.dbf"),
                           os.path.join("NHDPlusAttributes", "PlusFlowlineVAA.dbf"),
                           os.path.join("NHDPlusCatchment", "featureidgridcode.dbf")]

            table_fields = [plus_flow_fields, flowline_fields, vaa_fields, gridcode_fields]
            hydro_table = None
            for table, fields in zip(table_paths, table_fields):
                new_table = \
                    read_dbf(os.path.join(region_path, table))[fields.old].rename(
                        columns=fields.convert)
                if hydro_table is None:
                    hydro_table = new_table
                else:
                    hydro_table = hydro_table.merge(new_table, on='comid')
            np.savez_compressed(outfile, data=hydro_table.as_matrix(), columns=hydro_table.columns.values)


def main():
    root_dir = os.path.join("P:\\", "GIS_Data")
    ssurgo_inpath = os.path.join(root_dir, "ssurgo", "gSSURGO_2016")
    ssurgo_outpath = os.path.join(root_dir, "SAM_2018", "bin", "Tables", "CustomSoil")
    nhd_inpath = os.path.join(root_dir, "NHDPlus2", "NHDPlusV2", "NHDPlus{}")
    nhd_outpath = os.path.join(root_dir, "SAM_2018", "bin", "Tables", "CustomNHD")

    process_ssurgo = True
    process_nhd = False
    overwrite = True

    if process_ssurgo:
        CustomSSURGO(ssurgo_inpath, ssurgo_outpath, overwrite)

    if process_nhd:
        CustomNHDPlus(nhd_inpath, nhd_outpath, overwrite)


nhd_regions = ['01', '02', '03N', '03S', '03W', '04', '05', '06', '07', '08', '09', '10L', '10U', '11', '12', '13',
               '14', '15', '16', '17', '18']

states = ['AL', 'AR', 'AZ', 'CA', 'CO', 'CT', 'DC', 'DE', 'FL', 'GA', 'IA', 'ID', 'IL', 'IN', 'KS', 'KY', 'LA', 'MA',
          'MD', 'ME', 'MI', 'MN', 'MO', 'MS', 'MT', 'NC', 'ND', 'NE', 'NH', 'NJ', 'NM', 'NV', 'NY', 'OH', 'OK', 'OR',
          'PA', 'RI', 'SC', 'SD', 'TN', 'TX', 'UT', 'VA', 'VT', 'WA', 'WI', 'WV', 'WY']

main()
