import arcpy
from arcpy.sa import Combine
import os
import csv

from collections import OrderedDict

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True
arcpy.env.tileSize = "128 128"


def write_to_csv(out_table, in_raster):
    header = [field.name for field in arcpy.ListFields()]
    with open(out_table, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter=',', lineterminator='\n')
        writer.writeheader()
        for row in arcpy.da.SearchCursor(in_raster):
            writer.writerow(dict(zip(header, row)))


def main():

    nhd_states = OrderedDict((('01', {"ME", "NH", "VT", "MA", "CT", "RI", "NY"}),
                              ('02', {"VT", "NY", "PA", "NJ", "MD", "DE", "WV", "DC", "VA"}),
                              ('03N', {"VA", "NC", "SC", "GA"}),
                              ('03S', {"FL", "GA"}),
                              ('03W', {"FL", "GA", "TN", "AL", "MS"}),
                              ('04', {"WI", "MN", "MI", "IL", "IN", "OH", "PA", "NY"}),
                              ('05', {"IL", "IN", "OH", "PA", "WV", "VA", "KY", "TN"}),
                              ('06', {"VA", "KY", "TN", "NC", "GA", "AL", "MS"}),
                              ('07', {"MN", "WI", "SD", "IA", "IL", "MO", "IN"}),
                              ('08', {"MO", "KY", "TN", "AR", "MS", "LA"}),
                              ('09', {"ND", "MN", "SD"}),
                              ('10U', {"MT", "ND", "WY", "SD", "MN", "NE", "IA"}),
                              ('10L', {"CO", "WY", "MN", "NE", "IA", "KS", "MO"}),
                              ('11', {"CO", "KS", "MO", "NM", "TX", "OK", "AR", "LA"}),
                              ('12', {"NM", "TX", "LA"}),
                              ('13', {"CO", "NM", "TX"}),
                              ('14', {"WY", "UT", "CO", "AZ", "NM"}),
                              ('15', {"NV", "UT", "AZ", "NM", "CA"}),
                              ('16', {"CA", "OR", "ID", "WY", "NV", "UT"}),
                              ('17', {"WA", "ID", "MT", "OR", "WY", "UT", "NV"}),
                              ('18', {"OR", "NV", "CA"})))

    nhd_path = "C:/Users/shelly/Documents/SAM/nhd_processing/state_raster_gdbs/{0}_nhdplus.gdb/{0}_mosaic"
    weather_path = "C:/Users/shelly/Documents/SAM/weather_data/weather_30m_states/redo/{}_weather30m"
    soil_path = "C:/Users/shelly/Documents/SAM/ssurgo_2016/soil_and_aggregation_oct2017/{}_agg.tif"
    mask_path = 'C:/Users/shelly/Documents/SAM/nhd_processing/mask_files'
    cdl_path = os.path.join("C:/", "data", "CDL", "{0}_30m_cdls", "{0}_30m_cdls.img")
    output_path = os.path.join("..", "NewTables", "Combos", "pre_recipe_{}_cdl{}.txt")

    years = ['2010', '2011', '2012', '2013', '2014', '2015', '2016']

    for region, states in nhd_states.items():

        for state in map(lambda x: x.lower(), states):

            print("Processing State: " + state)
            # Fix this when it works
            arcpy.env.mask = os.path.join(mask_path, state + '_boundary.shp')

            state_rasters = [weather_path.format(state), soil_path.format(state), nhd_path.format(state)]

            # Calculate the combos of weather, landcover, and soil
            for year in years:
                cdl_raster = cdl_path.format(year)

                layers = [arcpy.Raster(p) for p in state_rasters + [cdl_raster]]

                # Overlay
                combined = Combine(layers)

                # Clean up
                for layer in zip(*layers)[0]:
                    arcpy.Delete_management(layer)

                # Write to file
                write_to_csv(output_path.format(state, year), combined)

main()
