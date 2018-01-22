import arcpy
import os
from arcpy.sa import *
from collections import OrderedDict

arcpy.CheckOutExtension("Spatial")


# Using hardwired local paths for now - big big datasets
# Can this be de-ESRIfied?

def project_raster(in_path, out_path, sample_raster, overwrite, resample_soil=False):
    if overwrite or not os.path.exists(out_path):
        if arcpy.Exists(out_path):
            arcpy.Delete_management(out_path)
        else:
            print(out_path)
            if not os.path.exists(os.path.dirname(out_path)):
                os.makedirs(os.path.dirname(out_path))
        raster = arcpy.Raster(in_path)
        if resample_soil:
            arcpy.ProjectRaster_management(raster, out_path, sample_raster, "NEAREST", cell_size=30)
        else:
            arcpy.ProjectRaster_management(raster, out_path, sample_raster, "NEAREST")


def main():
    #from Preprocessing.utilities import nhd_states, states

    # Set input paths
    nhd_path = os.path.join(r"C:\\", "data", "NHDPlusV2", "NHDPlus{}", "NHDPlusCatchment", "cat")
    cdl_path = os.path.join(r"C:\\", "data", "CDL", "{0}_30m_cdls", "{0}_30m_cdls.img")
    met_path = \
        os.path.join(r"C:\\", "Users", "shelly", "Documents", "SAM", "weather_data", "weather_30m_US", "weather_30m")
    soil_path = os.path.join(r"C:\\", "data", "CustomSSURGO", "{0}", "{0}")

    # Set output paths
    root_path = os.path.join(r"C:\\", "Users", "shelly", "Documents", "SAM", "Trip", "Output")
    projected_met_path = os.path.join(root_path, "ProjectedLayers", "weather_grid")
    projected_cdl_path = os.path.join(root_path, "ProjectedLayers", "cdl{}")
    projected_nhd_path = os.path.join(root_path, "NHDCatchments", "region{}")
    projected_soil_path = os.path.join(root_path, "SSURGO", "{}")

    years = range(2010, 2017)
    overwrite_nhd = False
    overwrite_cdl = False
    overwrite_met = False
    overwrite_ssurgo = False

    # Pick a sample SSURGO raster to use as template
    state_soil_raster = arcpy.Raster(soil_path.format('AL'))
    arcpy.env.snapRaster = state_soil_raster
    arcpy.env.outputCoordinateSystem = state_soil_raster

    # Project weather grid
    print("Projecting weather grid...")
    project_raster(met_path, projected_met_path, state_soil_raster, overwrite_met)

    # Project CDL rasters
    for year in years:
        print("Projecting CDL for {}...".format(year))
        project_raster(cdl_path.format(year), projected_cdl_path.format(year), state_soil_raster, overwrite_cdl)

    for region, states in nhd_states.items():
        print("Projecting catchments for region {}...".format(region))
        project_raster(nhd_path.format(region), projected_nhd_path.format(region), state_soil_raster, overwrite_nhd)

    for state in all_states:
        print("Projecting SSURGO for {}...".format(state))
        project_raster(soil_path.format(state), projected_soil_path.format(state), state_soil_raster, overwrite_ssurgo,
                       resample_soil=True)


# NHD regions and the states that overlap
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

# All states
all_states = sorted(set().union(*nhd_states.values()))

main()
