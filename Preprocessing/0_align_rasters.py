import arcpy
import os
from arcpy.sa import *

arcpy.CheckOutExtension("Spatial")


# Can this be de-ESRIfied?

def project_raster(in_path, out_path, sample_raster, overwrite):
    if overwrite or not os.path.exists(out_path):
        if arcpy.Exists(out_path):
            arcpy.Delete_management(out_path)
        else:
            print(out_path)
        raster = arcpy.Raster(in_path)
        arcpy.ProjectRaster_management(raster, out_path, sample_raster, "NEAREST")


def main():
    from Tool.params import nhd_states

    # Set input paths
    nhd_path = os.path.join(r"C:\\", "data", "NHDPlusV2", "NHDPlus{}", "NHDPlusCatchment", "cat")
    cdl_path = os.path.join(r"C:\\", "data", "CDL", "{0}_30m_cdls", "{0}_30m_cdls.img")
    met_path = \
        os.path.join(r"C:\\", "Users", "shelly", "Documents", "SAM", "weather_data", "weather_30m_US", "weather_30m")
    soil_path = os.path.join(r"C:\\", "data", "CustomSSURGO", "{0}", "{0}")

    # Set output paths
    projected_met_path = os.path.join("Output", "ProjectedLayers", "weather_grid")
    projected_cdl_path = os.path.join("Output", "ProjectedLayers", "cdl{}")
    projected_nhd_path = os.path.join("Output", "NHDCatchments", "region{}")

    years = range(2010, 2017)
    overwrite_nhd = True
    overwrite_cdl = False
    overwrite_met = False


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


main()
