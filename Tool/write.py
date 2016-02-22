import os
from collections import OrderedDict
import numpy as np


def daily(output_file, recipe_id, year, dates, total_flow, baseflow, total_runoff, total_conc, runoff_conc, runoff_mass,
          aqconc_avg1=np.array([]), aqconc_avg2=np.array([]), aq1store=np.array([])):

    # Create the output directory if it doesn't exist
    if not os.path.isdir(output_file.dir):
        os.mkdir(output_file.dir)
    output_file = output_file.format(recipe_id, year)

    # Create output
    output_fields = OrderedDict([
        ("Date", dates),                    # Dates
        ("Conc(ug/L)", total_conc),         # Total pesticide concentration
        ("RMass(kg)", runoff_mass),         # Mass of pesticide in runoff
        ("Runoff(m)", total_runoff),        # Quantity of water in runoff
        ("RConc(ug/L)", runoff_conc),       # Concentration of pesticide in runoff
        ("WC_Conc(kg/m3)", aqconc_avg1),    # Concentration of pesticide in water column
        ("Ben_Conc(kg/m3)", aqconc_avg2),   # Concentration of pesticide in benthic sediments
        ("PeakWC(kg/m3)", aq1store),        # Peak concentration in water column (JCH - is that right?)
        ("TotalFlow(m3)", total_flow),      # Total discharge of water
        ("Baseflow(m3)", baseflow)          # Baseflow
        ])

    # Filter out fields that don't have data
    for key, value in output_fields.items():
        if not (key == "Date" or value.any()):
            del output_fields[key]

    # Initialize header and field types.  Field type is automatically float for all but dates
    header = ",".join(output_fields.keys())
    out_values = np.array(list(output_fields.values())).T
    fmt = ["%s"] + ["%1.4e"] * (len(output_fields.values()) - 1)

    # Write to file
    np.savetxt(output_file, out_values, delimiter=",", fmt=fmt, header=header, newline="\n", comments="")

    if __name__ == "__main__":
        print("This is a library. Run pesticide_calculator.py or travel_time.py")