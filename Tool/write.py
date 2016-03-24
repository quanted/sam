import os
from collections import OrderedDict
import numpy as np

def daily(output_file, recipe_id, year, dates, total_flow, baseflow, total_runoff,
          runoff_mass, runoff_conc, aqconc_avg_wb, aqconc_avg1,aqconc_avg2, aq_peak1=np.array([])):

    # Create the output directory if it doesn't exist
    if not os.path.isdir(output_file.dir):
        os.mkdir(output_file.dir)
    output_file = output_file.format(recipe_id, year)

    # Create output
    raw_output = OrderedDict([
        ("Date", dates),                    # Dates
        ("TotalFlow(m3)", total_flow),      # Total discharge of water
        ("Baseflow(m3)", baseflow)          # Baseflow
        ("Runoff(m)", total_runoff),        # Runoff flow
        ("RMass(kg)", runoff_mass),         # Mass of pesticide in runoff
        ("RConc(ug/L)", runoff_conc),       # Concentration of pesticide in runoff
        ("WB_Conc(ug/L)", aqconc_avg_wb)    # Daily avg aqueous conc in water column of water body, no benthic
        ("WC_Conc(ug/L)", aqconc_avg1),     # Daily avg aqueous conc of pesticide in water column of water body
        ("Ben_Conc(ug/L)", aqconc_avg2),    # Daily avg aqueous conc of pesticide in benthic region of water body
        ("WC_Peak(ug/L)", aq_peak1),        # Daily peak aqueous conc of pesticide in water column of water body

        ])

    # Filter out fields that don't have data
    final_output = OrderedDict()
    for key, value in raw_output.items():
        if (key == "Date" or value.size):
            final_output[key] = value

    # Initialize header and field types.  Field type is automatically float for all but dates
    header = ",".join(final_output.keys())
    out_values = np.array(list(final_output.values())).T
    fmt = ["%s"] + ["%1.4e"] * (len(final_output.values()) - 1)

    # Write to file
    np.savetxt(output_file, out_values, delimiter=",", fmt=fmt, header=header, newline="\n", comments="")

    if __name__ == "__main__":
        print("This is a library. Run pesticide_calculator.py or travel_time.py")