import os
from collections import OrderedDict
import numpy as np
import pandas as pd


def daily(output_file, recipe_id, year, dates, total_flow, baseflow, runoff_and_erosion, aqconc_avg_wb, runoff_conc,
          transported_mass, aqconc_avg, aq_peak):

    # Create the output directory if it doesn't exist
    if not os.path.isdir(output_file.dir):
        os.mkdir(output_file.dir)
    output_file = output_file.format(recipe_id, year)

    # Modify dates

    total_runoff, total_erosion = runoff_and_erosion
    runoff_mass, erosion_mass = transported_mass

    # Create output
    raw_output = [
        ("TotalFlow(m3)", total_flow),  # Total discharge of water
        ("Baseflow(m3)", baseflow),  # Baseflow
        ("Runoff(m)", total_runoff),  # Runoff flow
        ("RMass(kg)", runoff_mass),  # Mass of pesticide in runoff
        ("RConc(ug/L)", runoff_conc),  # Concentration of pesticide in runoff
        ("Conc(ug/L)", aqconc_avg_wb),
        # Daily avg aqueous conc in water column of water body, no benthic  # JCH - to be called WB_Conc eventually
        ("WC_Conc(ug/L)", aqconc_avg[0]),  # Daily avg aqueous conc of pesticide in water column of water body
        ("Ben_Conc(ug/L)", aqconc_avg[1]),  # Daily avg aqueous conc of pesticide in benthic region of water body
        ("WC_Peak(ug/L)", aq_peak[0]),  # Daily peak aqueous conc of pesticide in water column of water body
    ]

    # Filter out fields that don't have data
    raw_output = filter(lambda x: x[1].size > 1, raw_output)  # JCH - these should all be uniform arrays eventually

    header, out_fields = zip(*raw_output)
    output_array = np.vstack(out_fields).T
    df = pd.DataFrame(output_array, dates, header, dtype=float)

    df.to_csv(output_file, index_label="Date")

    if __name__ == "__main__":
        print("This is a library. Run pesticide_calculator.py or travel_time.py")


def daily_tot(output_path, sam_output_id, reach_id, mode, dates,
              total_flow, mass_and_runoff, baseflow, conc, runoff_conc,
              aqconc_avg=None, aq_peak=None):
    # Create the output directory if it doesn't exist
    if not os.path.isdir(output_path.dir):
        os.mkdir(output_path.dir)

    output_file = output_path.format(sam_output_id, reach_id, mode)

    runoff_mass, total_runoff = mass_and_runoff

    # Create output
    raw_output = OrderedDict([
        ("Date", dates),  # Dates
        ("TotalFlow(m3)", total_flow),  # Total discharge of water
        ("Baseflow(m3)", baseflow),  # Baseflow
        ("Runoff(m)", total_runoff),  # Runoff flow
        ("RMass(kg)", runoff_mass),  # Mass of pesticide in runoff
        ("RConc(ug/L)", runoff_conc),  # Concentration of pesticide in runoff
        ("Conc(ug/L)",
         conc)])  # Daily avg aqueous conc in water column of water body, no benthic  # JCH - to be called WB_Conc eventually
    if aqconc_avg:
        raw_output.update(OrderedDict([
            ("WC_Conc(ug/L)", aqconc_avg[0]),  # Daily avg aqueous conc of pesticide in water column of water body
            ("Ben_Conc(ug/L)", aqconc_avg[1])]))  # Daily avg aqueous conc of pesticide in benthic region of water body
    if aq_peak:
        raw_output.update(OrderedDict([
            ("WC_Peak(ug/L)", aq_peak[0])]))  # Daily peak aqueous conc of pesticide in water column of water body

    # Initialize header and field types.  Field type is automatically float for all but dates
    header = ",".join(raw_output.keys())
    out_values = np.array(list(raw_output.values())).T
    fmt = ["%s"] + ["%1.4e"] * (len(raw_output.values()) - 1)

    # Write to file
    np.savetxt(output_file, out_values, delimiter=",", fmt=fmt, header=header, newline="\n", comments="")


if __name__ == "__main__":
    print("This is a library. Run pesticide_calculator.py or travel_time.py")