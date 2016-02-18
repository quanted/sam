import os
import numpy as np

# JCH - made benthic and water column output optional for now, in case we need to run without.
#       For diagnostic purposes now, 'family of SAM' later?
def daily(output_file, dates, total_flow, baseflow, total_runoff, total_conc, runoff_conc, runoff_mass,
          aqconc_avg1=None, aqconc_avg2=None, aq1store=None):

    # Create the output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Create output
    output_fields = {"Date":            dates,          # Dates
                     "Conc(ug/L)":      total_conc,     # Total pesticide concentration
                     "RMass(kg)":       runoff_mass,    # Mass of pesticide in runoff
                     "Runoff(m)":       total_runoff,   # Quantity of water in runoff
                     "RConc(ug/L)":     runoff_conc,    # Concentration of pesticide in runoff
                     "WC_Conc(kg/m3)":  aqconc_avg1,    # Concentration of pesticide in water column
                     "Ben_Conc(kg/m3)": aqconc_avg2,    # Concentration of pesticide in benthic sediments
                     "PeakWC(kg/m3)":   aq1store,       # Peak concentration in water column (JCH - is that right?)
                     "TotalFlow(m3)":   total_flow,     # Total discharge of water
                     "Baseflow(m3)":    baseflow}       # Baseflow

    # Filter out fields that don't have data
    output_fields = {field: values for field, values in output_fields.items() if values}

    # Initialize header and field types.  Field type is automatically float for all but dates
    header = ",".join(output_fields.keys())
    out_values = np.array(output_fields.values()).T
    fmt = ["%s"] + ["%1.4e"] * (len(output_fields.values()) - 1)

    # Write to file
    np.savetxt(output_file, out_values, delimiter=",", fmt=fmt, header=header, newline="\n", comments="")


if __name__ == "__main__":
    print("This is a library. Run pesticide_calculator.py or travel_time.py")