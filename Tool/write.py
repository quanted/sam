import os
import numpy as np

def daily(output_file, total_conc, runoff_conc, runoff_mass, dates, aqconc_avg1, aqconc_avg2, aq1store, total_flow, baseflow,
          total_runoff, year=None):

    output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Create output
    header = "Date,Conc(ug/L),RMass(kg),Runoff(m),RConc(ug/L),WC_Conc(kg/m3), Ben_Conc(kg/m3), PeakWC (kg/m3), TotalFlow(m3),baseflow(m3)"
    fmt = ["%s"] + ["%1.4e"] * 6
    out_values = np.array([dates, total_conc, runoff_mass, total_runoff, runoff_conc, aqconc_avg1, aqconc_avg2, aq1store, total_flow, baseflow]).T

    np.savetxt(output_file, out_values, delimiter=",", fmt=fmt, header=header, newline="\n", comments="")


if __name__ == "__main__":
    print("This is a library. Run pesticide_calculator.py or travel_time.py")