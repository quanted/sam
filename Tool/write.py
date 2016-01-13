import os
import numpy as np

def daily(output_path, reach, total_conc, runoff_conc, runoff_mass, dates, total_flow, baseflow,
          total_runoff, year=None):

    output_dir = os.path.dirname(output_path)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    output_file = output_path.format(reach, year)

    # Create output
    header = "Date,Conc(ug/L),RMass(kg),Runoff(m),RConc(ug/L),TotalFlow(m3),baseflow(m3)"
    fmt = ["%s","%f","%f","%f","%f","%f","%f"]
    out_values = np.array([dates, total_conc, runoff_mass, total_runoff, runoff_conc, total_flow, baseflow]).T

    np.savetxt(output_file, out_values, delimiter=",", fmt=fmt, header=header, newline="\n", comments="")


if __name__ == "__main__":
    print("This is a library. Run pesticide_calculator.py or travel_time.py")