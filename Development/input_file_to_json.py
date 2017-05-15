import json
import os

def read_input_file(input_file):

    def fetch(reader):
        line = reader.readline()
        return line.split("!")[0].strip()

    with open(input_file, 'r') as f:
        p = {
            "eco_or_dw": fetch(f),			                # eco or dw
            "start_count": int(fetch(f)),		        	# num_record of simulation start date, since 1/1/1961
            "startjul": int(fetch(f)),	        	    	# Simulation Start Julian date
            "endjul": int(fetch(f)),		             	# Simulation End Julian date
            "firstyear": int(fetch(f)),	        	    	# First year of simulation
            "firstmon": int(fetch(f)),	        	    	# First month of simulation
            "firstday": int(fetch(f)),		              	# First day of simulation
            "lastyear": int(fetch(f)),		             	# Last year of simulation
            "numberyears": int(fetch(f)),		        	# Number of years in simulation
            "ndates": int(fetch(f)),		    	        # Total # days in simulation
            "jul_dates": list(map(int, fetch(f).split())),  # Julian dates of simulation days
            "sdates": fetch(f),			                    # Actual date strings
            "chem": fetch(f),			                    # Chemical name
            "number_crop_ids": int(fetch(f)),               # Total # crops
            "cropdesired": list(map(int, fetch(f).split())),# Crop IDs
            "koc": float(fetch(f)),			                # Koc, read as mL/g
            "kflag": int(fetch(f)),			                # Koc=1, Kd=2
            "soil_halfLife": float(fetch(f)),               # Soil half life
            "appflag": int(fetch(f)),			            # Application by Crop Stage (1) or User-defined (2)
            "distribflag": int(fetch(f)),		        	# Application distribution flag (1=unif, 2=unif step, 3=triang)
            "cropstage": int(fetch(f)),		                # Crop stage for app (pl=1,emer=2,mat=3,harv=4) or 0
            "stagedays": int(fetch(f)),			            # Days after/before crop stage, or 0
            "stageflag": int(fetch(f)),			            # after(1) or before(2) crop stage, or 0
            "napps": int(fetch(f)),			                # Total Number of Applications, =0 cropstage app
            "twindow1": int(fetch(f)),			            # Application time window1 (d), applicable to Unif, Unif Step, Triangular
            "twindow2": int(fetch(f)),			            # Application time window2 (d), applicable to Unif Step
            "pct1": float(fetch(f)),			            # Percent of App during window1 (%), applicable to Unif, Unif Step
            "pct2": float(fetch(f)),			            # Percent of App during window2 (%), applicable to Unif Step
            "appnumrec_init": int(fetch(f)),			    # Application num_record, =0 cropstage app
            "appdate_init": int(fetch(f)),			        # Application Julian dates, =0 cropstage app
            "appmass_init": float(fetch(f)),                # Mass of each application (kg/ha, coverted to kg/m2 below)
            "appmethod_init": int(fetch(f)),			    # Application method (1=ground,2=foliar)
            "outtype": int(fetch(f)),			            # Output type (1=Daily,2=TimeAvg)
            "avgpd": int(fetch(f)),			                # Averaging period (days)
            "outputtype": int(fetch(f)),			        # Time-Avg Output Type (1=AvgConcs, 2=ToxExceed)
            "timeavg": int(fetch(f)),		            	# Time-Avg Conc Options Selected
            "threshold": int(fetch(f)),		            	# Threshold(ug/L)
            "thresoutput": int(fetch(f)),			        # Threshold Options Selected
            "outformat": int(fetch(f))			            # Output format (1=table,2=map,3=plot,4=download)
            }

        return p

def write_to_json(data):
    print("{\"inputs\": {)")
    for k, v in data.items():
        print("\t\"{}\": \"{}\",".format(k, v))
    print("}}")

def main():
    input_file = "SAM.inp"
    data = read_input_file(input_file)
    write_to_json(data)


main()