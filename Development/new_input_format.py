    data_map = \
             {"scenario_selection": int,
              "crop": lambda x: list(map(int, x.split())),
              "refine": str,
              "output_time_avg_conc": int,
              "application_rate": float,
              "crop_list_no": lambda x: list(map(int, x.split(","))),   # Crops used (cropdesired)
              "output_avg_days": int,
              "workers": int,
              "crop_number": int,
              "chemical_name": str,
              "soil_metabolism_hl": int,                                # Soil metabolism halflife (soil_halfLife)
              "refine_time_window2": int,                               # twindow2
              "refine_time_window1": int,                               # twindow1
              "coefficient": int,                                       #
              "sim_date_1stapp": parse_date,                            # Start date
              "output_tox_value": int,
              "output_format": int,
              "sim_date_start": parse_date,
              "sim_type": str,
              "output_time_avg_option": int,
              "output_tox_thres_exceed": int,
              "processes": int,
              "sim_date_end": parse_date,                               # End date (start date + ndates)
              "application_method": int,                                # Application method (appmethod)
              "region": str,
              "apps_per_year": int,
              "output_type": int,
              "refine_percent_applied2": int,
              "koc": int,
              "refine_percent_applied1": int,
              "run_type": str,
              "input_years": lambda x: list(map(int, x.split())),
              "process_benthic": lambda x: x == "True",
              "process_erosion": lambda x: x == "True",
              "write_daily_files": lambda x: x == "True",
              "convolution": lambda x: x == "True"}