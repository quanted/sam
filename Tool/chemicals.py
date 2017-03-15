chlorpyrifos = \
    {"inputs":
         {"chemical_name": "chlorpyrifos",
          "region": "07",
          "applications":
          #crop,event,offset,window1,pct_applied1,window2,pct_applied2,rate,method,refine
              "10,emergence,3,10,100,0,0,1.1,foliar,uniform_step\n"
              "10,emergence,23,10,100,0,0,1.1,foliar,uniform_step\n"
              "40,emergence,-7,10,100,0,0,1.1,ground,uniform_step\n"
              "40,emergence,3,10,100,0,0,1.1,foliar,uniform_step\n"
              "70,harvest,14,10,100,0,0,2.2,foliar,uniform_step\n"
              "70,maturity,14,10,100,0,0,2.2,foliar,uniform_step\n"
              "50,emergence,-7,10,100,0,0,1.1,ground,uniform_step\n"
              "50,emergence,21,10,100,0,0,1.1,foliar,uniform_step\n"
              "90,emergence,7,10,100,0,0,1.1,foliar,uniform_step\n"
              "90,emergence,17,10,100,0,0,1.1,foliar,uniform_step\n"
              "90,emergence,27,10,100,0,0,1.1,foliar,uniform_step\n"
              "110,emergence,-10,10,100,0,0,1.1,ground,uniform_step\n"
              "110,emergence,21,10,100,0,0,1.1,foliar,uniform_step",
          "soil_hl": "170",
          "wc_metabolism_hl": "91",
          "ben_metabolism_hl": "203",
          "aq_photolysis_hl": "30",
          "hydrolysis_hl": "0",
          "kd_flag": "1",
          "koc": "6040",
          "sim_date_start": "01/01/2000",
          "sim_date_end": "12/31/2014",
          "output_type": "2",
          "output_time_avg_conc": "1",
          "output_avg_days": "4",
          "output_tox_value": "4",
          "output_format": "3",
          "output_time_avg_option": "2",
          "output_tox_thres_exceed": "1",
          "workers": "1",
          "processes": "1"
          },
     "run_type": "single"}

atrazine = \
    {"inputs":
         {"chemical_name": "atrazine",
          "region": "07",
          "applications":
              #crop,event,offset,window1,pct_applied1,window2,pct_applied2,rate,method,refine
              "10,plant,14,50,100,0,0,1.32,ground,uniform_step",
          "soil_hl": "139",
          "wc_metabolism_hl": "277",
          "ben_metabolism_hl": "277",
          "aq_photolysis_hl": "168",
          "hydrolysis_hl": "0",
          "kd_flag": "1",
          "koc": "75",
          "sim_date_start": "01/01/2000",
          "sim_date_end": "12/31/2014",
          "output_type": "2",
          "output_time_avg_conc": "1",
          "output_avg_days": "4",
          "output_tox_value": "4",
          "output_format": "3",
          "output_time_avg_option": "2",
          "output_tox_thres_exceed": "1",
          "workers": "1",
          "processes": "1"
          },
     "run_type": "single"}
