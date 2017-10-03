atrazine = {'simulation_name': 'Mark Twain Atrazine 062217',
            'chemical_name': 'Atrazine',
            'region': 'mtb',
            'sim_date_start': '01/01/2000',
            'sim_date_end': '12/31/2015',
            'sim_type': 'eco',
            'applications': [
                {'refine': 'uniform', 'window1': 10, 'pct1': 100, 'offset': 3, 'crop': '10', 'method': 'foliar',
                 'effic': 1, 'window2': 0, 'event': 'emergence', 'pct2': 0, 'rate': 1.1},
                {'refine': 'uniform', 'window1': 10, 'pct1': 100, 'offset': 17, 'crop': '90', 'method': 'foliar',
                 'effic': 1, 'window2': 0, 'event': 'emergence', 'pct2': 0, 'rate': 1.1}],
            'aq_photolysis_hl': 168.0,
            'ben_metabolism_hl': 277.0,
            'hydrolysis_hl': 0.0,
            'soil_hl': 139.0,
            'wc_metabolism_hl': 277.0,
            'kd_flag': 1,
            'koc': 75.0,
            'endpoints': {
                'chronic_fw_inv': 60.0,
                'chronic_fw_fish': 0.5,
                'acute_em_fish': 1000.0,
                'overall_human': '',
                'overall_vasc_plant': 0,
                'acute_em_inv': 24.0,
                'overall_nonvasc_plant': 0,
                'chronic_vasc_plant': '',
                'chronic_em_fish': 0.5,
                'acute_fw_fish': 2650.0,
                'overall_em_fish': 0,
                'acute_nonvasc_plant': 1.0,
                'acute_human': 3.4,
                'chronic_em_inv': 80.0,
                'chronic_human': '',
                'overall_fw_fish': 0,
                'overall_em_inv': 0,
                'acute_vasc_plant': 4.6,
                'overall_fw_inv': 0,
                'acute_fw_inv': 360.0,
                'chronic_nonvasc_plant': ''}}


# This will no longer work because the format has changed.  Keeping it here so it can be mined for numbers
chlorpyrifos_old = \
    {"inputs":
         {"chemical_name": "chlorpyrifos",
          "applications":
          # crop,event,offset,window1,pct_applied1,window2,pct_applied2,rate,method,refine
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
          "sim_date_start": "2000-01-01",
          "sim_date_end": "2014-12-31",
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
