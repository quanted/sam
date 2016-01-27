import numpy as np



def qa_applications(stageflag, distribflag, appnumrec, twindow1, stagedays, appmass, appdate, startjul, appmethod,
                         pct1, pct2, twindow2, cropstage, plant_factor, start_count, rain, ndates,
                         soil_distribution_2cm, covmax, foliar_degradation, washoff_coeff):

    application_dates = get_dates(plant_factor, start_count)
    twindow2 = 10
    pct1, pct2 = 60, 40
    for distribflag in (1, 2, 3):
        for cropstage in (1, 2, 3, 4):
            for stageflag in (1, 2):
                new_appmass, new_appnumrec, new_appdate, new_appmethod = \
                    details(application_dates, stageflag, distribflag, appnumrec, twindow1, stagedays, appmass, appdate,
                            startjul, appmethod, pct1, pct2, twindow2)
                pesticide_mass_soil = \
                    process(new_appmass, new_appnumrec, new_appmethod, plant_factor, rain, start_count, ndates,
                                      soil_distribution_2cm, covmax, foliar_degradation, washoff_coeff)

                start_val = end_val = 0
                for i, v in enumerate(pesticide_mass_soil):
                    if v > 0.0:
                        if start_val == 0:
                            start_val = i
                        end_val = i
                print(start_val, end_val, sum(pesticide_mass_soil))

def qa_process():

    plant_factor = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.0, 1.0, 0.8, 0.0, 0.5, 1.0, 1.0, 1.0])

    appmass = np.array([1.1, 0.5, 1.8, 0.2, 2.3, 0])

    appnumrec = np.array([6, 7, 8, 13, 14, 0])

    appmethod = np.array([1, 2, 2, 2, 1, 0])

    rain = np.array([0, 0, 0, 0, 1.0, 1.1, 1.65, 0, 0, 0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.6, 0, 0]) # 20

    start_count = 0

    ndates = 20

    soil_distribution_2cm = 0.75

    foliar_degradation = 0.2

    washoff_coeff = 0.1

    covmax = 1.0

    pesticide_mass_soil = process_working(appmass, appnumrec, appmethod, plant_factor, rain, start_count, ndates,
                                  soil_distribution_2cm, covmax, foliar_degradation, washoff_coeff)
    print(1, list(pesticide_mass_soil))

    pesticide_mass_soil = process_new(appmass, appnumrec, appmethod, plant_factor, rain, start_count, ndates,
                                  soil_distribution_2cm, covmax, foliar_degradation, washoff_coeff)
    print(2, list(pesticide_mass_soil))


__author__ = 'Trip Hook'
