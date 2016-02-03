

def pesticide_parameters():
    delta_x = 0.02  # meters
    foliar_degradation = 0.0  # per day
    washoff_coeff = 0.1
    soil_distribution_2cm = 0.75  # REVISED for 1 COMPARTMENT - UNIFORM EXTRACTION
    runoff_effic = 0.266
    return delta_x, foliar_degradation, washoff_coeff, soil_distribution_2cm, runoff_effic
