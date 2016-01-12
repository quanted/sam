1/12/2016 - James (Trip) Hook, OPP/EFED

This directory contains the working Python version of SAM (Spatial Aquatic Model).

The two primary components of the model are the PESTICIDE CALCULATOR, which processes input scenarios and watershed
recipes to assess pesticide applications and transport of pesticides into waterways, and the TIME OF TRAVEL engine,
which take output from the pesticide calculator and perform a downstream convolution to simulate dispersion of pesticide
in water as the pulse moves downstream.

The file "pesticide_calculator.py" is the main script for running the pesticide calculator, and "time_of_travel.py" is
the script that runs the time of travel engine.  All other scripts in the \Tool directory are libraries used by these
two executable scripts.

The folder \Preprocessing contains scripts that are used to generate input files for the pesticide calculator and
time of travel engine.  At the moment, these tools are primarily for time of travel inputs, but may be expanded to
include other preprocessing routines as time goes on.

The folder \QA_Tools contains scripts used to assess the results of the primary SAM tools.  This library is still
weakly developed.