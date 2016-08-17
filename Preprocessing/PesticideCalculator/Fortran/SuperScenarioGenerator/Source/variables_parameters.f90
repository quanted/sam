module variables_parameters
implicit none

character(len=100),parameter :: inputpath = '..\..\..\scenarios\scenario_mar2016\'
character(len=100),parameter :: metpath = '..\..\..\Met_nc\Met_txt_linint2\1961-2014\'
character(len=100),parameter :: outpath = '..\..\..\binScenarios_all\binScenarios_withErosion\'

!real,parameter :: extraction_depth  = 0.02   !m, depth = 2 cm changed to just defining node below
integer, parameter :: maxyears= 65
integer, parameter :: maxdays = 19723         !for current met; before it was 23742 for 1/1/1948 to present (this is the # of days used for allocations, equal to # of days in met files)

real,parameter :: washoff_depth     = 0.02    !m, depth at which foliar washoff deposits

real,parameter :: cn_moisture_depth = 0.1     !m, depth at which moisture is checked for CN adjustment, BUT no longer using adjustment, mmf 3/2015

integer, parameter :: nuslec = 2              !number of cover management practices

integer, parameter :: nHoriz=2

real, parameter    :: thkn_layer1   = 2.0     !cm, thickness of first layer
real,parameter     :: thkn_layer2   = 100.0   !cm, thickness of 2nd layer
real,parameter     :: maxdepth = 1.02         !total depth in meters

integer, parameter :: increments_1  = 1       !number of increments in top 2-cm layer: 1 COMPARTMENT, UNIFORM EXTRACTION

integer, parameter :: increments_2  = 20      !number of increments in 2nd 100-cm layer (not used in extraction)

integer,parameter  :: number_soil_incr = 21   !total number of soil compartments: 1 top 2-cm compartment + 20 compartments below

!delta_x is in meters
real, dimension(number_soil_incr),parameter :: delta_x  = (/0.02,     &
                                              0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,               &
                                              0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05/)


end module variables_parameters