module variables
implicit none
save
    
    !***** Recipe Path ***************************
    character(len=100)  :: RecipePath 
    character(len=100)  :: DwRecipepath =   '..\..\..\dwRecipes\dwi2012recipes_second_batch_100114\' 
    character(len=100)  :: EcoRecipepath =  '..\..\..\EcoRecipes_nhd\OH_test_recipes_mar2016\'   
    
    character(len=100)  :: FlowPath
    character(len=100)  :: DwOutletpath =   '..\..\..\dwRecipes\dwi2012recipes_second_batch_100114\'
    character(len=100)  :: EcoOutletpath =  '..\..\..\OutletFlows_nhd\states\'  
    
    !***** binScenario Path ***************************
    character(len=100)  :: Scenario_Path =  '..\..\..\binScenarios_all\binScenarios_withErosion\'  !..\..\..\binScenarios_all
     
    !**** OUTPUT PATH **************************************
    character(len=100)  :: Outpath
    character(len=100)  :: EcoOutpath =  '..\..\..\EcoOutput_all\EcoOutput_test_OH\' 
    character(len=100)  :: DwOutpath =  '..\..\..\dwOutput_all\dwOutput_SG_CPC_rev\'
    character(len=4)    :: eco_or_dw       !specifies eco or dw path
    
    character(len=100)  :: recipe_name           !current recipe being processed
    character(len=100)  :: inpath = '..\..\..\'  
    character(len=100)  :: infile = '\scenarios\scenario_mar2016\OH_scenarios_2013_test.csv'
    
    character(len=100) :: in_recipe2010, in_recipe2011, in_recipe2012, in_recipe2013
    character(len=100) :: recipefilename2010, recipefilename2011, recipefilename2012, recipefilename2013
    character(len=100) :: scenariofilename2010, scenariofilename2011, scenariofilename2012, scenariofilename2013

    character(len=100)  :: metfile
    character(5)        :: metnumber
    character(len=100)  :: mukey_crop
    
    integer :: yr_recipe
    
    integer ::  eof   !end of file flag
    
    !Parameters
    integer, parameter :: maxdays = 24000   !# days used for allocations
    
    ! OUTPUT******************************************
    real, dimension(maxdays):: Total_Runoff2010, Total_Runoff2011, Total_Runoff2012, Total_Runoff2013        !m3 runoff from watershed
    real, dimension(maxdays):: Total_Erosion2010, Total_Erosion2011, Total_Erosion2012, Total_Erosion2013    !kg erosion loss from watershed
    
    real :: totalarea2010,totalarea2011,totalarea2012,totalarea2013            !m2, watershed area

    integer, parameter :: maxsoilinc = 50

   integer      :: num_records    
   integer      :: NumberOfYears
   real      :: pfac

   real,dimension(maxdays)                        :: rain              !only precip above zero C
   
   !met data - used previously for tests
   !real,dimension(maxdays)                        :: temp
   !real,dimension(maxdays)                        :: precip
   !real                                           :: evap
   !real,dimension(maxdays)                        :: potential_et
   !real,dimension(maxdays)                        :: Total_precip
   !real,dimension(maxdays)                        :: Total_potential_et
  
   real,dimension(maxdays)                        :: plant_factor
   real,dimension(maxsoilinc)                     :: depth             !vector of cumulative depth (m)

   real,dimension(maxdays)                        :: runoff   !in m
   real,dimension(maxdays)                        :: erosion  !in kg/ha
   integer,dimension (maxdays)                    :: day      !index of days when runoff occurs
      
    real,dimension(maxdays) :: leached_mass, stored_mass, runoff_mass, degraded_mass   
    real :: totalMassRunoff, totalApplied
    
    real               :: area2010,area2011,area2012,area2013  !area in m2
    
           
end module variables