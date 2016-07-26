module variables_inputs
    use variables_parameters,ONLY: nuslec, nhoriz

    implicit none
    save
    
    character(len=25)    :: scenid  
    integer              :: mukey  
    integer              :: crop
    
    character(len=120)   :: metfile  !name of metfile  
    real                 :: pfac
    real                 :: sfac     !snowmelt factor
    real                 :: anetd
    real                 :: USLEK
    real                 :: USLELS
    real                 :: USLEP
    real                 :: AFIELD
    integer              :: rain_distrib
    real                 :: slope
    real                 :: slp_length
    real                 :: cintcp
    real                 :: root_max       !max rooting depth for crop from PRZM
    real                 :: rootzone_max   !max root zone depth from SSURGO value-added parameters
    real                 :: covmax
                      
    integer,dimension(nuslec)  :: erosion_day 
    integer,dimension(nuslec)  :: erosion_month    
    real,dimension(nuslec)     :: USLEC, USLEC2 
    integer,dimension(nuslec)     :: cn_PRZM
       
    real    :: manning_n 
    
    !Crop stage inputs
    integer :: plant_jd, plant_beg, plant_end
    integer :: harvest_jd, harvest_beg, harvest_end
    integer :: emerg_beg, emerg_end
    integer :: bloom_beg, bloom_end
    integer :: mat_beg, mat_end
        
   !IRRIGATION INPUTS - updated mmf 9/2015
    real       :: irr_pct                 !percent of crop area irrigated - to be collected by Kurt
    integer    :: irr_type                !irrigation type 3 or 4 for over/undercanopy irrigation 
    real       :: irr_depletion           !allowed depletion (default) 
    real       :: fleach                  !extra water fraction (default) 
    character(len=5) :: crop_prac         !dominant crop practice used for deriving usle p factor - to be collected by Kurt
    
    real       :: irrigRate               !in surface_hydrology.f90, 0.2*((2540./cn)-25.4)/100.  !max rate causing no runoff, m/day based on WQTT Advisory Irrigation Guidance
    real       :: irrigation_depth        !default irrigation to root depth as in SWCC, mmf 9/2015
        
 !  real,dimension(2)    :: thkns
    real,dimension(2)    :: bd
    real,dimension(2)    :: FC
    real,dimension(2)    :: WP
    real                 :: OC
                                   
                    
end module variables_inputs
