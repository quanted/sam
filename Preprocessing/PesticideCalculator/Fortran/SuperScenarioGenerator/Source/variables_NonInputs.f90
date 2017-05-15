module variables_NonInputs
    use variables_parameters,ONLY: maxdays, maxyears, number_soil_incr

   implicit none
   save

   real,dimension(maxdays) :: temp          !temp as read from met file (C)
   real,dimension(maxdays) :: precip        !output: daily precipitation (m)
   real,dimension(maxdays) :: potential_et
   real,dimension(maxdays) :: rain          !only precip above zero C
   real,dimension(maxdays) :: rain_and_melt !effective "rain" for runoff and erosion
   real,dimension(maxdays) :: plant_factor
   real,dimension(maxdays) :: usle_klscp    !daily product of L*LS*C*P  
   real,dimension(maxdays) :: runoff        !m, 
   real,dimension(maxdays) :: erosion_loss  !kg/ha daily loss
   real,dimension(maxdays) :: snow_melt
   real,dimension(maxdays) :: cn  
   
  !curve number dimension: num_records
   real,dimension(number_soil_incr):: soil_water_m    !array of soil water in profile
   real,dimension(number_soil_incr):: bulk_density    !stored as kg/m3, read in as g/ml
   real,dimension(number_soil_incr):: depth           !vector of cumulative depth (m)
   real,dimension(number_soil_incr):: wilt_m          !meters in each segment
   real,dimension(number_soil_incr):: field_m         !meters in each segment

   real :: org_carbon  !fraction oc
   
   real    :: daily_max_irrigation ! m/day

   integer,dimension(maxyears):: emergence 
   integer,dimension(maxyears):: maturity
   integer,dimension(maxyears):: harvest

   !variables obtained from the met file
   integer :: NumberOfYears
   integer :: num_records

  !Time Variables as dictated by metfile
   integer :: startday
   integer :: firstyear
   integer :: lastyear

   !calculated and transferrred variables:
     
   !You might want to make these allocatable to prevent heap issue and overflow, for now Heap Array set to 0 in properties
   real,dimension(0:number_soil_incr,maxdays)  :: velocity_all       !daily vertical water velocity
   real,dimension(number_soil_incr,maxdays)    :: soil_water_m_all  

   real,dimension(0:number_soil_incr) :: velocity
   real,dimension(number_soil_incr)   :: available_water_m


end module variables_NonInputs