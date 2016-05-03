module metfile_processing_module

!this module reads the PRZM-type metfile
!it populates the values of variables from the metfile
!as stored in the "variables_weather" module,
    contains
    
    
subroutine count_met(baddata)
    !************************************************************************
    !read met file and count records for array allocation, returns num_records
    !stores num_records in weather_parameters_module
    !************************************************************************
     use variables_noninputs, ONLY: num_records, &
                                    NumberOfYears, & !output
                                    startday, &      !output
                                    firstyear,  &    !output
                                    lastyear         !output
                               
     use variables_Inputs,  ONLY: metfile            !input
     use utilities_module	
         
     implicit none
     integer :: ieof,ierror  !signal for no more records '/=0' = no record 
     integer :: firstmon
     integer :: firstday
     integer :: lastmon
     integer :: lastday 
     integer :: year
     logical,intent(out) :: baddata
  
     
     baddata = .false.
     open (UNIT=5, FILE=metfile,STATUS='old',ACTION='read',IOSTAT=ierror)

     if (ierror /=0) then
       baddata = .True.
       return
     else                      !file exists, count records 
        num_records=0
        read (5,*, IOSTAT=IEOF) firstmon,firstday, firstyear
        rewind (5)
        do      
           read (5,*, IOSTAT=IEOF)
           if (ieof /= 0)exit  !at end of file
           num_records=num_records+1
        end do
     end if
     close(UNIT=5)
         
     startday = jd (firstyear,firstmon,firstday)
       
     call get_date (startday+num_records-1, lastYEAR,lastMON,lastDAY)
     NumberOfYears = lastYear-firstYear + 1
 
end subroutine count_met


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_metfile(baddata)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !This subroutine reads the metfile and puts precip,wind evap,temp in arrays
  ! as well as startday into weather_parameters_module
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   use variables_inputs, ONLY: metfile,pfac            !input

                               
   use variables_noninputs,ONLY: num_records,&
                                 precip,temp,  &
                                 potential_et, &  !output
                                 startday         !output
  !**************************************************************************
   implicit none

   integer ::  eof             !end of file flag
   integer :: i                !do loop counter
   integer :: ierror
   real    :: evap             !output: daily evaporation (m)
   real :: dummy
   logical,intent(out) :: baddata
   
   baddata = .false.
   eof=0 
   ierror=0
   precip = 0.
   potential_et= 0.
   
   !%%%%%%%%%%%%%%  Check if Met file Exists  %%%%%%%%%%%%%%%%%%%%
   open (UNIT=5, FILE=metfile,STATUS='old',ACTION='read',IOSTAT=ierror)
   !order of met variables- mm,dd,yyyy,precip,ET,temp,wind,sr 
   !Read met file
      do i=1,num_records 
          read (5,*, IOSTAT=eof) dummy, dummy, dummy, precip(i),potential_et(i),temp(i), dummy, dummy   
         !potential_et(i) = evap*pfac  !no longer need conversion, ET is available directly in met file
      
      !%%%%%%% Get rid of any missing precip values %%%%%%%%%%%%%%
      if (precip(i) > 30000.) then
        precip(i) = 0.0
      end if

         if (eof /=0) then               !met file read error check
          baddata = .true.
          return
         end if
      end do

       !%%%% si conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       precip = precip/100.                  !convert to meters
       potential_et = potential_et/100.      !convert to meters 
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

       close (UNIT=5)
end subroutine read_metfile
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine rain_and_snow_processing
   !this subroutine calculates rain_and_melt and snow_melt 
  
   
    use variables_inputs,    ONLY: sfac            !input
    use variables_noninputs, ONLY: num_records,          & 
                                   temp,          & !input
                                   precip,        & !input
                                   rain,          & !input
                                   rain_and_melt, & !output
                                   snow_melt        !output      
   !*********************************************************************
   implicit none

   real,dimension(num_records):: snow  !snow days
   integer :: i
   real:: snow_accumulation			   !running tracker of snow accumulation

   !*** Find rain and snow days ***********************************************
   rain = precip               !initially set rain to precip
   
   where (temp <= 0.)rain = 0. !this takes snow out of the precip, "rain" is then only rain

   snow = 0.
   snow_accumulation =0.
   snow_melt=0.
   where (temp <= 0.)  snow = precip  !only the snow fall

   !****************************************
   do i=1, num_records
     if (temp(i) <= 0) then
         snow_accumulation = snow_accumulation +snow(i)
     else if (snow_accumulation > 0.0)then    !if here to speed calcs during warm periods
         snow_melt(i) = min(sfac*temp(i),snow_accumulation)  !snow melt for runoff calcs
         snow_accumulation = snow_accumulation - snow_melt(i)
     end if
   end do
   rain_and_melt = rain + snow_melt   !add rain and snow melt for use with curve number

end subroutine rain_and_snow_processing


end module metfile_processing_module