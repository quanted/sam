Program SuperScenario

!rules:
!All units are in s.i. when delivered to main routine except for Time, which is in Days (deg rates etc).
!Unit conversions must take place in input file processing
!dates are referenced to first day of simulation
!first day equals day 1
!date conversions must take place before passing to main routine
!vectorize where possible and efficient
!modules that contain sharable variable are stored in separate modules
!without procedures
!erosion dates (e.g.from record 9 in PRZM) must be ordered low to high

!procedure modules:
use erosion_module
use metfile_processing_module
use plant_growth_module
use SurfaceHydrology
use SuperScenarioModule
use ReadInputs_module
!use runoff_extraction_module  !not used 

use variables_inputs, ONLY: bd, FC
use variables_parameters, ONLY: inputpath
use variables_noninputs, ONLY: erosion_loss

implicit none
character(len=100) :: inputfile, inputname 
character(len=10) :: howmany   !howmany scenarios are in inputfile
integer :: number
integer  ierror, i
real :: t1, t2
integer :: count1, count2
real :: count_rate

integer :: error_flag
logical :: baddata

   call cpu_time (t1)  
   call SYSTEM_CLOCK(count1)

   call get_command_argument(1,inputname) !give filename on command line  !such as:  OH_scenarios_2013_test.txt, MO_scenarios_final.txt
   call get_command_argument(2,howmany)   !give # scenarios in inputfile on command line  !such as: 165156

   inputfile = trim(inputpath)//trim(inputname) 
   print *, "hi", howmany
   pause
   read(howmany,*) number
    print *, "hi"
  pause
   open(UNIT=10, FILE=inputfile,STATUS='old',ACTION='read')
   if (ierror /=0) then
       stop 'No input file. Program Terminated'
   end if  
   read(10,*)
  print *, "hi2"
  pause
   do i= 1, number
       call readinputs(baddata)
       if (baddata) then
          write(66,*) i, 'general bad data read'
          cycle
       end if
 
       call CheckInputs(baddata)
       
       if (baddata) then 
         write (66,*) i,  "Check FC and BD" 
         cycle
       end if
       
     
     call count_met(baddata)
     if (baddata) then 
         write (66,*) i,  "No met file" 
         cycle
     end if
     
     call read_metfile(baddata)
        if (baddata) then 
         write (66,*) i,  "Met File Read problem" 
         cycle
        end if

     call process_Inputs      
     call rain_and_snow_processing     
     call plant_growth   
     call surface_hydrology
     call erosion             !mmf 9/2015 - Erosion routine
    ! call process_extraction !not using extraction routine, single uniform compartment
     call CreateSuperScenario
   end do
    
      call cpu_time (t2)  
      CALL SYSTEM_CLOCK(count2,count_rate)
      write(*,*) 'cpu time= ', t2-t1
      write(*,*) 'clock time= ', (count2-count1)/count_rate

end program SuperScenario