!  FUNCTION:
!  SuperPRZM5hydro - Reads in Recipe files.  For each recipe, daily runoff and erosion are summed for entire watershed.
!  Daily runoff and erosion outputs for recipe watershed are written.

    program SuperPRZM5hydro

   use Variables                       
   use output
   use ReadScenarioFiles
   
   implicit none
   character(len=100) :: inputFile      !NHD flow inputs
   character(len=100) :: dummy
   character(len=100) :: recipes
      
   real               :: t1, t2
   integer            :: count1, count2, yyyy
   integer            :: crop           !CropID from file name
   integer            :: start, last    !Indices to locate CropID in filename 
   integer            :: io_status, ieof, ierror, status
   integer            :: mukey
   
   integer :: i
   integer :: number
   
   character(len=256) :: message   
   character(len=4) :: hydro_only_text
   logical :: hydro_only
   

   real :: badarea 
   integer:: badscenario
   real :: count_rate
   
   call cpu_time (t1)  
   call SYSTEM_CLOCK(count1)

   !Read Chemical Inputs and get Scenario Files
   call get_command_argument(1, eco_or_dw)  !At command line> SuperPRZM5hydro.exe eco
   
   inputFile = "OH_flows_test.csv"                                                                           !Eco
   !inputFile = 'DWI_Monthly_Flows_Reservoir_Only_metric.csv' !'DWI_Monthly_Flows_Flowing_Only_metric.csv'    !DW   
    
   if (eco_or_dw == "eco") then
     recipePath = EcoRecipePath 
     outpath    = EcoOutPath
     flowpath = EcoOutletPath
   else   !dw
     recipePath = DwRecipePath 
     outpath    = DwOutPath
     flowpath = DwOutletPath
   end if

   open (UNIT=10, FILE=trim(adjustl(flowpath))//trim(adjustl(inputFile)), STATUS = 'old', ACTION='read', IOSTAT=ierror)
   if (ierror /=0) then
       stop 'No input file. Program Terminated'
   end if  

  !List of Bad Scenario Files
   open(UNIT = 44, FILE = trim(adjustl(Outpath))// "BadScenarios.txt")
   
 read(10,*)  !header

 do    
     read(10,*) recipe_name    !Reads inputFile (flows) which includes recipe_names for particular state     
     
     !****Now get Scenario Files
         Total_Runoff2010 = 0.0
         Total_Runoff2011 = 0.0
         Total_Runoff2012 = 0.0
         Total_Runoff2013 = 0.0
         
         Total_Erosion2010 = 0.0
         Total_Erosion2011 = 0.0
         Total_Erosion2012 = 0.0
         Total_Erosion2013 = 0.0
         
         !Total_precip = 0.0       !extra
         !Total_potential_et = 0.0 !extra
         
         badarea = 0.0
         totalarea2010 = 0.0
         totalarea2011 = 0.0
         totalarea2012 = 0.0
         totalarea2013 = 0.0
         
         !Recipes 2010-2013        
         in_recipe2010 = "recipe_"//trim(adjustl(recipe_name))//"_cdl2010.txt"  
         in_recipe2011 = "recipe_"//trim(adjustl(recipe_name))//"_cdl2011.txt"  
         in_recipe2012 = "recipe_"//trim(adjustl(recipe_name))//"_cdl2012.txt"  
         in_recipe2013 = "recipe_"//trim(adjustl(recipe_name))//"_cdl2013.txt"  
                 
         open (UNIT=21, FILE=trim(adjustl(recipePath))//trim(adjustl(in_recipe2010)),IOSTAT=ierror, STATUS = 'OLD')
         open (UNIT=22, FILE=trim(adjustl(recipePath))//trim(adjustl(in_recipe2011)),IOSTAT=ierror, STATUS = 'OLD')       
         open (UNIT=23, FILE=trim(adjustl(recipePath))//trim(adjustl(in_recipe2012)),IOSTAT=ierror, STATUS = 'OLD')
         open (UNIT=24, FILE=trim(adjustl(recipePath))//trim(adjustl(in_recipe2013)),IOSTAT=ierror, STATUS = 'OLD')
         
         if (ierror /= 0)  then
             write(*,*) "Missing Recipe File: ", trim(recipe_name)
             cycle
         end if 
        
         if (is_iostat_end(io_status))  exit   !End of File
         write(44,*) trim(recipe_name)
         
         read(21,*) !read header
         read(22,*) 
         read(23,*) 
         read(24,*) 
         
         do 
                   read (21,*, IOSTAT = ierror) recipefilename2010, area2010   !area needs to be in m2 (confirm with Shelly)
                   if (is_iostat_end(ierror))  exit   !End of File
                     scenariofilename2010 = trim(adjustl(Scenario_Path)) // trim(recipefilename2010)            
                     totalarea2010 = totalarea2010 + area2010  
                     yr_recipe = 2010                          
                     call ReadHydroScenarios (scenariofilename2010, badscenario)
                     if(badscenario ==1) then
                         badarea = badarea +area2010
                         write(44,*) trim(recipefilename2010), badarea, "2010"
                     endif
         end do
        do
                   read (22,*, IOSTAT = ierror) recipefilename2011, area2011   !area in m2
                   if (is_iostat_end(ierror))  exit   !End of File
                     scenariofilename2011 = trim(adjustl(Scenario_Path)) // trim(recipefilename2011)            
                     totalarea2011 = totalarea2011 + area2011
                     yr_recipe = 2011
                     call ReadHydroScenarios (scenariofilename2011, badscenario)
                     if(badscenario ==1) then
                         badarea = badarea +area2011
                         write(44,*) trim(recipefilename2011), badarea, "2011"
                     endif
         end do
        do
                   read (23,*, IOSTAT = ierror) recipefilename2012, area2012   !area in m2
                   if (is_iostat_end(ierror))  exit   !End of File
                     scenariofilename2012 = trim(adjustl(Scenario_Path)) // trim(recipefilename2012)             
                     totalarea2012 = totalarea2012 + area2012
                     yr_recipe = 2012
                     call ReadHydroScenarios (scenariofilename2012, badscenario)
                     if(badscenario ==1) then
                         badarea = badarea +area2012
                         write(44,*) trim(recipefilename2012), badarea, "2012"
                     endif
         end do
        do
                   read (24,*, IOSTAT = ierror) recipefilename2013, area2013   !area in m2
                   if (is_iostat_end(ierror))  exit   !End of File
                     scenariofilename2013 = trim(adjustl(Scenario_Path)) // trim(recipefilename2013)             
                     totalarea2013 = totalarea2013 + area2013
                     yr_recipe = 2013
                     call ReadHydroScenarios (scenariofilename2013, badscenario)
                     if(badscenario ==1) then
                         badarea = badarea +area2013
                         write(44,*) trim(recipefilename2013), badarea, "2013"
                     endif
         end do 
         
         close(21)   !Close Recipes 
         close(22)
         close(23)
         close(24)
 
         call Write_Watershed_hydro_txt  ! Write out txt output
    
         write (44,'(A17, F14.0,A14 F14.0, A12,G10.3)') "Total Bad Area =", badarea, "Total Area =", totalarea2010, "2010 Fraction =", badarea/totalarea2010             
         write (44,'(A17, F14.0,A14 F14.0, A12,G10.3)') "Total Bad Area =", badarea, "Total Area =", totalarea2011, "2011 Fraction =", badarea/totalarea2011
         write (44,'(A17, F14.0,A14 F14.0, A12,G10.3)') "Total Bad Area =", badarea, "Total Area =", totalarea2012, "2012 Fraction =", badarea/totalarea2012             
         write (44,'(A17, F14.0,A14 F14.0, A12,G10.3)') "Total Bad Area =", badarea, "Total Area =", totalarea2013, "2013 Fraction =", badarea/totalarea2013
 end do
   
   call cpu_time (t2)  
   CALL SYSTEM_CLOCK(count2,count_rate)
   write(*,*) 'cpu time= ', t2-t1
   write(*,*) 'clock time= ', (count2-count1)/count_rate

     


    end program SuperPRZM5hydro
