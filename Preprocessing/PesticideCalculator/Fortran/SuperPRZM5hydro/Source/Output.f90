module output
implicit none

contains
   
subroutine Write_Watershed_hydro_txt
      use Variables 
      implicit none
      character(len= 100) :: filename
      integer :: j
      integer :: start_numrec, end_numrec
      
      if (eco_or_dw == "eco") then
          filename = trim(adjustl(recipe_name)) // "_hydro.txt"   !For EcoOutput txt files
      else
          filename = trim(adjustl(recipe_name)) // "_hydro.txt"   !For dwOutput txt files
      end if
      
      open(UNIT= 57, recl = 800, FILE = trim(Outpath) // filename)
    
    !Add 2014, 2015 when available, e.g. totalarea2014, totalarea2015
    write(57,*) "num_records,Area2010_m2,Area2011_m2,Area2012_m2,Area2013_m2"  
    write(57,*) num_records, totalarea2010, totalarea2011, totalarea2012, totalarea2013   !in m2
    write(57,*) "R2010_m3,R2011_m3,R2012_m3,R2013_m3,E2010_kg,E2011_kg,E2012_kg,E2013_kg"
    !Runoff in m3, Erosion in kg
    do j=1, num_records  
     !Add 2014, 2015 when available 
     write(57,*) Total_Runoff2010(j), Total_Runoff2011(j), Total_Runoff2012(j), Total_Runoff2013(j), Total_Erosion2010(j), Total_Erosion2011(j), Total_Erosion2012(j), Total_Erosion2013(j)
    end do 
    
     close(57)
       
   end subroutine Write_Watershed_hydro_txt   
   
end module output