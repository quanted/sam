Module ReadScenarioFiles
implicit none

contains
    

subroutine ReadHydroScenarios(filename,badscenario)
   !Reads scenarios and accumulates daily runoff and erosion
    use Variables 
    
    implicit none
    integer,intent(out) :: badscenario
    character(len=100),intent(in) :: filename
    integer :: ierror, reason, count, i, j, k, n, status
    logical :: baddata
    integer :: number 
    real :: dummy2
    
    baddata = .false.
    
    badscenario = 0
 
    open (UNIT=88, FILE=filename,IOSTAT=ierror, STATUS = 'OLD', FORM ='UNFORMATTED')   
    
    if (ierror /= 0) then 
      badscenario = 1
      runoff = 0.0
      erosion = 0.0
      return
    end if
        
    read(88,IOSTAT=reason) dummy2
    if (reason /= 0) then
        badscenario = 1
        return
    end if 
    
    read(88,IOSTAT=reason) num_records
    if (reason /= 0) then
        badscenario = 1
        return
    end if
    
    read(88,IOSTAT=reason) numberOfYears
    if (reason /= 0) then
        badscenario = 1
        return
    end if    
    
    read(88,IOSTAT=reason) count
    if (reason /= 0) then
        badscenario = 1
        return
    end if
    
    !runoff
    do i = 1, count
        read(88,IOSTAT=reason) day(i), runoff(i)     !m
        if (reason /= 0) then
            badscenario = 1
            return
        end if
    end do 
    
    !erosion
    do i = 1, count
        read(88,IOSTAT=reason) day(i), erosion(i)     !kg/d  (MUSLE area term brought in below: *(area of field/10000.)**.12)
        if (reason /= 0) then
            badscenario = 1
            return
        end if
    end do 
    
    close(88)
        
    if (yr_recipe == 2010) then
        Total_Runoff2010(day(1:count)) = Total_Runoff2010(day(1:count)) + runoff(1:count)*area2010                      !m3, daily
        Total_Erosion2010(day(1:count)) = Total_Erosion2010(day(1:count)) + erosion(1:count)*((area2010/10000.)**.12)   !10000. m2 / ha
    end if    
    if (yr_recipe == 2011) then
        Total_Runoff2011(day(1:count)) = Total_Runoff2011(day(1:count)) + runoff(1:count)*area2011                      !m3, daily
        Total_Erosion2011(day(1:count)) = Total_Erosion2011(day(1:count)) + erosion(1:count)*((area2011/10000.)**.12)   !kg, daily
    end if    
    if (yr_recipe == 2012) then
        Total_Runoff2012(day(1:count)) = Total_Runoff2012(day(1:count)) + runoff(1:count)*area2012
        Total_Erosion2012(day(1:count)) = Total_Erosion2012(day(1:count)) + erosion(1:count)*((area2012/10000.)**.12)
    end if    
    if (yr_recipe == 2013) then
        Total_Runoff2013(day(1:count)) = Total_Runoff2013(day(1:count)) + runoff(1:count)*area2013
        Total_Erosion2013(day(1:count)) = Total_Erosion2013(day(1:count)) + erosion(1:count)*((area2013/10000.)**.12)
    end if   
    
    !Total_Runoff(day(1:count)) = Total_Runoff(day(1:count)) + runoff(1:count)*area      !PREVIOUS, before multiple recipe years 
    
  
end subroutine ReadHydroScenarios
    
end Module ReadScenarioFiles