module SuperScenarioModule
implicit none 

contains
    subroutine CreateSuperScenario
    use variables_parameters, ONLY:  increments_1, outpath,delta_x
    use variables_noninputs,ONLY: runoff, erosion_loss, &
                                  num_records,      &
                                  numberOfYears,    &
                                  velocity_all,     &
                                  soil_water_m_all, &
                                  bulk_density,org_carbon,     &
                                  depth,    &
                                  rain,             &
                                  plant_factor
                                   
    use variables_Inputs,    ONLY:  crop, covmax , mukey, scenid, & 
                                    plant_beg,plant_end,harvest_beg,harvest_end, &
                                    emerg_beg,emerg_end,bloom_beg,bloom_end,mat_beg,mat_end
    
    integer :: i,count
    character(len=25) :: textscenid 
    character(len=25) :: textmukey
    character(len=20) :: textcrop   !crop ID
    
    integer :: ierror
        
    write(textscenid,*) scenid
    write(textmukey,*) mukey
    write(textcrop,*) crop
    
open (UNIT=88, FILE = trim( adjustl(outpath))//trim(adjustl(textscenid)) // ".csv",IOSTAT=ierror)  !//"_" //trim(adjustl(textcrop))

write(88) covmax

write(88) num_records

write(88) numberOfYears   

count= 0
do i = 1, num_records
   if (runoff(i) > 0.0) then
        count = count+1
   end if
end do

write(88) count

do i = 1, num_records
   if (runoff(i) > 0.0) then
        write(88) i, runoff(i)
   end if
end do

!Erosion loss
do i = 1, num_records
   if (runoff(i) > 0.0) then
        write(88) i, erosion_loss(i)
   end if
end do

write(88) soil_water_m_all(1:increments_1,1:num_records)

count= 0
do i = 1, num_records
   if (velocity_all (0,i) > 0.0) then
        count = count+1
   end if
end do

write(88) count

do i=1, num_records
    if (velocity_all (0,i) > 0.0) then 
      write(88) i, velocity_all (0:increments_1,i)
    end if
end do


write(88) org_carbon   
write(88) bulk_density(1)
write(88) delta_x(1)

write(88) rain(1:num_records)    

write(88) plant_factor(1:num_records)

!Save Julian crop dates to scenarios for now, 
! could be used later by Calculator to determine app dates
write(88) plant_beg
write(88) plant_end
write(88) harvest_beg
write(88) harvest_end
write(88) emerg_beg
write(88) emerg_end
write(88) bloom_beg
write(88) bloom_end
write(88) mat_beg
write(88) mat_end

close(88)

end subroutine CreateSuperScenario


end module SuperScenarioModule