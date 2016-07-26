module ReadInputs_module
implicit none

contains
subroutine readInputs(baddata)

   use variables_parameters, ONLY: nhoriz,metpath
   
   use variables_inputs  
                                      
   implicit none                              
   integer status , i   
   real dummy
   character :: dummy_var !variables not used in SAM scenarios, but in input matrix
   character(5)::  metnumber
   logical, intent(out) ::baddata
   integer :: crop_position
  
   baddata = .false.

   !**********SAP Pilot input matrix*********************   
   !read(10,*, iostat = status) mukey, metnumber,slope,dummy_var,USLEK,oc,bd(1),fc(1),wp(1),bd(2),fc(2),wp(2),plant_jd,harvest_jd,  &
         !cintcp,covmax,uslec(1),uslec(2),manning_n,sfac,root_max,uslels,cn_PRZM(1),cn_PRZM(2),rain_distrib,anetd
   !*****************************************************
   
   !Input matrix 4/2016
                               !scenario,mukey,cokey,state,cdl,weatherID,date,leachpot,hydgrp,
   read(10,*, iostat = status) scenid,mukey,dummy_var,dummy_var,crop,metnumber,dummy_var,dummy_var,dummy_var, & 
    !cn_ag,cn_fallow,orgC_5,bd_5,fc_5,wp_5,pH_5,sand_5,clay_5,
   cn_PRZM(1),cn_PRZM(2),oc,bd(1),fc(1),wp(1),dummy_var,dummy_var,dummy_var, &
    !orgC_20,bd_20,fc_20,wp_20,pH_20,sand_20,clay_20,                  
   dummy_var,bd(2),fc(2),wp(2),dummy_var,dummy_var,dummy_var, &  
    !orgC_50,bd_50,fc_50,wp_50,pH_50,sand_50,clay_50,
   dummy_var,dummy_var,dummy_var,dummy_var,dummy_var,dummy_var,dummy_var, &   
    !orgC_100,bd_100,fc_100,wp_100,pH_100,sand_100,clay_100,
   dummy_var,dummy_var,dummy_var,dummy_var,dummy_var,dummy_var,dummy_var, &  
    !kwfact,slope,slp_length,uslels,RZmax,sfac,rainfall,anetd
   USLEK,slope,slp_length,uslels,rootzone_max,sfac,rain_distrib,anetd, &
    !plntbeg,plntend,harvbeg,harvend,emrgbeg,emrgend,blmbeg,blmend,matbeg,matend,    
   plant_beg,plant_end,harvest_beg,harvest_end,emerg_beg,emerg_end,bloom_beg,bloom_end,mat_beg,mat_end, &
    !cintcp,covmax,amxdr,irr_pct,irr_type,deplallw,leachfrac,cropprac,cfact_fal,cfact_cov,ManningsN
   cintcp,covmax,root_max,irr_pct,irr_type,irr_depletion,fleach,crop_prac,uslec(1),uslec(2),manning_n

   
   If (status /= 0) then
       baddata = .TRUE.
       write(*,*) "baddata"
       return
   end if
       
   metfile =  trim(metpath)//trim(metnumber)//'_grid.wea'
   
end subroutine readInputs
!************************************************************************

subroutine CheckInputs(baddata)
    use variables_inputs, ONLY: bd, FC
    
    logical,intent(out) :: baddata
    baddata = .false.
    
    If (  any(FC < 0.0009)  .OR.  any(bd < 0.0009) ) then    
        baddata = .true.
    end if
    
end subroutine CheckInputs



!************************************************************************
subroutine process_inputs
   use utilities_module
   use variables_parameters, ONLY: nhoriz,maxyears,increments_1,increments_2,number_soil_incr,delta_x,maxdepth
   
   
   use variables_inputs,    ONLY: bd,oc,fc,wp,anetd,cintcp,afield,slope,&
                                  covmax, root_max,rootzone_max,irrigRate,sfac,nuslec,uslec,cn_PRZM, &
                                  uslek,uslels,uslep,erosion_day,erosion_month,plant_jd,harvest_jd, &
                                  plant_beg,plant_end,harvest_beg,harvest_end
                         
   use variables_Noninputs, ONLY: soil_water_m, &
                                  field_m, wilt_m, bulk_density,org_carbon, &
                                  daily_max_irrigation, &
                                  emergence, maturity, harvest, &
                                  NumberOfYears,num_records, &
                                  startday, firstyear,lastyear, depth,cn,usle_klscp
   implicit none
   integer :: status,i,j
   integer,dimension(maxyears) :: year
   
   integer :: last
   
   integer, allocatable,dimension(:) ::date_holder
   integer, allocatable,dimension(:) ::values
   
   integer, dimension(nuslec) :: index_day
   integer,dimension(nuslec) :: nuslecindex    
    
   real,dimension(num_records)::usle_c_factor
     
   integer ::holder_size
   integer,dimension(num_records):: erosion_vector  

   year = (/(i, i=firstyear,lastyear)/)

   !Dates used for plant and harvest, to derive emergence and maturity
   !can be modified in future to use plant,harvest,emerg,maturity dates from input matrix
   plant_jd = plant_beg
   harvest_jd = harvest_beg
   
   emergence = plant_jd +7 + jd(year,1,1)-startday+1
   maturity = (plant_jd + harvest_jd)/2 + jd(year,1,1)-startday+1
   harvest = harvest_jd + jd(year,1,1)-startday+1
   
   org_carbon = oc/100.
   bulk_density(1:increments_1) = bd(1)
   soil_water_m(1:increments_1) = fc(1)
   field_m(1:increments_1)      = fc(1)
   wilt_m (1:increments_1)      = wp(1)
   
   bulk_density(increments_1+1 :number_soil_incr) = bd(2)
   soil_water_m(increments_1+1 :number_soil_incr) = fc(2)
   field_m     (increments_1+1 :number_soil_incr) = fc(2)  
   wilt_m      (increments_1+1 :number_soil_incr) = wp(2)
   
  !********si conversions*****************************
  bulk_density = bulk_density*1000.     !now in kg/m3
  anetd = min( anetd/100., rootzone_max/100.) !meters - instead of min of anetd & maxdepth, now using anetd & rootzone_max
  cintcp = cintcp/100.                  !now in meters
  afield = 1.                           !hectare, PLACEHOLDER
  afield = afield*10000                 !now in m^2
  slope = slope/100.                    !fraction from percent
  covmax = covmax/100.
  root_max = min(root_max/100., rootzone_max/100.) !meters - instead of min of root_max & maxdepth, now using root_max & rootzone_max
  wilt_m  = wilt_m*delta_x              !meters in each delta_x
  field_m = field_m*delta_x             !meters in each delta_x
  soil_water_m = soil_water_m*delta_x   !meters
  sfac = sfac/100.                      !meters/degree C
  
  !****************************************************
  !the following generates a cumulative vector of depths for use in several routines

  call cumulative_depth (number_soil_incr, delta_x, depth)

  !now find the nodes that correspond to each depth 
  !determine the depth node corresponding to each application deposition in soil
  !allocate  (depi_node(naps), STAT= status)
  !forall(i=1:naps) depi_node(i) = find_node(number_soil_incr,depth,depi(i))
  !now YOU HAVE 2 VECTORS, A cummulative depth vector (depth) and the corresponding node(depi_node)
  !this returns daily values for manning n and product of K*LS*C*P
  !call input_array_processing(nuslec,gduslec,uslec,cn_PRZM,USLEK,USLELS,USLEP)
 
  erosion_vector= 1                    !if no erosion info, then default is first value 
  holder_size = nuslec*numberofyears   !potentially the number of erosion factors in the entire simulation
 
  allocate(date_holder(holder_size),stat=status)
  allocate(values(holder_size),stat=status)
  date_holder = 0

  forall(i=1:nuslec) nuslecindex(i)= i
  do i=1, numberofyears
	forall(j=1:nuslec)	index_day(j)= jd(firstyear+i-1,erosion_month(j),erosion_day(j)) - startday+1 !these are the days when switches occur
    date_holder(i*nuslec-nuslec+1:i*nuslec) = index_day
    values(i*nuslec-nuslec+1:i*nuslec) = nuslecindex
  end do
 
  forall(i=1:holder_size-1, date_holder(i)<num_records .AND.date_holder(i)>0)
    erosion_vector(date_holder(i): min(num_records,date_holder(i+1)))=values(i)
  end forall

  !get daily value of uslec factors and curve number 
  forall(i=1:num_records) usle_c_factor(i) = uslec(erosion_vector(i))
  forall(i=1:num_records) cn(i) = cn_PRZM(erosion_vector(i))
 ! forall(i=1:num_records) s(i) = max(0.,25.4/cn_PRZM(erosion_vector(i)) - 0.254)  !calculation for units of meters, max is here for cn=100 numerical stuff

  !USLE P factor - Default values for contouring practice, crop_prac = 1 (from PRZM, SWAT) 
  !Slope     uslep
  !1-2       0.6
  !2-7       0.5
  !7-12      0.6
  !12-18     0.8
  !18-24     0.9
  !>24       1.0

  if (slope <= 2.0) then 
      uslep = 0.6
  else if (slope > 2.0 .AND. slope <= 7.0) then 
      uslep = 0.5
  else if (slope > 7.0 .AND. slope <= 12.0) then
      uslep = 0.6
  else if (slope > 12.0 .AND. slope <= 18.0) then
      uslep = 0.8
  else if (slope > 18.0 .AND. slope <= 24.0) then
      uslep = 0.9
  else if (slope > 24.0) then
      uslep = 1.0
  end if
  
  usle_klscp = uslek*uslels*usle_c_factor*uslep
 !mmf Note: usle_c_factor, uslep  not in input scenario matrix explicitly

end subroutine process_inputs


end module ReadInputs_module