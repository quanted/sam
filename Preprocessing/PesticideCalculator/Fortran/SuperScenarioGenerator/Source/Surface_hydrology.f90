module SurfaceHydrology
!stores runoff, velocity_all, soil_water_m_all in module

contains
    !************************************************************
    !************************************************************
    subroutine surface_hydrology

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    use variables_NonInputs,  ONLY: num_records,      & !input
                                  rain,             &   !input
                                  rain_and_melt,    &   !input
                                  potential_et,     &   !input
                                  soil_water_m,daily_max_irrigation,&                                
                                  depth,            &
                                  field_m,          &
                                  wilt_m,           &
                                  plant_factor,    &
                                  cn ,           &
                                  runoff,           &   !output
                                  velocity_all,     &   !output
                                  soil_water_m_all, &   !output
                                  velocity,         &
                                  available_water_m      
                                  
    use variables_Inputs, ONLY:    anetd, cintcp, root_max, rootzone_max, irr_type,  &
                                   irrigation_depth,     &
                                   irr_depletion,      &
                                   fleach, irrigRate 

     use variables_parameters, ONLY: number_soil_incr,delta_x, cn_moisture_depth  

    use utilities_module
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    implicit none

    integer :: irrigation_node
    real    :: target_dryness
    real    :: total_fc           !total field capacity = amount of water after irrigation (m)
    real    :: leaching_factor
    real    :: available_soil_et
    real    :: overcanopy_irrigation
    real    :: undercanopy_irrigation
    real    :: delta_water 
    real    :: canopy_holdup
    real    :: available_canopy_gain
    real    :: leaching
    real    :: irrig_required,current_dryness

    real    :: effective_rain    !irrigation + snowmelt + rain
    integer :: day
    integer :: node
    real    :: s
    
    real    :: cn_1, cn_3
    
    real    :: et_from_canopy
    
    integer :: cn_moisture_node
    real    :: wp_plus_fc
    real    :: wp_plus_fc_over2
    real    :: antecedent_moisture
    
    real    :: curve_number                       !local CN, note: CN soil moisture adjustment removed below - mmf 9/2015

    real,dimension(number_soil_incr) :: fc_minus_wp

    integer :: i
    integer :: evapo_node
    integer,dimension(num_records):: et_node
    real,dimension(num_records)   :: et_depth

    real,dimension(number_soil_incr) :: et_factor !reduction of et with depth, could be dimensioned by max(root_max, anetd)
    real,dimension(number_soil_incr) :: soil_layer_loss
    real :: check_moisture_et 
    real :: target_moisture_et
    real :: water_level

    real :: xx

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    available_canopy_gain =0.
    canopy_holdup = 0.
    runoff= 0.           !by initialization, runoff is zero when rain < 0.2S
    velocity_all = 0.
    soil_water_m_all = 0.

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   !%% Some initial calculations that dont need to be included in loop %%%%%%%%%%%%%%%%
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    fc_minus_wp = field_m - wilt_m
      
    if (irr_type > 0) then
       irrigation_depth = root_max     !meters 
       irrigation_node = find_node(number_soil_incr,depth,irrigation_depth)
       target_dryness = sum( fc_minus_wp(1:irrigation_node)*irr_depletion + wilt_m(1:irrigation_node)  )
       total_fc = sum(field_m(1:irrigation_node))
       leaching_factor = fleach +1.
    end if

    where (plant_factor > 0.)
        et_depth = max(plant_factor*root_max,anetd)
    elsewhere
        et_depth = anetd
    endwhere

    evapo_node =  find_node(number_soil_incr,depth,anetd)  !node only for evaporation
    et_node = evapo_node                                   !initially set all to the minimum
   
    forall(i=1:num_records, et_depth(i) > anetd)
       et_node(i) = find_node(number_soil_incr,depth,et_depth(i))
    end forall
    
    cn_moisture_node = find_node(number_soil_incr,depth,cn_moisture_depth)        !removed CN soil moisture adjustment below

    wp_plus_fc = sum(field_m(1:cn_moisture_node))+sum(wilt_m(1:cn_moisture_node)) 
    wp_plus_fc_over2 = wp_plus_fc/2.                                              

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    do day = 1,num_records
      effective_rain = rain_and_melt(day)
      !Calculated Max Irrig Rate - causing no runoff (based on WQTT Advisory Irrigation Guidance) mmf 3/2016
      irrigRate = 0.2*((2540./cn(day)) - 25.4)/100.  !cm/day -> m/day  max rate of irrig water that will not cause runoff
      
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%%%%   check soil moisture for irrigation  %%%%%%%%%%%%%%%%%%
      !%%  irrigation only if flagged, its dry, and its not raining %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      overcanopy_irrigation=0.
      undercanopy_irrigation=0.

      if (irr_type > 0) then
          current_dryness = sum(soil_water_m(1:irrigation_node))
          if (current_dryness<target_dryness.AND. rain_and_melt(day) <= 0.) then !irrigation today
              
              irrig_required = (total_fc -current_dryness)*leaching_factor       !water to be added to bring irrigation zone to field capacity
              daily_max_irrigation = irrigRate    !m/day 
              
              select case (irr_type)
                   case (3)
                        !Maximum amount of irrigation is limited by "daily_max_irrigation"
                        overcanopy_irrigation = min(irrig_required,daily_max_irrigation)
                        effective_rain = overcanopy_irrigation
                   case (4)
                        undercanopy_irrigation = min(irrig_required,daily_max_irrigation)
                        effective_rain = undercanopy_irrigation
                   end select
          end if
      end if
       
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%%  Soil Moisture Curve Number Option   %%%%%%%%%%%%%%%%%    !removed CN soil moisture adjustment below, mmf 3/2015

      if(effective_rain >0.) then
              cn_1 = cn(day)/(2.281-0.01281*cn(day))  !no longer used, mmf 9/2015
              cn_3 = cn(day)/(0.427+0.00573*cn(day))  !no longer used, mmf 9/2015

              antecedent_moisture = sum(soil_water_m(1:cn_moisture_node))
              !if (antecedent_moisture <=wp_plus_fc_over2) then         !Removed CN soil moisture adjustment, mmf 9/2015
                 !curve_number = antecedent_moisture/wp_plus_fc_over2*(cn(day)-cn_1)+cn_1
              !else
                 !curve_number = (antecedent_moisture-wp_plus_fc_over2)/wp_plus_fc_over2*(cn_3-cn(day))+cn(day)
              !endif
             !s= 25.4/curve_number -.254
              
              s= 25.4/cn(day) -.254   !Revised, CN soil moisture adjustment removed, mmf 9/2015
              
      end if
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !###### Runoff by the Curve Number Method #################
        if(effective_rain > 0.2*s)then
            runoff(day) = (effective_rain - 0.2*s)**2/(effective_rain + 0.8*s)
        end if
      !##########################################################


      !%%%%%%%%%%   Canopy Calculations   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !the maximum available water that is available for canopy holdup of precip
      !assume that no snow accumulates on canopy; probably not important for crops
      !by curve number method, it is a given that "runoff" will occur and not be
      !a function of canopy hold up. Therefore, only that amount of precip that did not 
      !run off will be available for canopy hold up. In other words, canopy holdup is not a
      !given, but rather a maximum potential.  PRZM on the other hand modifies Ia by canopy hold up
      !which is inconsistent with the fundamentals of curve number method

      !****** Calculate Leaching into top layer ***************************************
      !Leaching is the amount of rain plus snowmelt, minus canopy holdup (prior to evap), & minus runoff

      if (rain(day) > 0. .OR. overcanopy_irrigation >0. )then !backing out the runoff contributions of rain and melt
           available_canopy_gain = (rain(day)+overcanopy_irrigation) *(1.- runoff(day)/effective_rain)
           delta_water = min(available_canopy_gain, cintcp*plant_factor(day)- canopy_holdup)
           canopy_holdup = canopy_holdup + delta_water
           leaching = effective_rain - runoff(day) - delta_water
      else
           leaching = effective_rain - runoff(day)  !this includes snowmelt,no need for canopy holdup update
      end if
  
      !%% update canopy holdup here
      et_from_canopy = canopy_holdup - potential_et(day)
      canopy_holdup     = max(0., et_from_canopy)
      available_soil_et = max(0.,-et_from_canopy)


!############## ET Calculations #####################################
      
  available_water_m(1:et_node(day)) = soil_water_m(1:et_node(day))-wilt_m(1:et_node(day))

      !Reduction factor below 0.6 available water
      check_moisture_et  = sum(soil_water_m(1:et_node(day))-wilt_m(1:et_node(day)))
      target_moisture_et = 0.6*sum(fc_minus_wp(1:et_node(day)))
      if (check_moisture_et <target_moisture_et ) then
          available_soil_et = available_soil_et*check_moisture_et/target_moisture_et
      end if

      !Reductions by depth and available water
      et_factor(1:et_node(day))= (depth(et_node(day))-depth(1:et_node(day)) +delta_x(1:et_node(day)) )*available_water_m(1:et_node(day))

      xx= sum(et_factor(1:et_node(day)))  !if completely dry, ET should be zero, prevents division by zero

      if (xx > 0.) then
        et_factor(1:et_node(day))= et_factor(1:et_node(day))/xx
      else
        et_factor(1:et_node(day))=0.
      endif

      soil_layer_loss =0.
      soil_layer_loss(1:et_node(day))=available_soil_et*et_factor(1:et_node(day))  !potential loss

!############################################################################
!********** Update soil water content due to leaching ***********************
!********** Evaporate the soil water ****************************************

      !Prepare to calculate the daily velocities throughout profile
      node=0
      velocity(0) =leaching

      do    !velocity loop
            !the "node" will be the node that loop ends on and will be partially filled.
            !if ends on last node, then this node could be filled entirely

            node=node+1  !update node

            water_level = velocity(node-1) -soil_layer_loss(node)+ soil_water_m(node)

            if(water_level > field_m(node)) then
                    velocity(node) = water_level-field_m(node)
                    soil_water_m(node) = field_m(node)         
            else
               velocity(node) = 0.
               soil_water_m(node) = max(water_level,wilt_m(node))             
            end if
              
            if (velocity(node) <=0. .AND. node> et_node(day) )then  !leave loop early because remaining calcs are unnecessary
                velocity(node:number_soil_incr) = 0.                !and soil water does not change
            exit        !velocity loop
            else
            end if
            
            if (node == number_soil_incr) exit !velocity loop if hit the bottom of soil
      end do            !velocity loop
           
        !Here's the output for this loop
        velocity_all(:,day) = velocity 
        soil_water_m_all(:,day) = soil_water_m
    end do

    end subroutine surface_hydrology

end module SurfaceHydrology