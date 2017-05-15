module erosion_module

contains

subroutine erosion

use variables_noninputs,  ONLY: num_records,rain_and_melt, runoff, cn,usle_klscp,erosion_loss
 use variables_inputs,  ONLY:rain_distrib, manning_n,slope,afield
implicit none

!********************************************************************************
! values are from Table F1 of TR-55, interpolated values are included to make arrays same size
! columns correspond to Ia/P = 
!              .1,      .15,     .2,      .25,     .3,      .35,     .4,       .45,    .5,
real, dimension(9,3) :: Type1(9,3)= &
      reshape((/2.30550, 2.27060, 2.23537, 2.18219, 2.10624, 2.00303, 1.87733, 1.76312, 1.67889, &
                -.51429, -.50908, -.50387, -.48488, -.45695, -.40769, -.32274, -.15644, -.0693,  &
                -.11750, -.10340, -.08929, -.06589, -.02835, -.01983, -.05754, -.00453, 0.0 /),  &
                (/9,3/))

real, dimension(9,3) :: Type1A(9,3)= &
      reshape((/2.03250, 1.97614, 1.91978, 1.83842, 1.72657, 1.70347, 1.68037,  1.65727, 1.63417, &
                -.31583, -.29899, -.28215, -.25543, -.19826, -.17145, -.14463,  -.11782, -.09100, &
                -.13748, -.10384, -.07020, -.02597, 0.02633, 0.01975, 0.01317,  0.00658, 0.0  /), &
                (/9,3/))                 

real, dimension(9,3) :: Type2(9,3)= &
      reshape((/2.55323, 2.53125, 2.50975, 2.48730,  2.46532, 2.41896, 2.36409, 2.29238, 2.20282, &
                -.61512, -.61698, -.61885, -.62071,  -.62257, -.61594, -.59857, -.57005, -.51599, &
                -.16403, -.15217, -.14030, -.12844,  -.11657, -.08820, -.05621, -.02281, -.01259/), &
                (/9,3/))    

real, dimension(9,3) :: Type3(9,3)= &
      reshape((/2.47317, 2.45395, 2.43473, 2.41550,  2.39628, 2.35477, 2.30726, 2.24876, 2.17772, &
                -.51848, -.51687, -.51525, -.51364,  -.51202, -.49735, -.46541, -.41314, -.36803, &
                -.17083, -.16124, -.15164, -.14205,  -.13245, -.11985, -.11094, -.11508, -.09525/), &
                (/9,3/))
!******************************************************************************

real :: L_sheet  !hydraulic length of sheet flow, max=300ft by TR-55
real :: L_shallow
real :: max_retention

real,dimension(num_records)    :: T_conc
real,dimension(num_records)    :: ia_over_p
real,dimension(num_records)    :: c_zero, c_one, c_two
real,dimension(9,3)            :: raintype
real,dimension(num_records)    :: temporary_variable
integer,dimension(num_records) :: lower_index
real,dimension(num_records)    :: remaining
real,dimension(num_records)    :: qp   !peak discharge mm/hr 

temporary_variable=0.
lower_index = 0
remaining=0.
erosion_loss = 0.

!Time of concentration, TR-55
!sheet flow + shallow concentrated flow
L_sheet = min((100.*sqrt(slope)/manning_n),100.)     !( McCuen and Spiess, 1995), slope = ft/ft   !min(hl*3.28084,100.) !in feet, 1 m = 3.2084 ft
L_shallow = (100.*sqrt(slope)/manning_n) - L_sheet   !in feet

ia_over_p = 0.


where (runoff>0.)
    !Time of conc is in hours by TR-55
        !Time of conc for shallow conc flow: T = L_shallow/(3600.*v)
        !v = average velocity, based on unpaved v = 16.1345(slope)^0.5
    !By Velocity method (units are in ft, inch, hour):
    T_conc = 0.007*(manning_n*L_sheet)**0.8/sqrt(rain_and_melt/0.0254)/(slope)**0.4   &
                +  L_shallow/58084.2/sqrt(slope)      !L_shallow/3600./(16.1345*sqrt(slope))
    
!***********************************************
!********rain and melt *************************
  ! ia_over_p = .2*s/rain_and_melt                    ! Ia = 0.2*s, both s and rain are in meters 
    ia_over_p = .0254*(200./cn-2.)/rain_and_melt      ! 0.2*s, in inches
end where

select case (rain_distrib)
    case (1)
        raintype = type1
    case (2)
        raintype = type1A
    case (3)
        raintype = type2
    case (4)
        raintype = type3
end select

!lower limit according to TR-55
where(runoff>0 .AND. ia_over_p <= 0.1)
    C_zero = raintype(1,1)
    C_one  = raintype(1,2)
    C_two  = raintype(1,3)
end where

!upper limit according to TR-55
where(runoff>0 .AND. ia_over_p >= 0.5)
    C_zero = raintype(9,1)
    C_one  = raintype(9,2)
    C_two  = raintype(9,3)
end where

!interpolation of intermediate values
where (runoff>0. .AND. ia_over_p < 0.5  .AND. ia_over_p > 0.1)
    temporary_variable = 20.*(ia_over_p - 0.05)
    lower_index = int(temporary_variable)
    remaining = mod(temporary_variable,1.)
    C_zero = raintype(lower_index,1) +remaining*(raintype(lower_index+1,1)-raintype(lower_index,1))
    C_one  = raintype(lower_index,2) +remaining*(raintype(lower_index+1,2)-raintype(lower_index,2))
    C_two  = raintype(lower_index,3) +remaining*(raintype(lower_index+1,3)-raintype(lower_index,3))
end where

temporary_variable =0.
!temporary variable here is the unit peak discharge
where (runoff >0.) temporary_variable =10.**(C_zero+C_one*log10(T_conc)+c_two*(log10(T_conc))**2)

! 1 sqmi = 258.998811 ha
! 1 sq meter = 3.86102e-7 sq mile
! 1 sq mi = 
! 1 meter = 39.37 inch
! 1 sq meter = 10.76391 sq ft
! 1 foot = 304.8 mm

!peak_discharge = temp_variable*(afield/2589988.11 sqmi)*(runoff*39.370079) /(Afield*10.7639104 ft2/m2)*(3600 sec/hr)*(304.8 mm/ft)
!               = temp_variable*runoff*3600.*304.8*39.370079/2589988.11/10.7639104

! qp(mm/hr) = 1.549587*runoff in meters* qu
!Fp assumed to be 1 (no swamp or pond on field)

where (runoff >0.) qp = 1.54958679*runoff*temporary_variable  !peak discharge (mm/hr)


!Erosion loss (kg) - MUSS equation (only for "small" watersheds) - not used here
!!where (runoff >0.) erosion_loss = 0.79*(runoff*1000.*qp)**.65*(afield/10000.)**.009*usle_klscp*1000.  !kg/d

!Erosion loss (kg) - MUSLE equation (Williams, 1975)  erosion_loss = 1.586*(runoff*1000.*qp)**.56*afield/10000.)**.12
where (runoff >0.) erosion_loss = 1.586*(runoff*1000.*qp)**.56*usle_klscp*1000.    !kg/d   !Area term of eqn moved to SuperPRZM5Hydro: *(afield/10000.)**.12
where (isnan(erosion_loss)) erosion_loss = 0.0   !remove NaN values

end subroutine erosion


end module erosion_module