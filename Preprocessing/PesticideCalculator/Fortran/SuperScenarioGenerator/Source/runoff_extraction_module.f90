module runoff_extraction_module
    integer,private::status

!Not currently used
    
contains
    !
    !subroutine process_extraction
    !
    !!this routine returns exrtraction array, which is the effective flow rate (meters) through each of the 
    !!extraction compartments
    !!Extraction relation is from PRZM manual
    !use variables_parameters,  ONLY: extraction_index
    !use variables_nonInputs,   ONLY: num_records,delta_x, number_soil_incr,depth,runoff, extraction_array 
    !implicit none
    !
    !!****Local Variables
    !real,dimension(extraction_index)::fractionx
    !integer :: i,j
    !real, dimension(0:extraction_index):: depth_cm
    !
    !real :: f,g,h,max_depth
    !real :: total_area
    !
    !
    !f = 0.7
    !g=2.
    !h=0.9
    !max_depth = 2.
    !
    !!*******calculate cumulative depth in cm at each delta x *********
    !
    !depth_cm(0)= 0.
    !depth_cm(1:extraction_index) = depth(1:extraction_index)*100.  !extraction relation is in
    !
    !!**********************************************************
    !!integral between
    !total_area = f*max_depth/(g*max_depth+h)/h
    !
    !!this gives array of fractional area for each delta x
    !fractionx = 0.1578*f*(depth_cm(1:extraction_index)- depth_cm(0:extraction_index-1))   &
    !           /(g*depth_cm(1:extraction_index)+h)/(g*depth_cm(0:extraction_index-1)+h)  &
    !           /total_area
    !
    !extraction_array =0.
    !
    !
    !
    !
    !
    !forall (i=1:extraction_index,j=1:num_records,runoff(j)>0.)
    !    extraction_array(i,j) = runoff(j)*fractionx(i)
    !end forall
    !
    !end subroutine process_extraction
    !
    !
    !
    !
end module runoff_extraction_module