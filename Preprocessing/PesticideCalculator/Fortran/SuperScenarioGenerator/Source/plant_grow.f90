module plant_growth_module
    implicit none
    integer,private::status

contains

subroutine plant_growth
   !********Description ***********************************************
   !linear growth between emergence (height = 0) to maturity (height = 1)
   !This subroutine returns the vector plant_factor which gives the daily
	!fraction of full growth. There are no restrictions on the number of
	!plantings or whether they cross calendar years. Evergreens, for instance
	!could have just one planting before the simulation and harvest after the
	!simulation. 
	!*****************************************************************
    use variables_noninputs, ONLY:num_records,numberOfYears  ,& 
                                emergence, &  !input
                                maturity,  &  !input
                                harvest,  &   !input  
                                plant_factor  !output
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! local variables
	integer :: i,j
	integer :: a_limit,b_limit,a_diff,b_diff,b_limit_2,g_limit
	integer :: d_size
	integer,dimension(numberOfYears) :: me, c
	real   ,dimension(numberOfYears) :: rme
	real,allocatable,dimension(:) :: d,fract
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
	plant_factor = 0.
	c=0

	me = maturity- emergence
	c = me+1     !number of days from emer to maturity (inclde first and last)
	rme =real(me)

	d_size = maxval(c)
	allocate(d(d_size),fract(d_size), STAT = status)

	forall(j=1:d_size) d(j) = real(j-1)  !create an array from zero to number of days between emerg and maturity

	!min and max functions are in here in case plant growth dates exceed simulation dates
	!This loop maps the emergence,maturity, harvest interval onto the same time refernce as the metfile
	do i=1,numberOfYears

		fract(1:c(i)) = d(1:c(i))/rme(i)       !array of fractional growth (i.e., 1,0.25,0.5,0.75,1.0)

		a_limit = max(emergence(i),1)          !index for plant factor cannot be negative, must start at 1
		b_limit = min(maturity(i),num_records) !index for plant factor cannot be >num_records

		a_diff = a_limit-emergence(i)+1        !amount to adjust index by
		b_diff = maturity(i) -b_limit+1

		if (a_limit <= num_records) then       !preclude the case where crop is out of range of simulation
        !these are the plant factors between emergence and harvest
        plant_factor(a_limit:b_limit) = fract(a_diff:b_diff)
		end if

		!do maturity to harvest (height =1)
		b_limit_2 = max(1,maturity(i))
		g_limit = min(num_records,harvest(i))

		if (g_limit >= 1) then                 !preclude the case where crop is out of range of simulation
        plant_factor(b_limit_2:g_limit)= 1
		end if

	end do
	
		
	deallocate (d,fract, STAT= status)

end subroutine plant_growth


end module plant_growth_module