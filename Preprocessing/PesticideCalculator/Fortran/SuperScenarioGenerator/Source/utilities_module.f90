module utilities_module

contains

	subroutine cumulative_depth (n, delta_x, depth)
		!generate a cumulative vector depth, given a vector
		!of incremental values delta_x

	    implicit none
		integer, intent(in) :: n	!size of vector delta_x
		real,dimension(n),intent(in) :: delta_x
		real,dimension(n),intent(out) ::depth  !cumulative depth array
		integer :: i

        !*** build array of cumulative depths ***
        depth = 0.
		depth(1) = delta_x(1)
        do i=2,n
            depth(i) = delta_x(i) + depth(i-1)
        end do
	end subroutine cumulative_depth



	!*********************************************************************
	 pure integer function find_node(n,depth,target_depth)
       !Given a vector depth of size n that contains cumulative depths,
	   !this routine gives the nearest node corresponding to target_depth
	   !minimum is one node
	   !maximum is n

		implicit none
		integer, intent(in) :: n  !size of vector depth
		real,intent(in),dimension(n) ::depth  !cumulative depth array
		real,intent(in) :: target_depth				!target depth for which node is to be found

		integer :: node					!the node output
		integer :: node_1(1)

		!***************************************************************
		if (target_depth >= depth(n)) then
			node = n
		else if (target_depth <= depth(1)) then
			node = 1
		else
			node_1 = minloc(depth,mask= depth >=target_depth )
			node = node_1(1)  
			if (depth(node) - target_depth  > target_depth-depth(node-1)) node = node-1
		endif
		
        !****************************************************************************   

		find_node = node
	end function find_node



   !*****************************************************************************
   pure elemental integer function jd (YEAR,MONTH,DAY)
     !calculate the days since 1/1/1900 given year,month, day, from Fliegel and van Flandern (1968)
	 !Fliegel, H. F. and van Flandern, T. C. (1968). Communications of the ACM, Vol. 11, No. 10 (October, 1968). 

     implicit none
     integer, intent(in) :: year,month,day

      JD= day-32075+1461*(year+4800+(month-14)/12)/4+367*(month-2-(month-14)/12*12) /12-3*((year+4900+(month-14)/12)/100)/4 -2415021

    end function jd
   !*****************************************************************************
 


	pure subroutine get_date (date1900, YEAR,MONTH,DAY)
!     computes THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY) given days since 1900
      implicit none

      integer,intent(out) :: YEAR,MONTH,DAY

	  integer,intent(in) :: date1900
	  integer :: L,n,i,j

      L= 2483590 + date1900

      n= 4*L/146097

      L= L-(146097*n+3)/4
      I= 4000*(L+1)/1461001
      L= L-1461*I/4+31
      J= 80*L/2447

      day= L-2447*J/80
      L= J/11
      month = J+2-12*L
      year = 100*(N-49)+I+L

    !  YEAR= I
    !  MONTH= J
    !  DAY= K

      end subroutine get_date

end module utilities_module