!*************************************************************
!*************************************************************
!Layne Price, University of Auckland, 3/9/2012
!*************************************************************
!*************************************************************

!This is a module that enables one to find portions of a data set, "success," that have a high density and are separated from any points in data set, "fail," by a predefined distance.  This is roughly-similar to some cluster finding techniques from standard data mining texts, e.g. DBSCAN.  However, it is implemented to be quick, rather than precise.  It will not be perfect in finding clusters, but will be able to identify if areas of overdensity exist.

!USAGE:
!FUNCTION	get_insulatedcorepts(succ,fail,metric,sizeneighb,dencrit)
!Returns all points of a table "succ" that are of a sufficient density "dencrit" and are not near points in the "fail" set.  This should work in any dimension, and also for sets embedded in a higher dimension.  This takes one of the metrics, defined below, as an argument. 

!If you only want to find sufficiently clustered points in "succ", without reference to a fail set, then input "fail" as one point sufficiently removed from the others.  (Input as null?? Not tested...)


MODULE fcluster
USE sorters
IMPLICIT NONE

CONTAINS



!Function which returns all the good points.  Sizeneigh and dencrit are optional arguments that have the defaults set to unity.

subroutine get_insulatedcorepts(core,succ,fail,metric,sizeneighb,dencrit)
implicit none

	double precision, dimension(:,:), intent(in) :: succ, fail
	interface
		pure function metric(pt1,pt2)
			implicit none
			double precision, dimension(:), intent(in) :: pt1, pt2
			double precision :: metric
		end function metric
	end interface
	double precision, dimension(:), intent(in), optional :: sizeneighb
	integer, intent(in), optional :: dencrit
	double precision, dimension(:,:), allocatable :: core
	logical, dimension(size(succ,1)) :: good
	integer :: i, length, cnt, density
	double precision, dimension(size(succ,2)) :: eps

	!Set the optional arguments.
	if (present(dencrit)) then
		density=dencrit
	else
		!Default density to having only one point in eps box.
		density=1
	end if

	if (present(sizeneighb)) then
		eps = sizeneighb
	else
		!Default size of eps set to unity
		eps=1D0
	end if

	!Returns logical vector good with .TRUE. for good points and .FALSE. for bad ones.
	good = goodpts(succ,fail,eps,density,metric)

	!Returns number of .TRUE. elmts in good.
	length = count(good)
	
	allocate(core(length,size(succ,2)))

	cnt=0
	do i=1,size(succ,1)
		if (good(i)) then
			cnt=cnt+1
			core(cnt,:)=succ(i,:)
		end if
	end do
	
end subroutine get_insulatedcorepts

!Functions which returns a Logical vector goodpts that has a .TRUE. for every point in succ that has 1) Neps(succ) > dencrit 2) Neps(fail)=0.

function goodpts(succ,fail,eps,dencrit,metric)
implicit none

	double precision, dimension(:,:), intent(in) :: succ, fail
	double precision, dimension(:), intent(in) :: eps
	integer, intent(in) :: dencrit
	interface
		pure function metric(pt1,pt2)
			implicit none
			double precision, dimension(:), intent(in) :: pt1, pt2
			double precision :: metric
		end function metric
	end interface
	logical, dimension(size(succ,1)) :: goodpts
	integer, dimension(size(succ,1)) :: epsnumb_succ, epsnumb_fail
	integer :: i

	
	!Get number of points in eps-Neigh for each succ point wrt succ set.
	call Neps(epsnumb_succ,succ,succ,eps,metric)
	!Get number of points in eps-Neigh for each succ point wrt fail set.
	call Neps(epsnumb_fail,succ,fail,eps,metric)

	!If more than crit numb of good pts in eps-Neigh & 0 fail points, then good point.
	goodpts(:) = ((epsnumb_succ(:) .ge. dencrit) .and. epsnumb_fail(:)==0)



end function goodpts



!Function which takes two sets, setA and setB, and returns the vector Neps that is the number of points in the epsilon-neighborhood of each element in setA with respect to setB.
subroutine Neps(output, setA, setB, eps, metric)
use omp_lib
implicit none

	double precision, dimension(:,:), intent(in) :: seta, setb
	double precision, dimension(:), intent(in) :: eps
	interface
		pure function metric(pt1,pt2)
			implicit none
			double precision, dimension(:), intent(in) :: pt1, pt2
			double precision :: metric
		end function metric
	end interface
	integer, dimension(size(seta,1)) :: output
	integer :: i, j

	!Parallelize

	!$OMP PARALLEL DEFAULT(NONE) &
	!$OMP& SHARED(setA,setB,eps, output)
	!$OMP DO SCHEDULE(STATIC)
	do i=1,SIZE(output)
		output(i)=eps_neigh(setA(i,:),setB, eps, metric)
	end do
	!$OMP END DO
	!$OMP END PARALLEL

end subroutine Neps


!Function to find the epsilon neighborhood of a point with respect to a set of D-dimensional points given in an NxD array (that has been *heapsorted*) and a given metric.
pure integer function eps_neigh(pt, set, eps, metric)
implicit none

	double precision, dimension(:), intent(in) :: eps
	interface
		pure function metric(pt1,pt2)
			implicit none
			double precision, dimension(:), intent(in) :: pt1, pt2
			double precision :: metric
		end function metric
	end interface
	double precision, dimension(:,:), intent(in) :: set
	double precision, dimension(:), intent(in) :: pt
	double precision :: x
	integer :: i, low, high

	!Find number of points in set that are within eps in ith dimension.
	low= location(set,pt(1)-eps(1))
	if (low==0) low=1
	high= location(set,pt(1)+eps(1))

	!Run through array.
	eps_neigh=0
	do i=low,high
		!Quick scan of chosen point.
		if(any(abs(set(i,:)-pt)>eps)) then
			cycle
			!If pass, then call in_box.
		else if (in_box(set(i,:),pt,eps,metric)) then
				eps_neigh=eps_neigh+1
		end if
	end do


end function eps_neigh


!Function to determine whether a point "pt1" is in the "eps" box of "pt2".
pure function in_box(pt1,pt2,eps,metric)
implicit none

	double precision, dimension(:), intent(in) :: pt1, pt2, eps
	logical :: in_box
	interface
		pure function metric(pt1,pt2)
			implicit none
			double precision, dimension(:), intent(in) :: pt1, pt2
			double precision :: metric
		end function metric
	end interface
	double precision :: dist
	integer :: i

	!Rescale the points so that eps-ball has unit radius.
	!Determine if in unit radius.
	in_box= (metric(pt1/eps,pt2/eps) .le. 1D0)

end function in_box

!Subroutine to set eps in each dimension according to the distribution of points in success and fail sets.  This will default eps to 10 times the average distance that one would expect between subsequent points in success.
pure subroutine set_eps(succ,eps,n)
implicit none

	double precision, dimension(:,:), intent(in) :: succ
	double precision, dimension(size(succ,2)), intent(out) :: eps
	double precision, dimension(size(succ,2)) :: maxim, minim
	double precision :: maxwork, minwork
	integer, optional, intent(in) :: n
	integer :: i, j

	!Find min and max of succ in each dimn.
	do i=1,size(succ,2)
		maxwork=succ(1,i)
		minwork=maxwork
		if (size(succ,1).le.1) cycle
		do j=2,size(succ,1)
			if (succ(j,i)>maxwork) maxwork=succ(j,i)
			if (succ(j,i)<minwork) minwork=succ(j,i)
		end do
		maxim(i)=maxwork
		minim(i)=minwork
	end do
	!Find avg spacing and set eps.
	if(present(n)) then
		eps=dble(n)*((maxim-minim)/dble(size(succ,1)))
	else
		eps=10D0*((maxim-minim)/dble(size(succ,1)))
	end if

end subroutine set_eps



!*******************************************************
!FUNCTION which will search using bisection a table table(n,m)  which has been previously ordered by its first column.  Given a value X it will return the value J such that X is between table(J,1) and table(J+1,1).  J=0 or J=SIZE(table,1) if X is out of range.

!This function is taken almost verbatim from Numerical Recipes pg 90.

!NOTE: there is a similar subroutine in sorters_d that does this too, so this one is called "location" where that one is called "locate".

pure integer function location(table,x)
implicit none

	double precision, dimension(:,:),intent(in) :: table
	double precision, intent(in) :: x
	integer :: jl, ju, n, jm
	integer :: i
	integer :: j

	n = size(table,1)

	jl = 0
	ju = n+1
	do while (ju-jl>1)
		jm=(ju+jl)/2
		if((table(n,1)> table(1,1)) .eqv. (x > table(jm,1))) then
			jl = jm
		else
			ju=jm
		end if
	end do
	j = jl
	location = j

end function location

!********************************************************
!Different metrics.  Computes distance between two D-dimensional points expressed as a vector.

pure double precision function manhattan(pt1,pt2)
implicit none

	double precision, dimension(:), intent(in) :: pt1, pt2

	manhattan=sum(abs(pt1-pt2))

end function manhattan

pure double precision function euclidean(pt1,pt2)
implicit none

	double precision, dimension(:), intent(in) :: pt1, pt2

	euclidean=sqrt(sum((pt1-pt2)*(pt1-pt2)))

end function euclidean

pure double precision function dist_n(pt1,pt2,n)
implicit none

	double precision, dimension(:), intent(in) :: pt1, pt2
	integer, intent(in) :: n

	dist_n = (sum(abs(pt1-pt2)**n))**(1d0/dble(n))


end function dist_n

!********************************************************



!Subroutine that reads the success set.
subroutine read_succ(success, fname, formt)
implicit none

	double precision, dimension(:,:), intent(inout) :: success
	character(len=*), intent(in) :: fname
	character(len=*), optional, intent(in) :: formt
	integer :: check, i, j, u
	integer :: length_s, width_s
	integer :: stat

	length_s = size(success,1)
	width_s = size(success,2)

	u=31415927

	if (present(formt)) then
		open(unit=u,status='old',file=fname,form=formt)
	else
		open(unit=u,status='old',file=fname)
	end if

	check = 0
	do i=1,length_s+1
		check = check + 1
		if (present(formt)) then
			read(u,iostat=stat) (success(i,j), j=1,width_s)
			if (is_iostat_end(stat)) exit
		else
			read(u,iostat=stat,fmt=*) (success(i,j), j=1,width_s)
			if (is_iostat_end(stat)) exit
		end if
		if(check==size(success,1)) exit
	end do

	close(unit=u)

end subroutine read_succ


!Subroutine that reads the fail set.
subroutine read_fail(fail, fname,formt)
implicit none

	double precision, dimension(:,:), intent(inout) :: fail
	character(len=*), intent(in) :: fname
	character(len=*), optional, intent(in) :: formt
	integer :: check, i, j, u
	integer :: length_f, width_f

	length_f = size(fail,1)
	width_f = size(fail,2)

	u=31415927

	if (present(formt)) then
		open(unit=u,status='old',file=fname,form=formt)
	else
		open(unit=u,status='old',file=fname)
	end if

	check = 0
	do i=1,length_f+1
		check = check + 1
		if (present(formt)) then
			read(u,end=10) (fail(i,j), j=1,width_f)
		else
			read(u,end=10,fmt=*) (fail(i,j), j=1,width_f)
		end if
		if(check==size(fail,1)) exit
	end do

10	close(unit=u)


end subroutine read_fail


!Subroutine that generates uniform noise to compare our success and fail sets against.  Normalized versus the maximum and minimum elements in each dimension for the comparison set.
subroutine noise(set, mnm, maxm)
use rng
implicit none

	double precision, dimension(:,:), intent(out) :: set
	double precision, dimension(:), intent(in) :: mnm, maxm
	integer :: i

	!Get random numbers.
	call random_number(set)

	!Normalize it.
	do i=1,size(set,2)
		set(:,i)=set(:,i)*(maxm(i)-mnm(i)) + mnm(i)
	end do	


end subroutine noise




end module fcluster



















!Pullback of the Euclidean metric onto the equal energy density constraint surface.
!DOUBLE PRECISION FUNCTION pullback_eucl(pt1,pt2)
!IMPLICIT NONE

!	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
!	DOUBLE PRECISION, DIMENSION(SIZE(pt1),SIZE(pt2)) :: metric
!	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: diff

!	ALLOCATE(diff(SIZE(pt1)))
!	diff=pt1-pt2

!	metric = 0D0

!	pullback_eucl=0D0

	

!END FUNCTION pullback_eucl



!******************************************************************



!************************************************************
!NOT READY FOR PRIMETIME: Currently implemented in Python.
!************************************************************
!FUNCTION that returns the number of points, N, that are in the eps-neighb of the nearest point in the data set to an argument point, pt.  This is basically a nearest-neighbor interpolation on the data set with respect to sample density.  Will be used for MCMC once a cluster has been found at a suitable value for eps.

!Needs the set to be pre-sorted by first column.  REQUIRES: numbeps=Neps(set,set,eps,metric) as input!!!!

!INTEGER FUNCTION nearest_neighbor_density(pt,set,numbeps,eps,metric)
!IMPLICIT NONE

!	DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: set
!	INTERFACE
!		FUNCTION metric(pt1,pt2)
!			IMPLICIT NONE
!			DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
!			DOUBLE PRECISION :: metric
!		END FUNCTION metric
!	END INTERFACE
!	DOUBLE PRECISION, INTENT(IN) :: eps
!	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt
!	INTEGER, INTENT(IN) :: overdens
!	INTEGER, DIMENSION(:), INTENT(IN) :: numbeps
!	INTEGER :: i, low, high

	!Find the nearest neighbor to pt.
!	low = location(set,pt(1))
!	IF (low==0) low=1

!END FUNCTION nearest_neighbor_density





























