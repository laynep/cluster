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


module fcluster
  use types, only : dp
  use features, only : newunit
  implicit none

contains

!Function which returns all the good points.  Sizeneigh and dencrit are optional arguments that have the defaults set to unity.
subroutine get_insulatedcorepts(core,succ,fail,metric,sizeneighb,dencrit)
	implicit none

	real(dp), dimension(:,:), intent(in) :: succ, fail
	interface
		pure function metric(pt1,pt2)
		  use types, only : dp
    	implicit none
		  real(dp), dimension(:), intent(in) :: pt1, pt2
			real(dp) :: metric
		end function metric
	end interface
	real(dp), dimension(:), intent(in), optional :: sizeneighb
	integer, intent(in), optional :: dencrit
	real(dp), dimension(:,:), allocatable :: core
	logical, dimension(size(succ,1)) :: good
	integer :: i, length, cnt, density
	real(dp), dimension(size(succ,2)) :: eps

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
		eps=1_dp
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

!Function which returns a Logical vector goodpts that has a .TRUE. for every point in succ that has 1) Neps(succ) > dencrit 2) Neps(fail)=0.
function goodpts(succ,fail,eps,dencrit,metric)
	implicit none

	real(dp), dimension(:,:), intent(in) :: succ, fail
	real(dp), dimension(:), intent(in) :: eps
	integer, intent(in) :: dencrit
	interface
		pure function metric(pt1,pt2)
		  use types, only : dp
    	implicit none
			real(dp), dimension(:), intent(in) :: pt1, pt2
			real(dp) :: metric
		end function metric
	end interface
	logical, dimension(size(succ,1)) :: goodpts
	integer, dimension(size(succ,1)) :: epsnumb_succ, epsnumb_fail
  logical :: failstop

  failstop=.true.
	
	!Get number of points in eps-Neigh for each succ point wrt succ set.
	call Neps(epsnumb_succ,succ,succ,eps,metric)
	!Get number of points in eps-Neigh for each succ point wrt fail set.
  !Set flag so that will stop counting when it hits one failure points.
	call Neps(epsnumb_fail,succ,fail,eps,metric, failstop)

	!If more than crit numb of good pts in eps-Neigh & 0 fail points, then good point.
	goodpts(:) = ((epsnumb_succ(:) .ge. dencrit) .and. epsnumb_fail(:)==0)

end function goodpts

!Function which takes two sets, setA and setB, and returns the vector Neps that is the number of points in the epsilon-neighborhood of each element in setA with respect to setB.  Also has an optional input "failstop" that will stop the counting of nearby points at one.  We can then use this when comparing against the fail set: if there's one or more points adjacent to a successful point, we can stop counting and discard that point from the protected cluster.
subroutine Neps(output, setA, setB, eps, metric, failstop)
	use omp_lib
	implicit none

	real(dp), dimension(:,:), intent(in) :: seta, setb
	real(dp), dimension(:), intent(in) :: eps
	interface
		pure function metric(pt1,pt2)
		  use types, only : dp
    	implicit none
			real(dp), dimension(:), intent(in) :: pt1, pt2
			real(dp) :: metric
		end function metric
	end interface
	integer, dimension(size(seta,1)) :: output
	integer :: i
  logical, intent(in), optional :: failstop

	!Parallelize

	!$OMP PARALLEL DEFAULT(NONE) &
	!$OMP& SHARED(setA,setB,eps, output, failstop)
	!$OMP DO SCHEDULE(STATIC)
	do i=1,size(output)
    if (present(failstop)) then
  		output(i)=eps_neigh(setA(i,:),setB, eps, metric, failstop)
    else
    	output(i)=eps_neigh(setA(i,:),setB, eps, metric)
    end if
	end do
	!$OMP END DO
	!$OMP END PARALLEL

end subroutine Neps


!Function to find the epsilon neighborhood of a point with respect to a set of D-dimensional points given in an NxD array (that has been *heapsorted*) and a given metric.  Also has an optional input "failstop" that will stop the counting of nearby points at one.  We can then use this when comparing against the fail set: if there's one or more points adjacent to a successful point, we can stop counting and discard that point from the protected cluster.
integer function eps_neigh(pt, set, eps, metric, failstop)
	implicit none

	real(dp), dimension(:), intent(in) :: eps
	interface
		pure function metric(pt1,pt2)
		  use types, only : dp
    	implicit none
			real(dp), dimension(:), intent(in) :: pt1, pt2
      !Get the corepoints at the minimum value for epsilong.
			real(dp) :: metric
		end function metric
	end interface
	real(dp), dimension(:,:), intent(in) :: set
	real(dp), dimension(:), intent(in) :: pt
	integer :: i, low, high
  logical, intent(in), optional :: failstop

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
		!If pass, then rescale points so that an epsilon ball has unit radius and
    !see if the two points are within a unit of each other.
		else if (metric(pt/eps,set(i,:)/eps) .le. 1_dp) then
				eps_neigh=eps_neigh+1
        if (present(failstop)) exit
		end if
	end do

end function eps_neigh


!Subroutine to set eps in each dimension according to the distribution of points in success and fail sets.  This will default eps to 10 times the average distance that one would expect between subsequent points in success.
pure subroutine set_eps(succ,eps,n)
	implicit none

	real(dp), dimension(:,:), intent(in) :: succ
	real(dp), dimension(size(succ,2)), intent(out) :: eps
	real(dp), dimension(size(succ,2)) :: maxim, minim
	real(dp) :: maxwork, minwork
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
		eps=10_dp*((maxim-minim)/dble(size(succ,1)))
	end if

end subroutine set_eps

!*******************************************************
!FUNCTION which will search using bisection a table table(n,m)  which has been previously ordered by its first column.  Given a value X it will return the value J such that X is between table(J,1) and table(J+1,1).  J=0 or J=SIZE(table,1) if X is out of range.

!This function is taken almost verbatim from Numerical Recipes pg 90.

!NOTE: there is a similar subroutine in sorters_d that does this too, so this one is called "location" where that one is called "locate".
pure integer function location(table,x)
	implicit none

	real(dp), dimension(:,:),intent(in) :: table
	real(dp), intent(in) :: x
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

pure real(dp) function manhattan(pt1,pt2)
	implicit none

	real(dp), dimension(:), intent(in) :: pt1, pt2

	manhattan=sum(abs(pt1-pt2))

end function manhattan

pure real(dp) function euclidean(pt1,pt2)
	implicit none

	real(dp), dimension(:), intent(in) :: pt1, pt2

	euclidean=sqrt(sum((pt1-pt2)*(pt1-pt2)))

end function euclidean

pure real(dp) function dist_n(pt1,pt2,n)
	implicit none

	real(dp), dimension(:), intent(in) :: pt1, pt2
	integer, intent(in) :: n

	dist_n = (sum(abs(pt1-pt2)**n))**(1_dp/dble(n))


end function dist_n

!********************************************************


!Subroutine that reads the success set.
subroutine read_succ(success, fname, formt)
	implicit none

	real(dp), dimension(:,:), intent(inout) :: success
	character(len=*), intent(in) :: fname
	character(len=*), optional, intent(in) :: formt
	integer :: check, i, j, u
	integer :: length_s, width_s
	integer :: stat

	length_s = size(success,1)
	width_s = size(success,2)

	!u=31415927

	if (present(formt)) then
		open(unit=newunit(u),status='old',file=fname,form=formt)
	else
		open(unit=newunit(u),status='old',file=fname)
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
subroutine read_fail(fail, fname, formt)
	implicit none

	real(dp), dimension(:,:), intent(inout) :: fail
	character(len=*), intent(in) :: fname
	character(len=*), optional, intent(in) :: formt
	integer :: check, i, j, u
	integer :: length_f, width_f

	length_f = size(fail,1)
	width_f = size(fail,2)

	!u=31415927

	if (present(formt)) then
		open(unit=newunit(u),status='old',file=fname,form=formt)
	else
		open(unit=newunit(u),status='old',file=fname)
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
	implicit none

	real(dp), dimension(:,:), intent(out) :: set
	real(dp), dimension(:), intent(in) :: mnm, maxm
	integer :: i

	!Get random numbers.
	call random_number(set)

	!Normalize it.
	do i=1,size(set,2)
		set(:,i)=set(:,i)*(maxm(i)-mnm(i)) + mnm(i)
	end do	


end subroutine noise

!A subroutine that prints the corepoints and deallocates the insulatedpoints
!array.
subroutine print_corepoints(success, insulatedpts, printing,eps,k)
  implicit none

  logical, intent(in) :: printing
  real(dp), dimension(:), intent(in) :: eps
  integer, intent(in), optional :: k
  real(dp), dimension(:,:), allocatable, intent(inout) :: insulatedpts, success
  integer :: u, i, j, kk
	character(len=18) :: corename
  real :: ratio

  if (present(k)) then
    kk=k
  else
    kk=1
  end if

  if (printing) print*,"Printing core points for eps=",eps
  !Name file.
	write(corename,'(a,i4.4,a)')'corepoints',(kk+1),".bin"
	!Open file.
	open(unit=newunit(u),file=corename,form='unformatted')
	do i=1,size(insulatedpts,1)
		write(unit=u), (insulatedpts(i,j),j=1,size(insulatedpts,2))
!    print*,(insulatedpts(i,j),j=1,size(insulatedpts,2))
	end do
	close(unit=u)		
	!Ratio of points in cluster to points in success.
  ratio=real(size(insulatedpts))/real(size(success))
	if (printing) print*,ratio, "Percent of total are core points"

end subroutine print_corepoints

!A function that removes a subset from a set.  Copies the original
!set, which makes this pretty memory intensive for large arrays.
subroutine complement(comp, set, subset, tol)
	use sorters, only : heapsort, locate
	implicit none

	real(dp), dimension(:,:), intent(in) :: set
	real(dp), dimension(:,:), intent(inout) :: subset
	real(dp), optional, intent(in) :: tol
	real(dp), dimension(:,:), allocatable, intent(out) :: comp
	integer :: i, j, k, counter, start
	logical :: same
	logical, dimension(size(set,1)) :: take
	real(dp) :: dt

	if (present(tol)) then
		dt=tol
	else
		dt=1e-10_dp
	end if

	!Presort the subset so can use locate later.
	call heapsort(subset)

	!Parallelize

	!$OMP PARALLEL &
	!$OMP& SHARED(set,subset,take)&
	!$OMP& PRIVATE(same)
	!$OMP DO 

doi:	do i=1,size(set,1)
    		call locate(subset,set(i,1),start)
        if (start==0) start=1
doj:		do j=start,size(subset,1)
	  		same=.false.
dok:			do k=1,size(subset,2)
	  			if (abs(set(i,k)-subset(j,k))<dt) then
	  				same=.true.
	  			else
	  				same=.false.
	  				exit dok
	  			end if
	  		end do dok	
	  		if (same) then
	  			!Don't take this row for complement.
	  			take(i)=.false.
	  			exit doj
	  		end if
	  	end do doj
	  	if (.not. same) then
	  		take(i)=.true.
	  	end if
  	end do doi

	!$OMP END DO
	!$OMP END PARALLEL

  !Make the complement array.
  if (.not. allocated(comp)) allocate(comp(count(take),size(set,2)))

	counter=0
	do i=1,size(set,1)
			if (take(i)) then
			  counter=counter+1
        if (counter>size(comp,1)) exit
        comp(counter,:)=set(i,:)
		end if
	end do

end subroutine complement

end module fcluster

