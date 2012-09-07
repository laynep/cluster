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

FUNCTION get_insulatedcorepts(succ,fail,metric,sizeneighb,dencrit)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: succ, fail
	INTERFACE
		pure FUNCTION metric(pt1,pt2)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
			DOUBLE PRECISION :: metric
		END FUNCTION metric
	END INTERFACE
	DOUBLE PRECISION, INTENT(IN), OPTIONAL :: sizeneighb
	INTEGER, INTENT(IN), OPTIONAL :: dencrit
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: get_insulatedcorepts
	LOGICAL, DIMENSION(:), ALLOCATABLE :: good
	INTEGER :: i, length, cnt, density
	DOUBLE PRECISION :: eps

	!Set the optional arguments.
	IF (PRESENT(dencrit)) THEN
		density=dencrit
	ELSE
		!Default density to having only one point in eps box.
		density=1
	END IF

	IF (PRESENT(sizeneighb)) THEN
		eps = sizeneighb
	ELSE
		!Default size of eps set to unity
		eps=1D0
		!eps = ((10./DBLE(SIZE(succ,1)))**(1D0/DBLE(SIZE(succ,2))))*&
		!&(MAXVAL(succ)-MINVAL(succ))
	END IF

	!Returns logical vector good with .TRUE. for good points and .FALSE. for bad ones.
	good = goodpts(succ,fail,eps,density,metric)

	!Returns number of .TRUE. elmts in good.
	length = COUNT(good)
	
	ALLOCATE(get_insulatedcorepts(length,SIZE(succ,2)))

	cnt=0
	DO i=1,SIZE(succ,1)
		IF (good(i)) THEN
			cnt=cnt+1
			get_insulatedcorepts(cnt,:)=succ(i,:)
		END IF
	END DO

	IF (ALLOCATED(good)) DEALLOCATE(good)
	
END FUNCTION get_insulatedcorepts

!Functions which returns a Logical vector goodpts that has a .TRUE. for every point in succ that has 1) Neps(succ) > dencrit 2) Neps(fail)=0.

FUNCTION goodpts(succ,fail,eps,dencrit,metric)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: succ, fail
	DOUBLE PRECISION, INTENT(IN) :: eps
	INTEGER, INTENT(IN) :: dencrit
	INTERFACE
		pure FUNCTION metric(pt1,pt2)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
			DOUBLE PRECISION :: metric
		END FUNCTION metric
	END INTERFACE
	LOGICAL, DIMENSION(SIZE(succ,1)) :: goodpts
!	integer, dimension(size(succ,1)) :: epsnumb_succ, epsnumb_fail
	INTEGER, DIMENSION(:), ALLOCATABLE :: epsnumb_succ, epsnumb_fail
	INTEGER :: i

	allocate (epsnumb_succ(size(succ,1)),epsnumb_fail(size(succ,1)))

	!Get number of points in eps-Neigh for each succ point wrt succ set.
	epsnumb_succ=Neps(succ,succ,eps,metric)
	!Get number of points in eps-Neigh for each succ point wrt fail set.
	epsnumb_fail=Neps(succ,fail,eps,metric)

	!If more than critical numb of good points in eps-Neigh and zero fail points, then consider this to be a good point.
	goodpts(:) = ((epsnumb_succ(:) .ge. dencrit) .AND. epsnumb_fail(:)==0)

	IF (ALLOCATED(epsnumb_succ)) DEALLOCATE(epsnumb_succ)
	IF (ALLOCATED(epsnumb_fail)) DEALLOCATE(epsnumb_fail)

END FUNCTION goodpts



!Function which takes two sets, setA and setB, and returns the vector Neps that is the number of points in the epsilon-neighborhood of each element in setA with respect to setB.
FUNCTION Neps(setA, setB, eps, metric)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: setA, setB
	DOUBLE PRECISION, INTENT(IN) :: eps
	INTERFACE
		pure FUNCTION metric(pt1,pt2)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
			DOUBLE PRECISION :: metric
		END FUNCTION metric
	END INTERFACE
	INTEGER, DIMENSION(SIZE(setA,1)) :: Neps
	INTEGER :: i, j
	INTEGER :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM


	!Parallelize

	!$OMP PARALLEL DEFAULT(NONE) &
	!$OMP& SHARED(setA,setB,eps,Neps)
	!$OMP DO SCHEDULE(STATIC)

	DO i=1,SIZE(setA,1)
		Neps(i)=eps_neigh(setA(i,:),setB, eps, metric)
	END DO

	!$OMP END DO
	!$OMP END PARALLEL

END FUNCTION Neps


!Function to find the epsilon neighborhood of a point with respect to a set of D-dimensional points given in an NxD array (that has been *heapsorted*) and a given metric.
pure INTEGER FUNCTION eps_neigh(pt, set, eps, metric)
IMPLICIT NONE

	DOUBLE PRECISION, INTENT(IN) :: eps
	INTERFACE
		pure FUNCTION metric(pt1,pt2)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
			DOUBLE PRECISION :: metric
		END FUNCTION metric
	END INTERFACE
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: set
	double precision, dimension(:), intent(in) :: pt
	DOUBLE PRECISION :: x
	DOUBLE PRECISION, DIMENSION(SIZE(pt)) :: pt1, pt2
	INTEGER :: i, low, high

	pt1=pt

	!Find number of points in set that are within eps in ith dimension.
	low= location(set,pt(1)-eps)
	IF (low==0) low=1
	high= location(set,pt(1)+eps)

	eps_neigh=0
	DO i=low,high
		pt2=set(i,:)
		IF (metric(pt1,pt2) .le. eps) THEN
			eps_neigh=eps_neigh+1
		END IF
	END DO

END FUNCTION eps_neigh



!*******************************************************
!FUNCTION which will search using bisection a table table(n,m)  which has been previously ordered by its first column.  Given a value X it will return the value J such that X is between table(J,1) and table(J+1,1).  J=0 or J=SIZE(table,1) if X is out of range.

!This function is taken almost verbatim from Numerical Recipes pg 90.

!NOTE: there is a similar subroutine in sorters_d that does this too, so this one is called "location" where that one is called "locate".

pure INTEGER FUNCTION location(table,x)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:),INTENT(IN) :: table
	DOUBLE PRECISION, INTENT(IN) :: x
	INTEGER :: jl, ju, n, jm
	INTEGER :: i
	INTEGER :: j

	n = SIZE(table,1)

	jl = 0
	ju = n+1
	DO WHILE (ju-jl>1)
		jm=(ju+jl)/2
		IF((table(n,1)> table(1,1)) .EQV. (x > table(jm,1))) THEN
			jl = jm
		ELSE
			ju=jm
		END IF
	END DO
	j = jl
	location = j

END FUNCTION location

!********************************************************
!Different metrics.  Computes distance between two D-dimensional points expressed as a vector.

pure DOUBLE PRECISION FUNCTION manhattan(pt1,pt2)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
	INTEGER :: i

	manhattan = 0D0

	DO i=1,SIZE(pt1)
		manhattan = manhattan + ABS(pt1(i)-pt2(i))
	END DO

	manhattan = ABS(manhattan)

END FUNCTION manhattan

pure DOUBLE PRECISION FUNCTION euclidean(pt1,pt2)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
	INTEGER :: i

	euclidean = 0D0
	DO i=1,SIZE(pt1)
		euclidean = euclidean + (pt1(i)-pt2(i))**2
	END DO

	euclidean = SQRT(euclidean)


END FUNCTION euclidean

pure DOUBLE PRECISION FUNCTION dist_N(pt1,pt2,N)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
	INTEGER, INTENT(IN) :: N
	INTEGER :: i

	dist_N=0D0	
	DO i=1,SIZE(pt1)
		dist_N = dist_N + (ABS(pt1(i)-pt2(i)))**N
	END DO

	dist_N = (dist_N)**(1D0/DBLE(N))

END FUNCTION dist_N

!Pullback of the Euclidean metric onto the equal energy density constraint surface.
DOUBLE PRECISION FUNCTION pullback_eucl(pt1,pt2)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
	DOUBLE PRECISION, DIMENSION(SIZE(pt1),SIZE(pt2)) :: metric
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: diff

	ALLOCATE(diff(SIZE(pt1)))
	diff=pt1-pt2

	metric = 0D0

	pullback_eucl=0D0

	

END FUNCTION pullback_eucl



!******************************************************************

!Function to make a cluster centered at "center."
FUNCTION make_blob(center,N,sigma)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: center
	DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: sigma
	DOUBLE PRECISION :: std
	INTEGER, INTENT(IN) :: N
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: make_blob
	INTEGER :: i,j

	!Optional argument for standard deviation.
	IF (PRESENT(sigma)) THEN
		std = sigma
	ELSE
		std = 1D0
	END IF

	CALL init_random_seed()

	ALLOCATE(make_blob(N,SIZE(center)))

	DO j=1,SIZE(center)
		DO i=1,N
			make_blob(i,j) = normal(center(j),std)
		END DO
	END DO


END FUNCTION make_blob

!************************************************************************
!This returns a normal distribution.  Taken from http://www.sdsc.edu/~tkaiser/f90.html .

function normal(mean,sigma) 
        implicit none
	DOUBLE PRECISION, parameter :: pi = 3.141592653589793239D0
        DOUBLE PRECISION :: normal,tmp, ran1, ran2
        DOUBLE PRECISION ::  mean,sigma
        INTEGER :: flag
        DOUBLE PRECISION :: fac,gsave,rsq,r1,r2
        save flag,gsave
        data flag /0/
        if (flag.eq.0) then
        rsq=2.0D0
        do while(rsq.ge. 1.0D0 .or. rsq.eq. 0.0D0)
		CALL random_number(ran1)
		CALL random_number(ran2)
                r1=2.0D0*ran1-1.0D0
                r2=2.0D0*ran2-1.0D0
                rsq=r1*r1+r2*r2
        enddo
            fac=sqrt(-2.0D0*log(rsq)/rsq)
            gsave=r1*fac
            tmp=r2*fac
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp*sigma+mean
      	RETURN
end function normal

!Subroutine to initialize random_seed according to CPU time.

SUBROUTINE init_random_seed()
IMPLICIT NONE
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37* (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE


!************************************************************
!NOT READY FOR PRIMETIME: Currently implemented in Python.
!************************************************************
!FUNCTION that returns the number of points, N, that are in the eps-neighb of the nearest point in the data set to an argument point, pt.  This is basically a nearest-neighbor interpolation on the data set with respect to sample density.  Will be used for MCMC once a cluster has been found at a suitable value for eps.

!Needs the set to be pre-sorted by first column.  REQUIRES: numbeps=Neps(set,set,eps,metric) as input!!!!

INTEGER FUNCTION nearest_neighbor_density(pt,set,numbeps,eps,metric)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: set
	INTERFACE
		FUNCTION metric(pt1,pt2)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt1, pt2
			DOUBLE PRECISION :: metric
		END FUNCTION metric
	END INTERFACE
	DOUBLE PRECISION, INTENT(IN) :: eps
	DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: pt
!	INTEGER, INTENT(IN) :: overdens
	INTEGER, DIMENSION(:), INTENT(IN) :: numbeps
	INTEGER :: i, low, high

	!Find the nearest neighbor to pt.
	low = location(set,pt(1))
	IF (low==0) low=1
	


	
	

END FUNCTION nearest_neighbor_density


subroutine read_succ(success, fname, formt)
implicit none

	double precision, dimension(:,:), intent(inout) :: success
	character(len=*), intent(in) :: fname
	character(len=*), optional, intent(in) :: formt
	integer :: check, i, j, u
	integer :: length_s, width_s

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
			read(u,end=10) (success(i,j), j=1,width_s)
		else
			read(u,end=10,fmt=*) (success(i,j), j=1,width_s)
		end if
		if(check==size(success,1)) exit
	end do

10	close(unit=u)

end subroutine read_succ



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


END MODULE fcluster





































