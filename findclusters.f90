PROGRAM findclusters
USE fcluster
USE sorters
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: success, fail
	INTEGER :: i,j,k,check
	INTEGER :: length_s, length_f, width_s, width_f
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: insulatedpts
	DOUBLE PRECISION :: eps
	INTEGER :: dencrit
	logical :: printing

	NAMELIST /tablel/ length_s, length_f, width_s, width_f

	printing = .true.

	!Reads file sizes from input file "setsizes.txt".
	OPEN(unit=1000, file="setsizes.txt", status="old", delim="apostrophe")
	READ(unit=1000, nml=tablel)
	CLOSE(unit=1000)
	ALLOCATE(success(length_s,width_s),fail(length_f,width_f))

	!Read succ and fail sets from file.
	if (printing) print*, "Reading files."
	call read_succ(success, "totalsucc.bin","unformatted")
	call read_fail(fail, "totalfail.bin","unformatted")

	!Sort the success and fail sets by value in first column.
	if (printing) print*, "Sorting files."
	call heapsort(success)
	call heapsort(fail)

	if (printing) PRINT*,"Getting core points."
	eps=.5D0
	dencrit=2
	call get_insulatedcorepts(insulatedpts,success,fail,euclidean,eps,dencrit)

	if (printing) then
		PRINT*,"Printing core points."
		open(unit=3,file="corepoints.bin",form='UNFORMATTED')
		DO i=1,SIZE(insulatedpts,1)
!print*, "printing",i, (insulatedpts(i,j),j=1,SIZE(insulatedpts,2))
		WRITE(unit=3), (insulatedpts(i,j),j=1,SIZE(insulatedpts,2))
		END DO
		CLOSE(unit=3)
	end if

	if(allocated(insulatedpts)) deallocate(insulatedpts)
	if(allocated(success)) deallocate(success)
	if(allocated(fail)) deallocate(fail)



END PROGRAM findclusters

































