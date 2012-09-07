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

	NAMELIST /tablel/ length_s, length_f, width_s, width_f

	!Reads file sizes from input file "setsizes.txt".
	OPEN(unit=1000, file="setsizes.txt", status="old", delim="apostrophe")
	READ(unit=1000, nml=tablel)
	CLOSE(unit=1000)
	ALLOCATE(success(length_s,width_s),fail(length_f,width_f))

	!Read succ and fail sets from file.
!call read_succ(success, "success.bin","unformatted")
!call read_fail(fail, "fail.bin","unformatted")
	call read_succ(success, "totalsucc.bin","unformatted")
	call read_fail(fail, "totalfail.bin","unformatted")

	!Sort the success and fail sets by value in first column.
	call heapsort(success)
	call heapsort(fail)

	PRINT*,"Getting core points"
	eps=.5D0
	dencrit=2
	insulatedpts=get_insulatedcorepts(success,fail,euclidean,eps,dencrit)
if(allocated(insulatedpts)) then
	print*,size(insulatedpts,1),size(insulatedpts,2)
else
	print*,"insulatedpts isn't allocated??"
end if
	PRINT*,"Printing core points"


	open(unit=3,file="corepoints.bin",status="new",form='UNFORMATTED')
	DO i=1,SIZE(insulatedpts,1)
print*, "printing",i, (insulatedpts(i,j),j=1,SIZE(insulatedpts,2))
!		WRITE(unit=3), (insulatedpts(i,j),j=1,SIZE(insulatedpts,2))
	END DO
!	CLOSE(unit=3)

	if(allocated(insulatedpts)) deallocate(insulatedpts)
	if(allocated(success)) deallocate(success)
	if(allocated(fail)) deallocate(fail)

	


END PROGRAM findclusters

































