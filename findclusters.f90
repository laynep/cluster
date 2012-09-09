program findclusters
use fcluster
use sorters
implicit none

	double precision, dimension(:,:), allocatable :: success, fail
	integer :: i,j,k,check
	integer :: length_s, length_f, width_s, width_f
	double precision, dimension(:,:), allocatable :: insulatedpts
	double precision, dimension(:), allocatable :: eps
	integer :: dencrit
	logical :: printing, auto

	namelist /tablel/ length_s, length_f, width_s, width_f, printing, auto

	!Reads file sizes from input file "setsizes.txt".
	open(unit=1000, file="setsizes.txt", status="old", delim="apostrophe")
	read(unit=1000, nml=tablel)
	close(unit=1000)
	allocate(success(length_s,width_s),fail(length_f,width_f))

	!Read succ and fail sets from file.
	if (printing) print*, "Reading files."
	call read_succ(success, "totalsucc.bin","unformatted")
	call read_fail(fail, "totalfail.bin","unformatted")

	!Sort the success and fail sets by value in first column.
	if (printing) print*, "Sorting files."
	call heapsort(success)
	call heapsort(fail)

	!Get the core points.
	if (printing) PRINT*,"Getting core points."
	!Set eps in each dimn.
	allocate(eps(size(success,2)))
	if (auto) then
		call set_eps(success,eps,10)
		if (printing) print*, "Epsilon is", eps
		dencrit=2	!No other points req in eps-ball.
	else
		eps=.5D0
		dencrit=2	!At least one other point in eps-ball.
	end if
	call get_insulatedcorepts(insulatedpts,success,fail,euclidean,eps,dencrit)

	!Print the core points
	if (printing) print*,"Printing core points."
	open(unit=3,file="corepoints.bin",form='unformatted')
	do i=1,size(insulatedpts,1)
		write(unit=3), (insulatedpts(i,j),j=1,size(insulatedpts,2))
	end do
	close(unit=3)

	if(allocated(eps)) deallocate(eps)
	if(allocated(insulatedpts)) deallocate(insulatedpts)
	if(allocated(success)) deallocate(success)
	if(allocated(fail)) deallocate(fail)

end program findclusters

































