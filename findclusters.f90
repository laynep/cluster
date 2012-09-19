program findclusters
  use fcluster
  use sorters, only : heapsort
  use rng, only : shuffle_cut
  use types, only : dp
  use features, only : newunit
  implicit none

	real(dp), dimension(:,:), allocatable :: success, fail
	integer :: i,j,k, kend, check
	integer :: length_s, length_f, width_s, width_f
	real(dp), dimension(:,:), allocatable :: insulatedpts
	real(dp), dimension(:), allocatable :: eps
	integer :: dencrit, u
	logical :: printing, auto, shuffling, find_min
	character(len=18) :: corename
	real :: ratio
  real(dp) :: energy_scale, mplanck

	namelist /tablel/ length_s, length_f, width_s, width_f, printing, auto, shuffling, find_min
  namelist /phy_param/ energy_scale, mplanck

	!Reads file sizes from input file "setsizes.txt".
	open(unit=newunit(u), file="setsizes.txt", status="old", delim="apostrophe")
	read(unit=u, nml=tablel)
  read(unit=u, nml=phy_param)
	close(unit=u)
	allocate(success(length_s,width_s),fail(length_f,width_f))

	!Read succ and fail sets from file.
	if (printing) print*, "Reading files."
	call read_succ(success, "totalsucc.bin","unformatted")
	call read_fail(fail, "totalfail.bin","unformatted")

	!If comparing vs a shuffled set, then shuffle.
	if (shuffling) then
		if (printing) print*,"Shuffling the data set."
		call shuffle_cut(success,fail)
	end if
	!Sort the success and fail sets by value in first column.
	if (printing) print*, "Sorting files."
	call heapsort(success)
	call heapsort(fail)

	!Get the core points.
	if (printing) print*,"Getting core points."
	!Set eps in each dimn.
	allocate(eps(size(success,2)))
	!Set loop end st numb pts in box ~ size(success,1)/20
	kend=(size(success,1)/100)/20

  if (find_min) then
    !Set eps to order of Hubble parameter, i.e. the size of quantum fluctuations
    !in the fields.
    eps=(energy_scale**2)/mplanck
    dencrit=1
    if (printing) print*,"Printing core points for eps=",eps
		!Name file.
		write(corename,'(a,i4.4,a)')'corepoints',(k+1),".bin"
		!Open file.
		open(unit=newunit(u),file=corename,form='unformatted')
		do i=1,size(insulatedpts,1)
			write(unit=u), (insulatedpts(i,j),j=1,size(insulatedpts,2))
		end do
		close(unit=u)		
		!Ratio of points in cluster to points in success.
    ratio=real(size(insulatedpts))/real(size(success))
		if (printing) print*,ratio, "Percent of total are core points"
		!Deallocate the insulatedpts array for next loop.
		if(allocated(insulatedpts)) deallocate(insulatedpts)
  else
  	do k=0,kend
  		if (auto) then
  			!Auto set eps to n times avg spatial distance.
  			call set_eps(success,eps,(100*k+1))
  			if (printing) print*, "Epsilon is", eps
  			dencrit=1	!No other points required in eps-ball.
  		else
  			eps=.5_dp
  			dencrit=2	!At least one other point in eps-ball.
  		end if
  		call get_insulatedcorepts(insulatedpts,success,fail,&
  			&euclidean,eps,dencrit)
  		!Print the core points
  		if (printing) print*,"Printing core points for eps=",eps
  		!Name file.
  		write(corename,'(a,i4.4,a)')'corepoints',(k+1),".bin"
  		!Open file.
  		open(unit=newunit(u),file=corename,form='unformatted')
  		do i=1,size(insulatedpts,1)
  			write(unit=u), (insulatedpts(i,j),j=1,size(insulatedpts,2))
  		end do
  		close(unit=u)		
  		!Ratio of points in cluster to points in success.
      ratio=real(size(insulatedpts))/real(size(success))
  		if (printing) print*,ratio, "Percent of total are core points"
  		!Deallocate the insulatedpts array for next loop.
  		if(allocated(insulatedpts)) deallocate(insulatedpts)
  		!Leave if we're not getting any corepoints after 5 iterations.
  		if (k>5 .and. ratio<1E-6) exit
  	end do
  end if

	if(allocated(eps)) deallocate(eps)
	if(allocated(success)) deallocate(success)
	if(allocated(fail)) deallocate(fail)

end program findclusters

