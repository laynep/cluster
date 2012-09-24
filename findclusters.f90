!This is a program that will find clusters in a data set that is loaded into the
!file totalsucc.bin.  The clusters will be referenced against how closely they
!cluster and how close they are to elements of the fail set contained in
!totalfail.bin.

!Options:  
!1.   If we want to print to stdout specify printing=.true.
!3.   If we want to shuffle the success and fail sets and do a run, then
!shuffling=.true.
!4.   If we want to find clusters based off of eps~H, ie the approximate size of
!a quantum fluctuation, then find_min=.true.
!5.   If we want to go from a large eps to a smaller eps, at each stage removing
!the previously found elements from the set, then reduce=.true.

program findclusters
  use fcluster
  use sorters, only : heapsort
  use rng, only : shuffle_cut
  use types, only : dp
  use features, only : newunit
  implicit none

	real(dp), dimension(:,:), allocatable :: success, fail, work
	integer :: i,j,k, kend, check
	integer :: length_s, length_f, width_s, width_f
	real(dp), dimension(:,:), allocatable :: insulatedpts
	real(dp), dimension(:), allocatable :: eps, scaling
	integer :: dencrit, u
	logical :: printing, shuffling, find_min, reduce
	real :: ratio
  real(dp) :: energy_scale, mplanck

	namelist /tablel/ length_s, length_f, width_s, width_f, printing, &
      &shuffling, find_min, reduce
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
	if (printing) print*,"*********************************"
  if (printing) print*,"     Getting core points."
	if (printing) print*,"*********************************"
  !Set eps in each dimn.
	allocate(eps(size(success,2)))
	!Set loop end st numb pts in box ~ size(success,1)/20
	kend=(size(success,1)/100)/20

  !Find clusters where eps is set to the minimum size ~ on order of quant fluct.
  if (find_min) then
    !Set eps to order of Hubble parameter, i.e. the size of quantum fluctuations
    !in the fields.
    eps=(energy_scale**2)/mplanck
    dencrit=1   !No other nearby points necessary.
    !Get the corepoints at the minimum value for epsilong.
		call get_insulatedcorepts(insulatedpts,success,fail,&
  			&euclidean,eps,dencrit)
    call print_corepoints(success, insulatedpts, printing,eps)
    deallocate(insulatedpts)
 else
    !Find clusters in a top-down approach: start with very large clusters, then
    !remove these progressively from the data set until we reach the smallest
    !possible eps~H.

    !Auto set eps to n times avg spatial distance.
    call set_eps(success,eps,(100*kend+1))
    !How much to scale every step by
    allocate(scaling(size(eps)))
    scaling=(eps-(energy_scale**2)/mplanck)/(dble(kend)-1_dp)
 		dencrit=1	!No other points required in eps-ball.
  	do k=kend,1,-1
   		if (printing) print*, "Epsilon is", eps
   		call get_insulatedcorepts(insulatedpts,success,fail,&
  			&euclidean,eps,dencrit)
      call print_corepoints(success, insulatedpts, printing,eps,k)
      !Remove the corepoints from success set.  Removes points that are within
      !1e-10 of a point in the insulatedpts array.
      if (printing) print*, "Taking complement of success set."
      if (allocated(insulatedpts)) then
        call complement(work, success,insulatedpts)
        if(.not. allocated(work)) then
          stop
        end if
        if(allocated(insulatedpts)) deallocate(insulatedpts)
        if(allocated(success)) deallocate(success)
        allocate(success(size(work,1),size(work,2)))
        success=work
        if(allocated(work)) deallocate(work)
      end if
      !Shrink eps.
      eps=eps-scaling
      if (printing) print*,"---------------------"
  	end do
  end if

	if(allocated(eps)) deallocate(eps)
	if(allocated(success)) deallocate(success)
	if(allocated(fail)) deallocate(fail)

end program findclusters

