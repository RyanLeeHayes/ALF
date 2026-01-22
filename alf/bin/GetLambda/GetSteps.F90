program GetSteps
! Adapted from msld_readld and msld_readld_header in charmm/source/pert/lambdadyn.F90

  integer,parameter :: chm_int4 = selected_int_kind(9)
  integer,parameter :: chm_real4 = selected_real_kind(6,37)

  ! My i/o variables
  ! https://fortran-lang.org/en/learn/best_practices/file_io/
  logical fileExists
  integer :: argc, a
  character(len=1000), dimension(:), allocatable :: argv
  integer :: iin, iout
  integer :: io

  ! Charmm-like variables
  real(chm_real4), allocatable, dimension(:) :: bielamb, bxlamb2, thetamld
  real(chm_int4), allocatable, dimension(:) :: isitemld
  integer(chm_int4), dimension(20) :: icntrl2
  character(len=80) :: title
  character(len=4) :: hdrr
  real(chm_real4) :: delta4, trajtemp
  integer(chm_int4) :: natom

  ! Other variables
  integer :: i, j, nblock, ifile, nfile, steps
  character(len=80) :: fmt86

  argc = command_argument_count()
  allocate(argv(argc))
  do a = 1, argc
    call get_command_argument(a,argv(a))
  end do

  ! write(io, *) a, b
  inquire(file=argv(1), exist=fileExists)
  if (fileExists) then
    open(newunit=iin, file=argv(1), status="old", form="unformatted", action="read")

    read(iin,iostat=io) hdrr, icntrl2
    if (io==0) then
    nblock = icntrl2(7)
    if (.not.allocated(bielamb)) allocate(bielamb(nblock))
    if (.not.allocated(bxlamb2)) allocate(bxlamb2(nblock))
    if (.not.allocated(thetamld)) allocate(thetamld(nblock))
    if (.not.allocated(isitemld)) allocate(isitemld(nblock))
    write(fmt86,'(A,I7,A)') '(',nblock,'(1X,F8.6))'
    read(iin,iostat=io) delta4
    read(iin,iostat=io) title
    read(iin,iostat=io) natom           ! number of variable biases
    read(iin,iostat=io)                 ! variable_bias parameters
    read(iin,iostat=io) (isitemld(i),i=1,nblock)
    read(iin,iostat=io) trajtemp
    read(iin,iostat=io) (bielamb(i),i=1,nblock)
    if (io==0) then
      nfile=icntrl2(1)
    else
      nfile=0
    endif
    do ifile = 1, nfile
      ! read in lambda values
      read(iin,iostat=io) (bxlamb2(i),i=1,nblock)
      read(iin,iostat=io) (thetamld(i),i=2,nblock)

      if (.not.io==0) then
        exit
      endif
    enddo
    steps = ifile-1

    else
      steps = 0
    endif
    close(iin)
  else
    steps = 0
  endif

  write(*,fmt='(I12)') steps
end program GetSteps
