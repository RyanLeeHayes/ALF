program GetLambda
! Adapted from msld_readld and msld_readld_header in charmm/source/pert/lambdadyn.F90

  integer,parameter :: chm_int4 = selected_int_kind(9)
  integer,parameter :: chm_real4 = selected_real_kind(6,37)

  ! My i/o variables
  ! https://fortran-lang.org/en/learn/best_practices/file_io/
  logical fileExists
  integer :: argc, a
  character(len=1000), dimension(:), allocatable :: argv
  integer :: iin, iout

  ! Charmm-like variables
  real(chm_real4), allocatable, dimension(:) :: bielamb, bxlamb2, thetamld
  real(chm_int4), allocatable, dimension(:) :: isitemld
  integer(chm_int4), dimension(20) :: icntrl2
  character(len=80) :: title
  character(len=4) :: hdrr
  real(chm_real4) :: delta4, trajtemp
  integer(chm_int4) :: natom

  ! Other variables
  integer :: i, j, nblock, ifile, nfile
  character(len=80) :: fmt86

  argc = command_argument_count()
  allocate(argv(argc))
  do a = 1, argc
    call get_command_argument(a,argv(a))
  end do

  inquire(file=argv(1), exist=fileExists)
  if (fileExists) then
    open(newunit=iout, file=argv(1), status="replace", form="formatted", action="write")
  else
    open(newunit=iout, file=argv(1), status="new", form="formatted", action="write")
  endif

  ! write(io, *) a, b
  do a = 2, argc
    open(newunit=iin, file=argv(a), status="old", form="unformatted", action="read")

    read(iin) hdrr, icntrl2
    nblock = icntrl2(7)
    if (.not.allocated(bielamb)) allocate(bielamb(nblock))
    if (.not.allocated(bxlamb2)) allocate(bxlamb2(nblock))
    if (.not.allocated(thetamld)) allocate(thetamld(nblock))
    if (.not.allocated(isitemld)) allocate(isitemld(nblock))
    write(fmt86,'(A,I7,A)') '(',nblock,'(1X,F8.6))'
    read(iin) delta4
    read(iin) title
    read(iin) natom           ! number of variable biases
    read(iin)                 ! variable_bias parameters
    read(iin) (isitemld(i),i=1,nblock)
    read(iin) trajtemp
    read(iin) (bielamb(i),i=1,nblock)
    nfile=icntrl2(1)
    do ifile = 1, nfile
      ! read in lambda values
      read(iin) (bxlamb2(i),i=1,nblock)
      read(iin) (thetamld(i),i=2,nblock)

      write(iout,fmt=fmt86) (bxlamb2(i),i=2,nblock)
    enddo

    close(iin)
  end do

  close(iout)
end program GetLambda
