module instance

  implicit none

  logical, parameter :: debug = .false.
  
  ! nclients is the number of clients
  ! npmax is the maximum number of points or vertices per client
  ! np(k) is the number of points of client k
  ! p(1:2,i,k) Cartesian coordinates of point i of client k
  ! maxlen is the maximum allowed length of a route in a solution
  ! dist(i1,k1,i2,k2) corresponds to the distance from point i1 of client k1 to the point i2 of client k2

  ! By definition client k=0 corresponds to the depot, np(0) = 1, and
  ! p(1:2,1,0) corresponds to the Cartesian coordinates of the depot.
  
  type instance_type
     integer :: nclients,npmax
     real(kind=8) :: maxlen,maxcap,scalefactor
     integer, allocatable :: np(:)
     real(kind=8), allocatable :: p(:,:,:),dist(:,:,:,:),dem(:)
  end type instance_type

  public :: instance_type,CreateInstance,DeleteInstance

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine CreateInstance(inst,filename)

    implicit none

    ! SCALAR ARGUMENTS
    type(instance_type), intent(inout) :: inst

    ! ARRAY ARGUMENTS
    character(len=80), intent(in) :: filename
    
    ! LOCAL SCALARS
    integer :: allocerr,i1,i2,j,k,k1,k2

    inst%maxlen = huge(1.0d0)
        
    open(10,file=trim(adjustl(filename)))

    read(10,*) inst%scalefactor
    read(10,*) inst%maxcap
    read(10,*) inst%nclients
    read(10,*) inst%npmax

    allocate(inst%dem(1:inst%nclients), &
         inst%np(0:inst%nclients), &
         inst%p(2,inst%npmax,0:inst%nclients), &
         inst%dist(1:inst%npmax,0:inst%nclients,1:inst%npmax,0:inst%nclients), &
         stat=allocerr)

    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error.'
       stop
    end if

    ! Demand of each client
    do k = 1,inst%nclients
       read(10,*) inst%dem(k)
    end do
    
    ! Points or vertices of each client
    do k = 0,inst%nclients
       read(10,*) inst%np(k)
       do j = 1,inst%np(k)
          read(10,*) inst%p(1:2,j,k)
       end do
    end do

    ! Distance matrix between all pairs of points    
    do k1 = 0,inst%nclients
       do i1 = 1,inst%np(k1)
          do k2 = k1 + 1,inst%nclients
             do i2 = 1,inst%np(k2)
                inst%dist(i2,k2,i1,k1) = norm2( inst%p(1:2,i2,k2) - inst%p(1:2,i1,k1) )
                inst%dist(i1,k1,i2,k2) = inst%dist(i2,k2,i1,k1)
             end do
          end do
       end do
    end do
    
    close(10)

    if ( debug ) then
       write(*,*) 'maxlen = ',inst%maxlen
       write(*,*) 'maxcap = ',inst%maxcap
       write(*,*) 'nclients = ',inst%nclients
       write(*,*) 'Demand of each client: '
       do k = 1,inst%nclients
          write(*,*) inst%dem(k)
       end do
       write(*,*) 'nclients = ',inst%npmax
       write(*,*) 'Points or vertices of each client: '
       do k = 0,inst%nclients
          write(*,*) 'client: ',k
          write(*,*) 'number of points: ',inst%np(k)
          do j = 1,inst%np(k)
             write(*,*) 'point ',j,' coordinates = ',inst%p(1:2,j,k)
          end do
       end do
       write(*,*) 'Distance matrix between all pairs of points: '
       do k1 = 0,inst%nclients
          do i1 = 1,inst%np(k1)
             do k2 = k1 + 1,inst%nclients
                do i2 = 1,inst%np(k2)
                   write(*,*) 'k1, i1, k2, i2: ',k1,i1,k2,i2,' dist = ',inst%dist(k1,i1,k2,i2)
                   write(*,*) 'k2, i2, k1, i1: ',k2,i2,k1,i1,' dist = ',inst%dist(k2,i2,k1,i1)
                end do
             end do
          end do
       end do
    end if
       
    end subroutine CreateInstance
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine DeleteInstance(inst)

    implicit none

    ! SCALAR ARGUMENTS
    type(instance_type), intent(inout) :: inst

    ! LOCAL SCALARS
    integer :: allocerr

    deallocate(inst%dem,inst%np,inst%p,inst%dist,stat=allocerr)
    
    if ( allocerr .ne. 0 ) then
       write(*,*) 'Deallocation error.'
       stop
    end if

  end subroutine DeleteInstance
  
end module instance
