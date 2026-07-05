module blockcoordinatedescent

  use instance, only: instance_type
  use solution, only: solution_type
  use statistics, only: statistics_type
    
  implicit none

  ! PARAMETERS
  logical, parameter :: debug = .false.
  real(kind=8), parameter :: macheps12 = sqrt( epsilon( 1.0d0 ) )

  public :: BCD
  
contains
  
  ! ****************************************************************
  ! ****************************************************************

  subroutine BCD(stat,inst,rnc,route,pathcoo,f)
    
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: rnc
    real(kind=8), intent(out) :: f
    type(instance_type), intent(in) :: inst
    type(statistics_type) :: stat

    ! ARRAY ARGUMENTS
    integer, intent(in) :: route(rnc)
    real(kind=8), intent(inout) :: pathcoo(2,rnc)

    ! LOCAL SCALARS
    integer :: iter,istop

    call BCDb(inst,rnc,route,pathcoo,f,iter,istop)

    if ( .not. ( 0 .le. istop .and. istop .le. 1 ) ) then
       write(*,*) 'There is an error in BCD. istop must be between 0 and 1.'
       stop
    end if

    stat%bcdncalls = stat%bcdncalls + 1
    stat%bcdtotit = stat%bcdtotit + iter
    stat%bcdmaxit = max( stat%bcdmaxit, iter )
    stat%bcdstp(istop) = stat%bcdstp(istop) + 1
    
  end subroutine BCD
  
  ! ****************************************************************
  ! ****************************************************************

  subroutine BCDb(inst,rnc,route,pathcoo,f,iter,istop)
    
    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: rnc
    integer, intent(out) :: iter,istop
    real(kind=8), intent(out) :: f
    type(instance_type), intent(in) :: inst
    
    ! ARRAY ARGUMENTS
    integer, intent(in) :: route(rnc)
    real(kind=8), intent(inout) :: pathcoo(2,rnc)

    ! It is assumed that pathcoo(1:2,1:rnc) is such that
    ! pathcoo(1:2,i) \in P_route(i) for i = 1, ..., rnc.

    ! istop = 0: lack of progress
    ! istop = 1: maximum of iterations
    
    ! PARAMETERS
    integer, parameter :: maxit = 1000
    real(kind=8), parameter :: macheps12 = sqrt( epsilon( 1.0d0 ) )
    
    ! LOCAL SCALARS
    integer :: i,j
    real(kind=8) :: diffnrm,fprev

    ! LOCAL ARRAY
    integer :: istar(rnc)
    real(kind=8) :: diff(2)
    character(len=7) :: lstar
    
    iter = 0
    
    ! Objective at initial guess
    
    f = 0.0d0
    
    diff(1:2) = inst%p(1:2,1,0) - pathcoo(1:2,1)
    diffnrm = norm2( diff(1:2) )
    f = f + diffnrm
    
    do j = 2,rnc
       diff(1:2) = pathcoo(1:2,j) - pathcoo(1:2,j-1)
       diffnrm = norm2( diff(1:2) )
       f = f + diffnrm
    end do

    diff(1:2) = pathcoo(1:2,rnc) - inst%p(1:2,1,0)
    diffnrm = norm2( diff(1:2) )
    f = f + diffnrm

    if ( debug ) write(*,*) 'iter = ',iter,' f = ',f,' route = ',route(1:rnc)

100 continue
    
    ! =======
    ! Iterate
    ! =======
    
    iter = iter + 1
          
    do i = 1,rnc
       call PolMinimize(inst,rnc,route,pathcoo,i,pathcoo(1:2,i),istar(i),lstar)    
       if ( debug ) write(*,*) 'i = ',i,' pathcoo(1:2,i) = ',pathcoo(1:2,i),'istar(i) = ',istar(i),' lstar = ',lstar
    end do
    
    fprev = f
    
    ! ==================
    ! Objective function
    ! ==================
    
    f = 0.0d0
    
    diff(1:2) = inst%p(1:2,1,0) - pathcoo(1:2,1)
    diffnrm = norm2( diff(1:2) )
    f = f + diffnrm
   
    do j = 2,rnc
       diff(1:2) = pathcoo(1:2,j) - pathcoo(1:2,j-1)
       diffnrm = norm2( diff(1:2) )
       f = f + diffnrm
    end do
    
    diff(1:2) = pathcoo(1:2,rnc) - inst%p(1:2,1,0)
    diffnrm = norm2( diff(1:2) ) 
    f = f + diffnrm

    ! =================
    ! Consistency check
    ! =================
    
    if ( f .gt. fprev * ( 1.0d0 + macheps12 ) ) then
       write(*,*) 'There is something wrong!'
       write(*,*) 'In BCD, current f = ',f,' fprev = ',fprev
       !stop
    end if
    
    if ( debug ) write(*,*) 'iter = ',iter,' f = ',f,' route = ',route(1:rnc)
    
    ! ===============
    ! Should we stop?
    ! ===============

    if ( f .ge. fprev * ( 1.0d0 - macheps12 ) ) then
       if (  debug ) write(*,*) 'Stopping criterion satisfied: There was not significant progress in the last iteration.'
       istop = 0
       return
    end if

    if ( iter .gt. maxit ) then
       if ( debug ) write(*,*) 'Stopping criterion satisfied: Too many iterations.'
       istop = 1
       return
    end if

    go to 100

  end subroutine BCDb
  
  ! ****************************************************************
  ! ****************************************************************
  
  subroutine PolMinimize(inst,rnc,route,pathcoo,pos,g,istar,lstar)

    ! This routine computes a point g for visiting the polygon at position pos
    ! in the route. The point g minimizes the sum of Euclidean distances to the
    ! previous point p and next point q in the path. If the segment connecting
    ! p and q intersects the polygon, the intersection point is used. Otherwise,
    ! the point g is selected from the polygon's vertices or along its edges to
    ! minimize the total distance to p and q.

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: pos,rnc
    integer, intent(out) :: istar
    type(instance_type), intent(in) :: inst

    ! ARRAY ARGUMENTS
    integer, intent(in) :: route(rnc)
    real(kind=8), intent(in) :: pathcoo(2,rnc)
    real(kind=8), intent(inout) :: g(2)
    character(len=7), intent(out) :: lstar

    ! LOCAL SCALARS
    logical :: pinpol,qinpol
    integer :: i,next_i,polind
    real(kind=8) :: alphap,alphaq,delta,deltatemp,dnorm2,lambda,pdist,qdist

    ! LOCAL ARRAYS
    real(kind=8) :: d(2),p(2),pdiff(2),q(2),qdiff(2),v(2)

    if ( pos .lt. 1 .or. pos .gt. rnc ) then
       write(*,*) 'There is something wrong!'
       write(*,*) 'pos = ',pos,' must be between 1 and rnc = ',rnc
       stop
    end if

    polind = route(pos)

    if ( inst%np(polind) .eq. 1 ) then
       return
    end if

    if ( pos .eq. 1 ) then
       p(1:2) = inst%p(1:2,1,0)
    else
       p(1:2) = pathcoo(1:2,pos-1)
    end if

    if ( pos .eq. rnc ) then
       q(1:2) = inst%p(1:2,1,0)
    else
       q(1:2) = pathcoo(1:2,pos+1)
    end if

    delta = norm2( p(1:2) - g(1:2) ) + norm2( q(1:2) - g(1:2) )
    lstar = 'initial'
    istar = 0
    
    do i = 1,inst%np(polind)    
       deltatemp = norm2( inst%p(1:2,i,polind) - p(1:2) ) + norm2( inst%p(1:2,i,polind) - q(1:2) )
       if ( deltatemp .lt. delta ) then
          g(1:2) = inst%p(1:2,i,polind)
          delta = deltatemp
          lstar = 'vertex '
          istar = i
       end if

       next_i = mod( i, inst%np(polind) ) + 1
       d(1:2) = inst%p(1:2,next_i,polind) - inst%p(1:2,i,polind)
       dnorm2 = dot_product( d(1:2), d(1:2) )

       if ( dnorm2 .gt. macheps12 ) then

          alphap = dot_product( p(1:2) - inst%p(1:2,i,polind), d(1:2) ) / dnorm2
          alphaq = dot_product( q(1:2) - inst%p(1:2,i,polind), d(1:2) ) / dnorm2

          pdiff(1:2) = p(1:2) - inst%p(1:2,i,polind) - alphap * d(1:2)
          qdiff(1:2) = q(1:2) - inst%p(1:2,i,polind) - alphaq * d(1:2)

          pdist = norm2( pdiff(1:2) )
          qdist = norm2( qdiff(1:2) )

          if ( pdist + qdist .gt. macheps12 ) then

             lambda = ( alphaq * pdist + alphap * qdist ) / ( pdist + qdist )
             v(1:2) = inst%p(1:2,i,polind) + lambda * d(1:2)

             if ( 0.0d0 .le. lambda .and. lambda .le. 1.0d0 ) then
                if ( dot_product( pdiff(1:2), qdiff(1:2) ) .le. 0.0d0 ) then
                   g(1:2) = v(1:2)
                   lstar = 'pqinter'
                   istar = - i
                   return
                end if

                deltatemp = norm2( v(1:2) - p(1:2) ) + norm2( v(1:2) - q(1:2) )
                if ( deltatemp .lt. delta ) then
                   g(1:2) = v(1:2)
                   delta = deltatemp
                   lstar = 'philzr '
                   istar = - i
                end if
             end if

          else
             if ( ( 0.0d0 .lt. alphap .and. alphap .lt. 1.0d0 ) .and. ( 0.0d0 .lt. alphaq .and. alphaq .lt. 1.0d0 ) ) then
                g(1:2) = 0.5d0 * ( p(1:2) + q(1:2) )
                lstar = 'pqcolin'
                istar = - i
                return
             end if
          end if

       end if
    end do

    ! p and q in the interior of the polygon
    call polygon_contains_point_2d(inst%np(polind),inst%p(1:2,1:inst%np(polind),polind),p(1:2),pinpol)
    if ( pinpol ) then
       call polygon_contains_point_2d(inst%np(polind),inst%p(1:2,1:inst%np(polind),polind),q(1:2),qinpol)
       if ( qinpol ) then
          g(1:2) = 0.5d0 * ( p(1:2) + q(1:2) )
          lstar = 'pqinpol'
          istar = - ( inst%np(polind) + 1 )
          return
       end if
    end if

    if ( debug ) then
       write(*,*) 'g     = ',g(1:2)
       write(*,*) 'delta = ',delta
       write(*,*) 'lstar = ',lstar
       write(*,*) 'istar = ',istar
    end if

  end subroutine PolMinimize
    
end module blockcoordinatedescent
