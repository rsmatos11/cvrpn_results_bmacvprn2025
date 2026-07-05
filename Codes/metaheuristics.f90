module metaheuristics

  use instance, only: instance_type
  use solution, only: solution_type,CreateSolution,DeleteSolution,CopySolution
  use statistics, only: statistics_type,hlenmax
  use neighborhoods, only: Shake,Neighbor
  
  implicit none

  logical, parameter :: debug = .false.
  real(kind=8), parameter :: macheps12 = sqrt( epsilon( 1.0d0 ) )
  
  ! EXTERNAL FUNCTIONS
  real(kind=8), external :: drand 
  
  public :: GVNS
  
contains

  ! *****************************************************************
  ! GVNS - Gerneral Variable Neigborhood Search
  ! *****************************************************************

  subroutine GVNS(stat,inst,sol,seed,otype,ptype,recomp,nrend,norder,qmax,flipallowed,tlim)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in) :: flipallowed,recomp
    integer, intent(in) :: qmax,nrend,norder
    real(kind=8), intent(in) :: tlim
    real(kind=8), intent(inout) :: seed
    character, intent(in) :: otype,ptype
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(inout) :: sol
    type(statistics_type), intent(inout) :: stat

    ! LOCAL SCALARS
    logical :: feasible,improve,tlimreached,recomptmp
    integer :: i,j,k,jnovo,jnovoend,ell,ncsolbar,ncsolhat,nnsolbar,nnsolhat,noimprovecycles,nssolbar,nssolhat,temp
    real(kind=8) :: fsol,fsolbar,tini,tsolbar,tsolhat
    type(solution_type) :: solbar,solhat

    ! LOCAL ARRAYS
    character(len=2) :: ptypetmp
    logical :: singleroute(9)
    integer :: order(9)

    stat%hlen = stat%hlen + 1
    if ( stat%hlen .gt. hlenmax ) then
       write(*,*) 'ERROR: Increase hlenmax in statistics.f90 and re-run.'
       stop
    end if
    stat%vtime(stat%hlen) = 0.0d0
    stat%vsolval(stat%hlen) = sum( sol%length(1:sol%nroutes) )
    
    call CreateSolution(inst,solbar)
    call CreateSolution(inst,solhat)

    ! type of neighborhoods:
    ! 1: 2opt, 2: 3opt, 3: oropt, 4: 1-0Exc, 5: 1-1Exc, 6: 1-2Exc, 7: ReExc, 8: 2optStar, 9: Cross 
    singleroute(1:9) = (/ .true., .true., .true., .false., .false., .false., .false., .false., .false. /)

    if ( .not. ( 1 .le. norder .and. norder .le. 4 ) ) then
       write(*,*) 'In GVNS, norder must be an integer between 1 and 4.'
       stop
    end if

    ! order of neighborhoods
    if ( norder .eq. 1 ) then
       order(1:9) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /)

    else if ( norder .eq. 2 ) then
       order(1:9) = (/ 4, 5, 6, 7, 8, 9, 1, 2, 3 /)

    else if ( norder .eq. 3 ) then
       order(1:9) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /)

    else if ( norder .eq. 4 ) then
       order(1:9) = (/ 4, 5, 6, 7, 8, 9, 1, 2, 3 /)
    end if

    if ( debug ) write(*,*) 'norder = ',norder,' order (initial) of neighborhoods = ',order(1:9)

    noimprovecycles = 0

    tlimreached = .false.

    call cpu_time(tini)

    do while ( .not. tlimreached )

       stat%ncycles = stat%ncycles + 1

       ! -------------------------
       ! Randomize the neighborhood order after a cycle without improvement
       ! -------------------------
       if ( norder .eq. 3 .or. norder .eq. 4 ) then
          if ( noimprovecycles .ge. 1 ) then

             ! Shuffle the neighborhood order
             do i = 9,2,-1
                j = min( 1 + floor( drand(seed) * i ), i )
                temp = order(i)
                order(i) = order(j)
                order(j) = temp
             end do

             if ( debug ) write(*,*) 'New random order of neighborhoods = ',order(1:9)
             noimprovecycles = 0
          end if
       end if

       k = 1
       do while ( k .le. 9 .and. .not. tlimreached )
          call CopySolution(inst,sol,solbar)

          ! Shake qmax times, randomly selecting one of the 9 neighborhood structures
          do ell = 1,qmax
             i = min( 1 + floor( 9 * drand(seed) ), 9 )

             call Shake(stat,inst,solbar,ptype,recomp,nrend,flipallowed,seed,singleroute(order(i)),feasible,order(i))
             if ( debug ) then
                write(*,*) 'Shake with k = ',order(i),' f = ',sum( solbar%length(1:solbar%nroutes) ),' feasible = ',feasible
             end if
          end do
          stat%nshakes = stat%nshakes + 1

          call cpu_time(tsolbar)
          tsolbar = tsolbar - tini
          ncsolbar = stat%ncycles
          nssolbar = stat%nshakes
          nnsolbar = stat%nneighs

          if ( feasible ) then
             if ( ptype .eq. 'N' ) then
                jnovoend = 27
             else
                jnovoend = 9
             end if
                
             jnovo = 1
             do while ( jnovo .le. jnovoend .and. .not. tlimreached )

                if ( ptype .eq. 'N' ) then
                   if ( 1 .le. jnovo .and. jnovo .le. 9 ) then
                      j = jnovo
                      ptypetmp = 'N'
                      recomptmp = .false.
                   else if ( 10 .le. jnovo .and. jnovo .le. 18 ) then
                      j = jnovo - 9
                      ptypetmp = 'G'
                      recomptmp = .true.
                   else ! if ( 19 .le. jnovo .and. jnovo .le. 27 ) then
                      j = jnovo - 18
                      ptypetmp = 'N'
                      recomptmp = .true.
                   end if

                else
                   j = jnovo   
                   ptypetmp = ptype
                   recomptmp = recomp
                end if
                
                call CopySolution(inst,solbar,solhat)

                call Neighbor(stat,inst,solhat,otype,ptypetmp,recomptmp,nrend,flipallowed,seed,singleroute(order(j)), &
                     improve,order(j))

                if ( debug ) then
                   write(*,*) 'Neighborhood with j = ',order(j),' f = ',sum( solhat%length(1:solhat%nroutes) )
                end if
                stat%nneighs = stat%nneighs + 1

                call cpu_time(tsolhat)
                tsolhat = tsolhat - tini
                ncsolhat = stat%ncycles
                nssolhat = stat%nshakes
                nnsolhat = stat%nneighs

                if ( tsolhat .gt. tlim ) then
                   tlimreached = .true.
                end if

                if ( improve ) then
                   if ( 1 .le. jnovo .and. jnovo .le. 9 ) then
                      jnovo = 1
                   else if ( 10 .le. jnovo .and. jnovo .le. 18 ) then
                      jnovo = 10
                   else ! if ( 19 .le. jnovo .and. jnovo .le. 27 ) then
                      jnovo = 19
                   end if
                   
                   call CopySolution(inst,solhat,solbar)
                   tsolbar = tsolhat
                   ncsolbar = ncsolhat
                   nssolbar = nssolhat
                   nnsolbar = nnsolhat
                else
                   jnovo = jnovo + 1
                end if
             end do

             fsol = sum( sol%length(1:sol%nroutes) )
             fsolbar = sum( solbar%length(1:solbar%nroutes) )
             
             if ( fsolbar .le. fsol - macheps12 * max( 1.0d0, abs( fsol ) ) ) then
                if ( debug ) write(*,*) 'An improved solution was obtained'

                k = 1
                call CopySolution(inst,solbar,sol)
                stat%tsol = tsolbar
                stat%ncsol = ncsolbar 
                stat%nssol = nssolbar 
                stat%nnsol = nnsolbar

                stat%hlen = stat%hlen + 1
                if ( stat%hlen .gt. hlenmax ) then
                   write(*,*) 'ERROR: Increase hlenmax in statistics.f90 and re-run.'
                   stop
                end if
                stat%vtime(stat%hlen) = tsolbar
                stat%vsolval(stat%hlen) = sum( solbar%length(1:solbar%nroutes) )
    
                noimprovecycles = 0
             else
                k = k + 1
                if ( k .gt. 9 ) noimprovecycles = noimprovecycles + 1
             end if
          else
             k = k + 1
             if ( k .gt. 9 ) noimprovecycles = noimprovecycles + 1
          end if

       end do
    end do

    call DeleteSolution(solbar)
    call DeleteSolution(solhat)

  end subroutine GVNS
  
end module metaheuristics
