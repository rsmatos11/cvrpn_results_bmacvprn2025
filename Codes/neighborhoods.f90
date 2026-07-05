module neighborhoods

  use instance, only: instance_type
  use solution, only: solution_type
  use statistics, only: statistics_type
  use dijkstra, only: DijkstraRoute
  use blockcoordinatedescent, only: BCD
  
  implicit none

  logical, parameter :: debug = .false.

  ! EXTERNAL FUNCTIONS
  real(kind=8), external :: drand 
  
  public :: Neighbor,Shake
  
contains
  
  !==============================!
  ! 1. INTRA-ROUTE NEIGHBORHOODS !
  !==============================!
  
  ! *****************************************************************
  ! 2-OPT NEIGHBORHOOD
  ! *****************************************************************
  
  subroutine Move2Opt(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,length,route,pathind,pathcoo)
    
    ! This routine implements the 2-Opt move for the GVRP (Generalized
    ! VRP) problem or the VRPN (VRP with Neighborhoods) problem,
    ! depending on whether the ptype parameter is equal to 'G' or
    ! 'N'. If the parameter otype is equal to 'S' a shake operation is
    ! performed. If the otype parameter is equal to 'B' the best
    ! neighbor of the 2-Opt neighborhood is calculated while if otype
    ! is equal to 'F' the routine ends when the first neighbor of the
    ! 2-Opt neighborhood that improves the current solution is found.
    
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: recomp
    integer, intent(in) :: i,rnc
    character, intent(in) :: otype,ptype
    real(kind=8), intent(inout) :: seed
    real(kind=8), intent(out) :: length
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat

    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients)
    
    ! LOCAL SCALARS
    integer :: j,k
    real(kind=8) :: ltrial
    
    ! LOCAL ARRAYS
    integer :: pindtrial(inst%nclients),rtrial(inst%nclients)
    real(kind=8) :: pcootrial(2,inst%nclients)

    if ( otype .eq. 'S' ) then
       j = min( floor( drand( seed ) * ( ( rnc - 2 ) + 1 ) ), rnc - 2 )
       k = min( j + 2 + floor( drand( seed ) * ( rnc - ( j + 2 ) + 1 ) ), rnc )
              
       call SingleRoute(stat,inst,sol,ptype,recomp,i,j,k,rnc,length,route,pathind,pathcoo)
       if ( length .gt. inst%maxlen ) then
          length = - 1.0d0
       end if
       
    else ! if ( otype .eq. 'B' .or. otype .eq. 'F' ) then
       length = sol%length(i)
       do j = 0,rnc - 2
          do k = j + 2,rnc
             call SingleRoute(stat,inst,sol,ptype,recomp,i,j,k,rnc,ltrial,rtrial,pindtrial,pcootrial)
             if ( ltrial .le. inst%maxlen .and. ltrial .lt. length ) then
                length = ltrial
                route(1:rnc) = rtrial(1:rnc)
                pathind(1:rnc) = pindtrial(1:rnc)
                pathcoo(1:2,1:rnc) = pcootrial(1:2,1:rnc)
                if ( otype .eq. 'F' ) return
             end if
          end do
       end do
    end if
    
  contains
    
    subroutine SingleRoute(stat,inst,sol,ptype,recomp,i,j,k,rnc,ltrial,rtrial,pindtrial,pcootrial)
      
      implicit none
      
      ! SCALAR ARGUMENTS
      logical, intent(in) :: recomp
      integer, intent(in) :: i,j,k,rnc
      character, intent(in) :: ptype
      real(kind=8), intent(out) :: ltrial
      type(instance_type), intent(in) :: inst
      type(solution_type), intent(in) :: sol
      type(statistics_type), intent(inout) :: stat
      
      ! ARRAY ARGUMENTS
      integer, intent(out) :: pindtrial(inst%nclients),rtrial(inst%nclients)
      real(kind=8), intent(out) :: pcootrial(2,inst%nclients)
    
      ! LOCAL SCALARS
      integer :: ell
      
      ! The block from j+1 to k of route i is inverted.
      rtrial(1:rnc) = (/ sol%route(1:j,i), sol%route(k:j+1:-1,i), sol%route(k+1:rnc,i) /)
      
      if ( recomp ) then 
         if ( ptype .eq. 'G' ) then
            call DijkstraRoute(inst,rnc,rtrial,pindtrial,ltrial)
            do ell = 1,rnc
               pcootrial(1:2,ell) = inst%p(1:2,pindtrial(ell),rtrial(ell))
            end do
         else ! if ( ptype .eq. 'N' ) then
            pcootrial(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:j,i), sol%pathcoo(1:2,k:j+1:-1,i), &
                 sol%pathcoo(1:2,k+1:rnc,i) /), (/ 2, rnc /) )
            call BCD(stat,inst,rnc,rtrial,pcootrial,ltrial)
            pindtrial(1:rnc) = - 1
         end if
         
      else ! if ( .not. recomp ) then
         pindtrial(1:rnc) = (/ sol%pathind(1:j,i), sol%pathind(k:j+1:-1,i), sol%pathind(k+1:rnc,i) /)                
         pcootrial(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:j,i), sol%pathcoo(1:2,k:j+1:-1,i), &
              sol%pathcoo(1:2,k+1:rnc,i) /), (/ 2, rnc /) )               
         
         if ( ptype .eq. 'G' ) then
            ltrial = inst%dist(pindtrial(1),rtrial(1),1,0)
            do ell = 2,rnc
               ltrial = ltrial + inst%dist(pindtrial(ell),rtrial(ell),pindtrial(ell-1),rtrial(ell-1))
            end do
            ltrial = ltrial + inst%dist(1,0,pindtrial(rnc),rtrial(rnc))
            
         else ! if ( ptype .eq. 'N' ) then
            ltrial = norm2( pcootrial(1:2,1) - inst%p(1:2,1,0) )
            do ell = 2,rnc
               ltrial = ltrial + norm2( pcootrial(1:2,ell) - pcootrial(1:2,ell-1) )
            end do
            ltrial = ltrial + norm2( inst%p(1:2,1,0) - pcootrial(1:2,rnc) )
         end if
      end if
      
      if ( debug ) write(*,*) 'j = ',j,' k = ',k,' rtrial = ',rtrial(1:rnc),' ltrial = ',ltrial
      
    end subroutine SingleRoute

  end subroutine Move2Opt

  ! *****************************************************************
  ! 3-OPT NEIGHBORHOOD
  ! *****************************************************************

  subroutine Move3Opt(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,length,route,pathind,pathcoo)
   
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: recomp
    integer, intent(in) :: i,rnc
    character, intent(in) :: otype,ptype
    real(kind=8), intent(inout) :: seed
    real(kind=8), intent(out) :: length
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat

    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients)

    ! PARAMETERS
    integer, parameter :: ntrialmax = 100
    
    ! LOCAL SCALARS
    logical :: twocons,threecons
    integer :: case,j,k,ntrial,s
    real(kind=8) :: ltrial
    
    ! LOCAL ARRAYS
    integer :: pindtrial(inst%nclients),rtrial(inst%nclients)
    real(kind=8) :: pcootrial(2,inst%nclients)
    
    if ( otype .eq. 'S' ) then
       ntrial = 0
       
100    continue

       ntrial = ntrial + 1
       
       case = min( 1 + floor( drand( seed ) * 4 ), 4 )
       
       j = min( floor( drand( seed ) * ( ( rnc - 3 ) + 1 ) ), rnc - 3 )
       k = min( j + 1 + floor( drand( seed ) * ( ( rnc - 1 ) - ( j + 1 ) + 1 ) ), rnc - 1 )
       s = min( k + 1 + floor( drand( seed ) * ( rnc - ( k + 1 ) + 1 ) ), rnc )
       
       twocons = j + 1 .eq. k .or. k + 1 .eq. s .or. ( s .eq. rnc .and. j .eq. 0 )       
       threecons = ( j + 1 .eq. k .and. k + 1 .eq. s ) .or. &
            ( j .eq. 0 .and. s .eq. rnc .and. ( j + 1 .eq. k .or. k + 1 .eq. s ) )
       
       if ( ( case .eq. 1 .and. twocons .and. .not. threecons ) .or. &
            ( ( 2 .le. case .and. case .le. 4 ) .and. .not. ( twocons .or. threecons ) ) ) then       
          call SingleRoute(stat,inst,sol,ptype,recomp,case,i,j,k,s,rnc,length,route,pathind,pathcoo)
          if ( length .gt. inst%maxlen ) then
             length = - 1.0d0
          end if
          
       else
          if ( ntrial .lt. ntrialmax ) then
             go to 100
          else
             length = - 1.0d0
             if ( debug ) write(*,*) 'In Move3Opt, with otype equal to S, ntrialmax reached, but no feasible neighbor was found!'
          end if
       end if
       
    else ! if ( otype .eq. 'B' .or. otype .eq. 'F' ) then
       length = sol%length(i)
       do j = 0,rnc - 3
          do k = j + 1,rnc - 1
             do s = k + 1,rnc
                twocons = j + 1 .eq. k .or. k + 1 .eq. s .or. ( s .eq. rnc .and. j .eq. 0 )
                threecons = ( j + 1 .eq. k .and. k + 1 .eq. s ) .or. &
                     ( j .eq. 0 .and. s .eq. rnc .and. ( j + 1 .eq. k .or. k + 1 .eq. s ) )
                
                do case = 1,4
                   if ( ( case .eq. 1 .and. twocons .and. .not. threecons ) .or. &
                        ( ( 2 .le. case .and. case .le. 4 ) .and. .not. ( twocons .or. threecons ) ) ) then
                      call SingleRoute(stat,inst,sol,ptype,recomp,case,i,j,k,s,rnc,ltrial,rtrial,pindtrial,pcootrial)
                      if ( ltrial .le. inst%maxlen .and. ltrial .lt. length ) then
                         length = ltrial
                         route(1:rnc) = rtrial(1:rnc)
                         pathind(1:rnc) = pindtrial(1:rnc)
                         pathcoo(1:2,1:rnc) = pcootrial(1:2,1:rnc)
                         if ( otype .eq. 'F' ) return
                      end if
                   end if
                end do
             end do
          end do
       end do
    end if
    
  contains
    
    subroutine SingleRoute(stat,inst,sol,ptype,recomp,case,i,j,k,s,rnc,ltrial,rtrial,pindtrial,pcootrial)
      
      implicit none
      
      ! SCALAR ARGUMENTS
      logical, intent(in) :: recomp
      integer, intent(in) :: case,i,j,k,s,rnc
      character, intent(in) :: ptype
      real(kind=8), intent(out) :: ltrial
      type(instance_type), intent(in) :: inst
      type(solution_type), intent(in) :: sol
      type(statistics_type), intent(inout) :: stat
      
      ! ARRAY ARGUMENTS
      integer, intent(out) :: pindtrial(inst%nclients),rtrial(inst%nclients)
      real(kind=8), intent(out) :: pcootrial(2,inst%nclients)
      
      ! LOCAL SCALARS
      integer :: ell
      
      if ( case .eq. 1 ) then
         rtrial(1:rnc) = (/ sol%route(1:j,i), sol%route(k+1:s,i), sol%route(j+1:k,i), sol%route(s+1:rnc,i) /)
      else if ( case .eq. 2 ) then
         rtrial(1:rnc) = (/ sol%route(1:j,i), sol%route(k+1:s,i), sol%route(k:j+1:-1,i), sol%route(s+1:rnc,i) /)
      else if ( case .eq. 3 ) then
         rtrial(1:rnc) = (/ sol%route(1:j,i), sol%route(k:j+1:-1,i), sol%route(s:k+1:-1,i), sol%route(s+1:rnc,i) /)
      else ! if ( case .eq. 4 ) then
         rtrial(1:rnc) = (/ sol%route(1:j,i), sol%route(s:k+1:-1,i), sol%route(j+1:k,i), sol%route(s+1:rnc,i) /)
      end if
      
      if ( recomp ) then 
         if ( ptype .eq. 'G' ) then
            call DijkstraRoute(inst,rnc,rtrial,pindtrial,ltrial)
            do ell = 1,rnc
               pcootrial(1:2,ell) = inst%p(1:2,pindtrial(ell),rtrial(ell))
            end do
         else ! if ( type .eq. 'N' ) then
            if ( case .eq. 1 ) then
               pcootrial(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:j,i), sol%pathcoo(1:2,k+1:s,i), &
                    sol%pathcoo(1:2,j+1:k,i), sol%pathcoo(1:2,s+1:rnc,i) /), (/ 2, rnc /) )               
            else if ( case .eq. 2 ) then
               pcootrial(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:j,i), sol%pathcoo(1:2,k+1:s,i), &
                    sol%pathcoo(1:2,k:j+1:-1,i), sol%pathcoo(1:2,s+1:rnc,i) /), (/ 2, rnc /) )
            else if ( case .eq. 3 ) then
               pcootrial(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:j,i), sol%pathcoo(1:2,k:j+1:-1,i), &
                    sol%pathcoo(1:2,s:k+1:-1,i), sol%pathcoo(1:2,s+1:rnc,i) /), (/ 2, rnc /) )               
            else ! if ( case .eq. 4 ) then
               pcootrial(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:j,i), sol%pathcoo(1:2,s:k+1:-1,i), &
                    sol%pathcoo(1:2,j+1:k,i), sol%pathcoo(1:2,s+1:rnc,i) /), (/ 2, rnc /) )               
            end if
            call BCD(stat,inst,rnc,rtrial,pcootrial,ltrial)
            pindtrial(1:rnc) = - 1
         end if
         
      else ! if ( .not. recomp ) then
         if ( case .eq. 1 ) then
            pindtrial(1:rnc) = (/ sol%pathind(1:j,i), sol%pathind(k+1:s,i), sol%pathind(j+1:k,i), &
                 sol%pathind(s+1:rnc,i) /)                
            pcootrial(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:j,i), sol%pathcoo(1:2,k+1:s,i), &
                 sol%pathcoo(1:2,j+1:k,i), sol%pathcoo(1:2,s+1:rnc,i) /), (/ 2, rnc /) )        
         else if ( case .eq. 2 ) then
            pindtrial(1:rnc) = (/ sol%pathind(1:j,i), sol%pathind(k+1:s,i), sol%pathind(k:j+1:-1,i), &
                 sol%pathind(s+1:rnc,i) /)                
            pcootrial(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:j,i), sol%pathcoo(1:2,k+1:s,i), &
                 sol%pathcoo(1:2,k:j+1:-1,i), sol%pathcoo(1:2,s+1:rnc,i) /), (/ 2, rnc /) )        
         else if ( case .eq. 3 ) then
            pindtrial(1:rnc) = (/ sol%pathind(1:j,i), sol%pathind(k:j+1:-1,i), sol%pathind(s:k+1:-1,i), &
                 sol%pathind(s+1:rnc,i) /)                
            pcootrial(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:j,i), sol%pathcoo(1:2,k:j+1:-1,i), &
                 sol%pathcoo(1:2,s:k+1:-1,i), sol%pathcoo(1:2,s+1:rnc,i) /), (/ 2, rnc /) )        
         else ! if ( case .eq. 4 ) then
            pindtrial(1:rnc) = (/ sol%pathind(1:j,i), sol%pathind(s:k+1:-1,i), sol%pathind(j+1:k,i), &
                 sol%pathind(s+1:rnc,i) /)                
            pcootrial(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:j,i), sol%pathcoo(1:2,s:k+1:-1,i), &
                 sol%pathcoo(1:2,j+1:k,i), sol%pathcoo(1:2,s+1:rnc,i) /), (/ 2, rnc /) )        
         end if
         
         if ( ptype .eq. 'G' ) then
            ltrial = inst%dist(pindtrial(1),rtrial(1),1,0)
            do ell = 2,rnc
               ltrial = ltrial + inst%dist(pindtrial(ell),rtrial(ell),pindtrial(ell-1),rtrial(ell-1))
            end do
            ltrial = ltrial + inst%dist(1,0,pindtrial(rnc),rtrial(rnc))
            
         else ! if ( ptype .eq. 'N' ) then
            ltrial = norm2( pcootrial(1:2,1) - inst%p(1:2,1,0) )
            do ell = 2,rnc
               ltrial = ltrial + norm2( pcootrial(1:2,ell) - pcootrial(1:2,ell-1) )
            end do
            ltrial = ltrial + norm2( inst%p(1:2,1,0) - pcootrial(1:2,rnc) )
         end if
      end if
      
      if ( debug ) then
         write(*,*) 'case =  ',case,' j = ',j,' k = ',k,' s = ',s,' rtrial = ',rtrial(1:rnc),' ltrial = ',ltrial            
      end if
      
    end subroutine SingleRoute
    
  end subroutine Move3Opt

  ! *****************************************************************
  ! OR-OPT NEIGHBORHOOD
  ! *****************************************************************

  subroutine MoveOrOpt(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,length,route,pathind,pathcoo,nrini,nrend,flipallowed)
     
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: flipallowed,recomp
    integer, intent(in) :: i,nrini,nrend,rnc
    character, intent(in) :: otype,ptype
    real(kind=8), intent(inout) :: seed
    real(kind=8), intent(out) :: length
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat

    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients)
    
    ! LOCAL SCALARS
    integer :: j,k,nr,nrendt,op1,op1max
    real(kind=8) :: ltrial
        
    ! LOCAL ARRAYS
    logical :: flip(2)
    integer :: pindtrial(inst%nclients),rtrial(inst%nclients)
    real(kind=8) :: pcootrial(2,inst%nclients)

    flip(1:2) = (/ .false., .true. /)
    
    if ( otype .eq. 'S' ) then

       nrendt = min( nrend, rnc - 1 )
       nr = min( nrini + floor( drand(seed) * ( nrendt - nrini + 1 ) ), nrendt )
       
       j = min( 1 + floor( drand( seed ) * ( rnc - nr + 1 ) ), rnc - nr + 1 )
       k = min( 1 + floor( drand( seed ) * ( rnc - nr ) ), rnc - nr )
       if ( k .ge. j ) k = k + 1       
       if ( .not. flipallowed ) then
          op1 = 1
       else
          op1 = min( 1 + floor( drand( seed ) * min( nr, 2 ) ), min( nr, 2 ) )
       end if
       if ( debug ) write(*,*) 'nr = ',nr,' j = ',j,' k = ',k,' op1 =',op1,' flip(op1) = ',flip(op1)
       
       call SingleRoute(stat,inst,sol,ptype,recomp,i,j,k,nr,flip(op1),rnc,length,route,pathind,pathcoo)
       
       if ( debug ) write(*,*) 'route = ',route(1:rnc),' length = ',length
       
       if ( length .gt. inst%maxlen ) then
          length = -1.0d0
       end if
 
    else ! if ( otype .eq. 'B' .or. otype .eq. 'F' ) then
       length = sol%length(i)      
       ! In or-opt we move nr consecutive clients from route i to all other positions on route i
       do nr = nrini,nrend       
          do j = 1,rnc - ( nr - 1 )
             do k = 1,rnc - ( nr - 1 )             
                if ( k .ne. j ) then
                   if ( .not. flipallowed ) then
                      op1max = 1
                   else
                      op1max = 2
                   end if
                   
                   do op1 = 1,min( nr, op1max )
                      call SingleRoute(stat,inst,sol,ptype,recomp,i,j,k,nr,flip(op1),rnc,ltrial,rtrial,pindtrial,pcootrial)
                      
                      if ( debug ) write(*,*) 'nr = ',nr,' j = ',j,' k = ',k,' flip = ',flip(op1),' rtrial = ',rtrial(1:rnc), &
                           ' ltrial = ',ltrial
                      
                      if ( ltrial .le. inst%maxlen .and. ltrial .lt. length ) then
                         length = ltrial
                         route(1:rnc) = rtrial(1:rnc)
                         pathind(1:rnc) = pindtrial(1:rnc)
                         pathcoo(1:2,1:rnc) = pcootrial(1:2,1:rnc)
                         if ( otype .eq. 'F' ) return
                      end if
                   end do
                end if
             end do
          end do
       end do
    end if
    
  contains
    
    subroutine SingleRoute(stat,inst,sol,ptype,recomp,i,j,k,nr,flip,rnc,ltrial,rtrial,pindtrial,pcootrial)
      
      implicit none
      
      ! SCALAR ARGUMENTS
      logical, intent(in) :: flip,recomp
      integer, intent(in) :: i,j,k,nr,rnc
      character, intent(in) :: ptype
      real(kind=8), intent(out) :: ltrial
      type(instance_type), intent(in) :: inst
      type(solution_type), intent(in) :: sol
      type(statistics_type), intent(inout) :: stat

      ! ARRAY ARGUMENTS
      integer, intent(out) :: pindtrial(inst%nclients),rtrial(inst%nclients)
      real(kind=8), intent(out) :: pcootrial(2,inst%nclients)
    
      ! LOCAL SCALARS
      integer :: ell

      ! Remove clients from j to j + nr - 1
      rtrial(1:rnc-nr) = (/ sol%route(1:j-1,i), sol%route(j+nr:rnc,i) /)
      ! Reallocate starting at position k
      if ( .not. flip ) then
         rtrial(1:rnc) = (/ rtrial(1:k-1), sol%route(j:j+nr-1,i), rtrial(k:rnc-nr) /)
      else
         rtrial(1:rnc) = (/ rtrial(1:k-1), sol%route(j+nr-1:j:-1,i), rtrial(k:rnc-nr) /)
      end if
         
      if ( recomp ) then 
         if ( ptype .eq. 'G' ) then
            call DijkstraRoute(inst,rnc,rtrial,pindtrial,ltrial)
            do ell = 1,rnc
               pcootrial(1:2,ell) = inst%p(1:2,pindtrial(ell),rtrial(ell))
            end do
         else ! if ( type .eq. 'N' ) then
            pcootrial(1:2,1:rnc-nr) = reshape( (/ sol%pathcoo(1:2,1:j-1,i), sol%pathcoo(1:2,j+nr:rnc,i) /), &
                 (/ 2, rnc-nr /) )
            if( .not. flip ) then
               pcootrial(1:2,1:rnc) = reshape( (/ pcootrial(1:2,1:k-1), sol%pathcoo(1:2,j:j+nr-1,i), &
                    pcootrial(1:2,k:rnc-nr) /), (/ 2, rnc /) )
            else
               pcootrial(1:2,1:rnc) = reshape( (/ pcootrial(1:2,1:k-1), sol%pathcoo(1:2,j+nr-1:j:-1,i), &
                    pcootrial(1:2,k:rnc-nr) /), (/ 2, rnc /) )
            end if
            call BCD(stat,inst,rnc,rtrial,pcootrial,ltrial)
            pindtrial(1:rnc) = - 1
         end if

      else ! if ( .not. recomp ) then
         pindtrial(1:rnc-nr) = (/ sol%pathind(1:j-1,i), sol%pathind(j+nr:rnc,i) /)
         if ( .not. flip ) then
            pindtrial(1:rnc) = (/ pindtrial(1:k-1), sol%pathind(j:j+nr-1,i), pindtrial(k:rnc-nr) /)
         else
            pindtrial(1:rnc) = (/ pindtrial(1:k-1), sol%pathind(j+nr-1:j:-1,i), pindtrial(k:rnc-nr) /)
         end if            
         pcootrial(1:2,1:rnc-nr) = reshape( (/ sol%pathcoo(1:2,1:j-1,i), sol%pathcoo(1:2,j+nr:rnc,i) /), (/ 2, rnc-nr /) )
         if ( .not. flip ) then
            pcootrial(1:2,1:rnc) = reshape( (/ pcootrial(1:2,1:k-1), sol%pathcoo(1:2,j:j+nr-1,i), pcootrial(1:2,k:rnc-nr) /), &
                 (/ 2, rnc /) )
         else
            pcootrial(1:2,1:rnc) = reshape( (/ pcootrial(1:2,1:k-1), sol%pathcoo(1:2,j+nr-1:j:-1,i), pcootrial(1:2,k:rnc-nr) /), &
                 (/ 2, rnc /) )
         end if
         
         if ( ptype .eq. 'G' ) then
            ltrial = inst%dist(pindtrial(1),rtrial(1),1,0)
            do ell = 2,rnc
               ltrial = ltrial + inst%dist(pindtrial(ell),rtrial(ell),pindtrial(ell-1),rtrial(ell-1))
            end do
            ltrial = ltrial + inst%dist(1,0,pindtrial(rnc),rtrial(rnc))
            
         else ! if ( ptype .eq. 'N' ) then
            ltrial = norm2( pcootrial(1:2,1) - inst%p(1:2,1,0) )
            do ell = 2,rnc
               ltrial = ltrial + norm2( pcootrial(1:2,ell) - pcootrial(1:2,ell-1) )
            end do
            ltrial = ltrial + norm2( inst%p(1:2,1,0) - pcootrial(1:2,rnc) )
         end if
      end if
             
    end subroutine SingleRoute

  end subroutine MoveOrOpt

  !==============================!
  ! 2. INTER-ROUTE NEIGHBORHOODS !
  !==============================!
    
  subroutine RemoveClients(stat,inst,sol,ptype,recomp,i,jini,jend,rnc,length,route,pathind,pathcoo)
      
    ! Remove clients from jini to jend of route i.
    
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: recomp
    integer, intent(in) :: i,jini,jend
    character, intent(in) :: ptype
    integer, intent(inout) :: rnc
    real(kind=8), intent(out) :: length
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat
    
    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients)
    
    if ( rnc .eq. 0 ) then
       if ( debug ) write(*,*) 'jini = ',jini,' jend = ',jend,' the route is empty'
       length = 0.0d0
       return
    end if
    
    call SubstituteClients(stat,inst,sol,ptype,recomp,i,jini,jini-1,.false.,i,jini,jend,rnc,length,route,pathind,pathcoo)
    
  end subroutine RemoveClients
  
  ! *****************************************************************
  ! *****************************************************************
  
  subroutine InsertClients(stat,inst,sol,ptype,recomp,isource,jini,jend,flip,idest,k,rnc,length,route,pathind,pathcoo)
    
    ! Insert in position k of route idest clients from jini to jend
    ! of route isource.
    
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: flip,recomp
    integer, intent(in) :: idest,isource,jini,jend,k
    integer, intent(inout) :: rnc
    character, intent(in) :: ptype
    real(kind=8), intent(out) :: length
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat
    
    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients)
    
    call SubstituteClients(stat,inst,sol,ptype,recomp,isource,jini,jend,flip,idest,k,k-1,rnc,length,route,pathind,pathcoo)
    
  end subroutine InsertClients
  
  ! *****************************************************************
  ! *****************************************************************
  
  subroutine SubstituteClients(stat,inst,sol,ptype,recomp,isource,jini,jend,flip,idest,kini,kend,rnc,length,route,pathind,pathcoo)
    
    ! Remove clients from kini to kend of route idest and then insert
    ! in position kini clients from jini to jend of route isource. If
    ! flip = true, then the segment to be inserted is reversed.
    
    ! If jend = jini - 1, this routine does only the removal part.
    ! If kend = kini - 1, this routine does only the insert part.
    
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: flip,recomp
    integer, intent(in) :: isource,idest,jini,jend,kini,kend
    integer, intent(inout) :: rnc
    character, intent(in) :: ptype
    real(kind=8), intent(out) :: length
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat
    
    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients)
    
    ! LOCAL SCALARS
    integer :: ell,rncorig
    
    rncorig = rnc
    rnc = rnc - ( kend - kini + 1 ) + ( jend - jini + 1 )
    if ( rnc .eq. 0 ) then
       if ( debug ) write(*,*) 'the route ',idest,' is empty.'
       length = 0.0d0
       return
    end if
    
    ! We remove clients from kini to kend from route idest and insert clients from jini to jend from route isource in its place.
    if ( .not. flip ) then
       route(1:rnc) = (/ sol%route(1:kini-1,idest), sol%route(jini:jend,isource),    sol%route(kend+1:rncorig,idest) /)
    else
       route(1:rnc) = (/ sol%route(1:kini-1,idest), sol%route(jend:jini:-1,isource), sol%route(kend+1:rncorig,idest) /)
    end if
    
    if ( recomp ) then 
       if ( ptype .eq. 'G' ) then
          call DijkstraRoute(inst,rnc,route,pathind,length)
          do ell = 1,rnc
             pathcoo(1:2,ell) = inst%p(1:2,pathind(ell),route(ell))
          end do
       else ! if ( type .eq. 'N' ) then
          if ( .not. flip ) then
             pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini-1,idest), sol%pathcoo(1:2,jini:jend,   isource), &
                  sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )
          else
             pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini-1,idest), sol%pathcoo(1:2,jend:jini:-1,isource), &
                  sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )
          end if
          call BCD(stat,inst,rnc,route,pathcoo,length)
          pathind(1:rnc) = - 1
       end if
       
    else ! if ( .not. recomp ) then
       if ( .not. flip ) then
          pathind(1:rnc) = (/ sol%pathind(1:kini-1,idest), sol%pathind(jini:jend,isource),    sol%pathind(kend+1:rncorig,idest) /)
          pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini-1,idest), sol%pathcoo(1:2,jini:jend,   isource), &
               sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )     
       else
          pathind(1:rnc) = (/ sol%pathind(1:kini-1,idest), sol%pathind(jend:jini:-1,isource), sol%pathind(kend+1:rncorig,idest) /)
          pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini-1,idest), sol%pathcoo(1:2,jend:jini:-1,isource), &
               sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )     
       end if
       
       if ( ptype .eq. 'G' ) then
          length = inst%dist(pathind(1),route(1),1,0)
          do ell = 2,rnc
             length = length + inst%dist(pathind(ell),route(ell),pathind(ell-1),route(ell-1))
          end do
          length = length + inst%dist(1,0,pathind(rnc),route(rnc))
          
       else ! if ( ptype .eq. 'N' ) then
          length = norm2( pathcoo(1:2,1) - inst%p(1:2,1,0) )
          do ell = 2,rnc
             length = length + norm2( pathcoo(1:2,ell) - pathcoo(1:2,ell-1) )
          end do
          length = length + norm2( inst%p(1:2,1,0) - pathcoo(1:2,rnc) )
       end if
       
    end if
    
  end subroutine SubstituteClients
  
  ! *****************************************************************
  ! 1-0 EXCHANGE NEIGHBORHOOD
  ! *****************************************************************

  subroutine Move10Exc(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
       i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2)
    
    ! This routine considers neighbors constructed by removing every
    ! segment of size 1 in routes i and inserting it in every possible
    ! position in route i2.
    
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: recomp
    integer, intent(in) :: i,i2
    integer, intent(inout) :: rnc,rnc2
    character, intent(in) :: otype,ptype
    real(kind=8), intent(inout) :: seed
    real(kind=8), intent(out) :: aggdem,aggdem2,length,length2
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat
    
    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients),pathind2(inst%nclients), &
         route2(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients),pathcoo2(2,inst%nclients)

    call MoveReExcPlus(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
       i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2,1,1,.false.)
    
  end subroutine Move10Exc

  ! *****************************************************************
  ! 1-1 EXCHANGE NEIGHBORHOOD
  ! *****************************************************************

  subroutine Move11Exc(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
       i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2)
    
    ! This routine considers neighbors constructed by interchanging
    ! every possible pair of segments of size 1 in routes i and i2.
    
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: recomp
    integer, intent(in) :: i,i2
    integer, intent(inout) :: rnc,rnc2
    character, intent(in) :: otype,ptype
    real(kind=8), intent(inout) :: seed
    real(kind=8), intent(out) :: aggdem,aggdem2,length,length2
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat
    
    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients),pathind2(inst%nclients), &
         route2(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients),pathcoo2(2,inst%nclients)

    ! PARAMETERS
    integer, parameter :: ncombsize = 1
    integer, parameter :: combsize(2,ncombsize) = reshape( (/ 1, 1 /), shape( combsize ) )
    
    call MoveCrossExc(stat,inst,sol,otype,ptype,recomp,seed,ncombsize,combsize,0,0,i,rnc,length,aggdem,route,pathind,pathcoo, &
       0,0,i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2)
    
  end subroutine Move11Exc

  ! *****************************************************************
  ! 1-2 EXCHANGE NEIGHBORHOOD
  ! *****************************************************************

  subroutine Move12Exc(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
       i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2,flipallowed)

    ! This routine considers neighbors constructed by interchanging
    ! every possible pair of segments of size 2 and 1 in routes i and
    ! i2, respectively.
    
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: flipallowed,recomp
    integer, intent(in) :: i,i2
    integer, intent(inout) :: rnc,rnc2
    character, intent(in) :: otype,ptype
    real(kind=8), intent(inout) :: seed
    real(kind=8), intent(out) :: aggdem,aggdem2,length,length2
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat
    
    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients),pathind2(inst%nclients), &
         route2(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients),pathcoo2(2,inst%nclients)

    ! PARAMETERS
    integer, parameter :: ncombsize = 1
    integer, parameter :: combsize(2,ncombsize) = reshape( (/ 2, 1 /), shape( combsize ) )

    ! LOCAL SCALARS
    integer :: fmax
    
    if ( .not. flipallowed ) then
       fmax = 0
    else
       fmax = 1
    end if
    
    call MoveCrossExc(stat,inst,sol,otype,ptype,recomp,seed,ncombsize,combsize,0,fmax,i,rnc,length,aggdem,route,pathind,pathcoo, &
       0,0,i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2)
       
  end subroutine Move12Exc

  ! *****************************************************************
  ! CROSS-EXCHANGE NEIGHBORHOOD
  ! *****************************************************************

  subroutine MoveCrossExc(stat,inst,sol,otype,ptype,recomp,seed,ncombsize,combsize,fmin,fmax,i,rnc,length,aggdem,route,pathind, &
       pathcoo,fmin2,fmax2,i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2)

    ! This routine considers neighbors constructed by interchanging
    ! segments of varying size in routes i and i2.

    ! The size combinations are given in array combsize of dimension
    ! 2 x ncombsize. Each column represents a combination of sizes.
    
    ! If fmin = fmax = 0, then the segment from route i is inserted in
    ! route i2 without flipping only.

    ! If fmin = fmax = 1, then the segment from route i is inserted in
    ! route i2 flipped only.
    
    ! If fmin = 0 and fmax = 1, then both insertions are considered,
    ! unless in the case of segments of unitary size.
    
    ! The same things follows for the segment from route i2 that is
    ! inserted in route i and parameters fmin2, fmax2.
    
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: recomp
    integer, intent(in) :: i,i2,fmin,fmin2,fmax,fmax2,ncombsize
    integer, intent(inout) :: rnc,rnc2
    character, intent(in) :: otype,ptype
    real(kind=8), intent(inout) :: seed
    real(kind=8), intent(out) :: aggdem,aggdem2,length,length2
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat
    
    ! ARRAY ARGUMENTS
    integer, intent(in) :: combsize(2,ncombsize)
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients),pathind2(inst%nclients), &
         route2(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients),pathcoo2(2,inst%nclients)

    ! PARAMETERS
    !integer, parameter :: ncombsizemax = 100
    
    ! LOCAL SCALARS
    logical :: flip,flip2
    integer :: f,g,j,csind,nfc,rnctrial,rnctrial2,s,size,size2
    real(kind=8) :: ltrial,ltrial2,adtrial,adtrial2,tmp1,tmp2
    
    ! LOCAL ARRAYS
    integer :: pindtrial(inst%nclients),rtrial(inst%nclients),pindtrial2(inst%nclients), &
         rtrial2(inst%nclients),fc(ncombsize)!,combsize(2,ncombsizemax)
    real(kind=8) :: pcootrial(2,inst%nclients),pcootrial2(2,inst%nclients)
    
    if ( any( combsize(1:2,1:ncombsize) .lt. 1 ) ) then
       write(*,*) 'In MoveCrossExc: Entries of matrix combisze must be greater than or equal to 1!'
       stop
    end if
    
    if ( otype .eq. 'S' ) then
       
       nfc = count( combsize(1,1:ncombsize) .le. rnc .and. combsize(2,1:ncombsize) .le. rnc2 )
       
       if ( nfc .eq. 0 ) then
          aggdem  = - 1.0d0
          aggdem2 = - 1.0d0
          length  = - 1.0d0
          length2 = - 1.0d0
          if ( debug ) write(*,*) 'In MoveCrossExc, with otype=S: Matrix combsize contains no feasible combinations!'

       else
          fc(1:nfc) = pack( (/ ( i, i=1,ncombsize ) /), combsize(1,1:ncombsize) .le. rnc .and. combsize(2,1:ncombsize) .le. rnc2 )
       
          csind = min( 1 + floor( drand(seed) * ( nfc + 1 ) ), nfc )
          size  = combsize(1,fc(csind))
          size2 = combsize(2,fc(csind))

          j = min( floor( drand( seed ) * ( rnc - size + 1 ) ), rnc - size )
          f = fmin + floor( drand( seed ) * ( fmax - fmin + 1 ) )
          
          s = min( floor( drand( seed ) * ( rnc2 - size2 + 1 ) ), rnc2 - size2 )
          g = fmin2 + floor( drand( seed ) * ( fmax2 - fmin2 + 1 ) )
          
          if ( debug ) write(*,*) 'csind = ',csind,' i, j, size, f, i2, s, size2, g = ',i,j,size,f,i2,s,size2,g
          
          tmp1 = sum(inst%dem(sol%route(j+1:j+size,i)))
          tmp2 = sum(inst%dem(sol%route(s+1:s+size2,i2)))
          aggdem  = sol%aggdem(i)  - tmp1 + tmp2
          aggdem2 = sol%aggdem(i2) + tmp1 - tmp2
             
          if ( aggdem .gt. inst%maxcap .or. aggdem2 .gt. inst%maxcap ) then
             aggdem  = - 1.0d0
             aggdem2 = - 1.0d0
             length  = - 1.0d0
             length2 = - 1.0d0
          
          else
             flip2 = g .ne. 0
             call SubstituteClients(stat,inst,sol,ptype,recomp,i2,s+1,s+size2,flip2,i,j+1,j+size,rnc,length,route,pathind,pathcoo)
             if ( debug ) write(*,*) 'route = ',route(1:rnc),' length = ',length,' aggdem = ',aggdem
             
             flip = f .ne. 0
             call SubstituteClients(stat,inst,sol,ptype,recomp,i,j+1,j+size,flip,i2,s+1,s+size2,rnc2,length2,route2,pathind2, &
                  pathcoo2)
             if ( debug ) write(*,*) 'route2 = ',route2(1:rnc2),' length2 = ',length2,' aggdem2 = ',aggdem2
             
             if ( length .gt. inst%maxlen .or. length2 .gt. inst%maxlen ) then
                length  = - 1.0d0
                length2 = - 1.0d0
             end if
          end if
       end if
          
    else !if ( otype .eq. 'B' .or. otype .eq. 'F' ) then
       ! The following three loops in j, k, and f generate all
       ! possible segments starting at client j+1 and ending at client
       ! k with size = k - j between tmin and tmax. When the size is 0
       ! or 1, it makes no sense to think in flipping the segment or
       ! not when inserting it in route i2. So, in this case f assumes
       ! the value 0, meaning that there will be no flip. When the
       ! size is strictly larger than 1, f assumes values 0 and 1 (0
       ! means no flip and 1 means flip).
       
       length = sol%length(i)
       length2 = sol%length(i2)

       do csind = 1,ncombsize
          size  = combsize(1,csind)
          size2 = combsize(2,csind)

          if ( size .le. rnc .and. size2 .le. rnc2 .and. ( size .ne. size2 .or. i2 .gt. i ) ) then
             do j = 0,sol%nclients(i) - size
                tmp1 = sum(inst%dem(sol%route(j+1:j+size,i)))
                ! The size of the segment to be removed from route i is size.
                ! If size <= 1, then we want flip = .false.
                ! If size  > 1, then we want flip \in { .false., .true. }
                do f = max( 0, fmin ),min( fmax, min( max( 0, size - 1 ), 1 ) )
                   flip = f .ne. 0
                   
                   do s = 0,sol%nclients(i2) - size2
                      tmp2 = sum(inst%dem(sol%route(s+1:s+size2,i2)))
                      do g = max( 0, fmin2 ),min( fmax2, min( max( 0, size2 - 1 ), 1 ) )
                         flip2 = g .ne. 0
                         
                         adtrial  = sol%aggdem(i)  - tmp1 + tmp2
                         adtrial2 = sol%aggdem(i2) + tmp1 - tmp2
                         
                         if ( adtrial .le. inst%maxcap .and. adtrial2 .le. inst%maxcap ) then
                            rnctrial = sol%nclients(i)
                            call SubstituteClients(stat,inst,sol,ptype,recomp,i2,s+1,s+size2,flip2,i,j+1,j+size,rnctrial,ltrial, &
                                 rtrial,pindtrial,pcootrial)
                            
                            if ( debug ) write(*,*) 'j = ', j,' size = ',size,' f = ',f,' s = ',s,' size2 = ',size2,' g = ',g, &
                                 ' rtrial = ',rtrial(1:rnctrial),' ltrial = ',ltrial,' adtrial = ',adtrial
                            
                            rnctrial2 = sol%nclients(i2)
                            call SubstituteClients(stat,inst,sol,ptype,recomp,i,j+1,j+size,flip,i2,s+1,s+size2,rnctrial2,ltrial2, &
                                 rtrial2,pindtrial2,pcootrial2)
                            
                            if ( debug ) write(*,*) 'j = ', j,' size = ',size,' f = ',f,' s = ',s,' size2 = ',size2,' g = ',g, &
                                 ' rtrial2 = ',rtrial2(1:rnctrial2),' ltrial2 = ',ltrial2,' adtrial2 = ',adtrial2
                            
                            if ( ltrial .le. inst%maxlen .and. ltrial2 .le. inst%maxlen ) then
                               if ( ltrial + ltrial2 .lt. length + length2 ) then
                                  length = ltrial
                                  aggdem = adtrial
                                  rnc  = rnctrial
                                  route(1:rnctrial) = rtrial(1:rnctrial)
                                  pathind(1:rnctrial) = pindtrial(1:rnctrial)
                                  pathcoo(1:2,1:rnctrial) = pcootrial(1:2,1:rnctrial)
                                  
                                  length2 = ltrial2
                                  aggdem2 = adtrial2
                                  rnc2 = rnctrial2
                                  route2(1:rnctrial2) = rtrial2(1:rnctrial2)
                                  pathind2(1:rnctrial2) = pindtrial2(1:rnctrial2)
                                  pathcoo2(1:2,1:rnctrial2) = pcootrial2(1:2,1:rnctrial2)
                                  if ( otype .eq. 'F' ) return
                               end if
                            end if
                         end if
                         
                      end do
                   end do
                   
                end do
             end do

          end if
       end do
    end if
    
  end subroutine MoveCrossExc
  
  ! *****************************************************************
  ! 2OPTSTAR NEIGHBORHOOD
  ! *****************************************************************

  subroutine Move2OptStar(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
       i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2)

    ! This routine considers neighbors constructed by interchanging
    ! segments of varying size in routes i and i2.

    ! First case: Exchanges the end of both routes, starting at client
    ! j+1 in route i and client k+1 in route i2.

    ! Second case: One route consists in the begining of route i, up
    ! to client j, plus the flipped begining of route i2, up to client
    ! k. The other route consists in the flipped end of route i,
    ! starting at client j+1, plus the end of route i2, starting at
    ! client k+1.
    
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: recomp
    integer, intent(in) :: i,i2
    integer, intent(inout) :: rnc,rnc2
    character, intent(in) :: otype,ptype
    real(kind=8), intent(inout) :: seed
    real(kind=8), intent(out) :: aggdem,aggdem2,length,length2
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat
    
    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients),pathind2(inst%nclients), &
         route2(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients),pathcoo2(2,inst%nclients)
    
    ! LOCAL SCALARS
    integer :: case,j,k,rnctrial,rnctrial2,rncorig,rncorig2
    real(kind=8) :: adtrial,adtrial2,ltrial,ltrial2
    
    ! LOCAL ARRAYS
    integer :: pindtrial(inst%nclients),rtrial(inst%nclients),pindtrial2(inst%nclients), &
         rtrial2(inst%nclients)
    real(kind=8) :: pcootrial(2,inst%nclients),pcootrial2(2,inst%nclients)
    
    rncorig = sol%nclients(i)
    rncorig2 = sol%nclients(i2)
    rnc = rncorig
    rnc2 = rncorig2
      
    if ( otype .eq. 'S' ) then
       case = min( 1 + floor( drand( seed ) * 2 ), 2 )
       
       j = min( floor( drand( seed ) * ( rnc  + 1 ) ), rnc  )
       k = min( floor( drand( seed ) * ( rnc2 + 1 ) ), rnc2 )
       if ( debug ) write(*,*) 'case = ',case,' i, j, i2, k = ',i,j,i2,k
       
       if ( case .eq. 1 ) then
          aggdem  = sum(inst%dem(sol%route(1:j,i))) + sum(inst%dem(sol%route(k+1:rnc2,i2)))
          aggdem2 = sum(inst%dem(sol%route(1:k,i2))) + sum(inst%dem(sol%route(j+1:rnc,i)))
          
          if ( aggdem .gt. inst%maxcap .or. aggdem2 .gt. inst%maxcap ) then
             aggdem  = - 1.0d0
             aggdem2 = - 1.0d0
             length  = - 1.0d0
             length2 = - 1.0d0
             
          else
             call SubstituteClients(stat,inst,sol,ptype,recomp,i2,k+1,rncorig2,.false.,i,j+1,rncorig,rnc,length,route,pathind, &
                  pathcoo)
             call SubstituteClients(stat,inst,sol,ptype,recomp,i,j+1,rncorig,.false.,i2,k+1,rncorig2,rnc2,length2,route2,pathind2, &
                  pathcoo2)
             
             if ( length .gt. inst%maxlen .or. length2 .gt. inst%maxlen ) then
                length  = - 1.0d0
                length2 = - 1.0d0
             end if
          end if
          
       else ! if (case .eq. 2 ) then
          aggdem  = sum(inst%dem(sol%route(1:j,i)))     + sum(inst%dem(sol%route(1:k,i2)))
          aggdem2 = sum(inst%dem(sol%route(j+1:rnc,i))) + sum(inst%dem(sol%route(k+1:rnc2,i2)))
          
          if ( aggdem .gt. inst%maxcap .or. aggdem2 .gt. inst%maxcap ) then
             aggdem  = - 1.0d0
             aggdem2 = - 1.0d0
             length  = - 1.0d0
             length2 = - 1.0d0
             
          else
             call SubstituteClients(stat,inst,sol,ptype,recomp,i2,1,k,.true.,i,j+1,rncorig,rnc,length,route,pathind,pathcoo)
             call SubstituteClients(stat,inst,sol,ptype,recomp,i,j+1,rncorig,.true.,i2,1,k,rnc2,length2,route2,pathind2,pathcoo2)

             if ( length .gt. inst%maxlen .or. length2 .gt. inst%maxlen ) then
                length  = - 1.0d0
                length2 = - 1.0d0
             end if
          end if
       end if
       
       if ( debug ) write(*,*) 'route = ',route(1:rnc),' length = ',length,' aggdem = ',aggdem
       if ( debug ) write(*,*) 'route2 = ',route2(1:rnc2),' length2 = ',length2,' aggdem2 = ',aggdem2
       
    else ! if ( otype .eq. 'B' .or. otype .eq. 'F' ) then
       length = sol%length(i)
       length2 = sol%length(i2)
       
       do j = 0,sol%nclients(i)
          do k = 0,sol%nclients(i2)
             do case = 1,2
                if ( case .eq. 1 ) then                   
                   ! Exchanges the end of both routes, starting at client j+1 in route i and client k+1 in route i2.
                   ! New route i : route(1:j,i)  + route(k+1:rncorig2,i2)
                   ! New route i2: route(1:k,i2) + route(j+1:rncorig,i)

                   adtrial  = sum(inst%dem(sol%route(1:j,i)))  + sum(inst%dem(sol%route(k+1:rncorig2,i2)))
                   adtrial2 = sum(inst%dem(sol%route(1:k,i2))) + sum(inst%dem(sol%route(j+1:rncorig,i)))
                   
                   if ( debug ) write(*,*) 'j = ', j,' k = ',k,' case = ',case,' adtrial = ',adtrial,' adtrial2 = ',adtrial2
                   
                   if ( adtrial .le. inst%maxcap .and. adtrial2 .le. inst%maxcap ) then
                      rnctrial = sol%nclients(i)
                      call SubstituteClients(stat,inst,sol,ptype,recomp,i2,k+1,rncorig2,.false.,i,j+1,rncorig,rnctrial,ltrial, &
                           rtrial,pindtrial,pcootrial)

                      if ( debug ) write(*,*) 'rtrial = ',rtrial(1:rnctrial),' ltrial = ',ltrial
                      
                      rnctrial2 = sol%nclients(i2)
                      call SubstituteClients(stat,inst,sol,ptype,recomp,i,j+1,rncorig,.false.,i2,k+1,rncorig2,rnctrial2,ltrial2, &
                           rtrial2, pindtrial2,pcootrial2)

                      if ( debug ) write(*,*) 'rtrial2 = ',rtrial2(1:rnctrial2),' ltrial2 = ',ltrial2
                      
                      if ( ltrial .le. inst%maxlen .and. ltrial2 .le. inst%maxlen ) then
                         if ( ltrial + ltrial2 .lt. length + length2 ) then
                            length = ltrial
                            aggdem = adtrial
                            rnc  = rnctrial
                            route(1:rnctrial) = rtrial(1:rnctrial)
                            pathind(1:rnctrial) = pindtrial(1:rnctrial)
                            pathcoo(1:2,1:rnctrial) = pcootrial(1:2,1:rnctrial)
                            
                            length2 = ltrial2
                            aggdem2 = adtrial2
                            rnc2 = rnctrial2
                            route2(1:rnctrial2) = rtrial2(1:rnctrial2)
                            pathind2(1:rnctrial2) = pindtrial2(1:rnctrial2)
                            pathcoo2(1:2,1:rnctrial2) = pcootrial2(1:2,1:rnctrial2)
                            if ( otype .eq. 'F' ) return
                         end if
                      end if
                   end if
                   
                else ! if ( case .eq. 2 ) then

                   ! One route consists in the begining of route i, up to client j, plus the flipped begining of route i2, up to client k.
                   ! The other route consists in the flipped end of route i, starting at client j+1, plus the end of route i2, starting at client k+1.
                   ! New route i: route(1:j,i) + flipped(route(1:k,i2)
                   ! New route i2: flipped(route(j+1:rncorig,i)) + route(k+1:rncorig2,i2)

                   adtrial  = sum(inst%dem(sol%route(1:j,i)))         + sum(inst%dem(sol%route(1:k,i2)))
                   adtrial2 = sum(inst%dem(sol%route(j+1:rncorig,i))) + sum(inst%dem(sol%route(k+1:rncorig2,i2)))

                   if ( debug ) write(*,*) 'j = ', j,' k = ',k,' case = ',case,' adtrial = ',adtrial,' adtrial2 = ',adtrial2
                   
                   if ( adtrial .le. inst%maxcap .and. adtrial2 .le. inst%maxcap ) then
                      rnctrial = sol%nclients(i)
                      call SubstituteClients(stat,inst,sol,ptype,recomp,i2,1,k,.true.,i,j+1,rncorig,rnctrial,ltrial,rtrial, &
                           pindtrial,pcootrial)
                      
                      if ( debug ) write(*,*) 'rtrial = ',rtrial(1:rnctrial),' ltrial = ',ltrial

                      rnctrial2 = sol%nclients(i2)
                      call SubstituteClients(stat,inst,sol,ptype,recomp,i,j+1,rncorig,.true.,i2,1,k,rnctrial2,ltrial2,rtrial2, &
                           pindtrial2, pcootrial2)
                   
                      if ( debug ) write(*,*) 'rtrial2 = ',rtrial2(1:rnctrial2),' ltrial2 = ',ltrial2                

                      if ( ltrial .le. inst%maxlen .and. ltrial2 .le. inst%maxlen ) then
                         if ( ltrial + ltrial2 .lt. length + length2 ) then
                            length = ltrial
                            aggdem = adtrial
                            rnc  = rnctrial
                            route(1:rnctrial) = rtrial(1:rnctrial)
                            pathind(1:rnctrial) = pindtrial(1:rnctrial)
                            pathcoo(1:2,1:rnctrial) = pcootrial(1:2,1:rnctrial)
                            
                            length2 = ltrial2
                            aggdem2 = adtrial2
                            rnc2 = rnctrial2
                            route2(1:rnctrial2) = rtrial2(1:rnctrial2)
                            pathind2(1:rnctrial2) = pindtrial2(1:rnctrial2)
                            pathcoo2(1:2,1:rnctrial2) = pcootrial2(1:2,1:rnctrial2)
                            if ( otype .eq. 'F' ) return
                         end if
                      end if
                   end if
                end if
             end do
          end do
       end do
       
    end if
    
  end subroutine Move2OptStar

  ! *****************************************************************
  ! RELOCATESTAR NEIGHBORHOOD
  ! *****************************************************************

  subroutine MoveReExcPlus(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
       i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2,nrini,nrend,flipallowed)
    
    ! This routine considers neighbors constructed by removing every
    ! segment of size between nrini and nrend in routes i and
    ! inserting it in every possible position in route i2,
        
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: flipallowed,recomp
    integer, intent(in) :: i,i2,nrini,nrend
    integer, intent(inout) :: rnc,rnc2
    character, intent(in) :: otype,ptype
    real(kind=8), intent(inout) :: seed
    real(kind=8), intent(out) :: aggdem,aggdem2,length,length2
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(statistics_type), intent(inout) :: stat
    
    ! ARRAY ARGUMENTS
    integer, intent(out) :: pathind(inst%nclients),route(inst%nclients),pathind2(inst%nclients), &
         route2(inst%nclients)
    real(kind=8), intent(out) :: pathcoo(2,inst%nclients),pathcoo2(2,inst%nclients)

    ! LOCAL SCALARS
    integer :: j,k,s,op1,op2,op1max,op2max,nr,nrendt,rnctrial,rnctrial2
    real(kind=8) :: adtrial,adtrial2,ltrial,ltrial2,tmp
    
    ! LOCAL ARRAYS
    logical :: flip(2),plus(2)
    integer :: pindtrial(inst%nclients),rtrial(inst%nclients),pindtrial2(inst%nclients), &
         rtrial2(inst%nclients)
    real(kind=8) :: pcootrial(2,inst%nclients),pcootrial2(2,inst%nclients)

    flip(1:2) = (/ .false., .true. /)
    plus(1:2) = (/ .false., .true. /)
    
    if ( otype .eq. 'S' ) then
       nrendt = min( nrend, rnc )
       nr  = min( nrini + floor( drand( seed ) * ( nrendt - nrini + 1 ) ), nrendt )       
       j = min( 1 + floor( drand( seed ) * ( rnc + 1 - nr ) ), rnc + 1 - nr )
       
       if ( rnc2 .gt. 0 ) then
          k = min( floor( drand( seed ) * rnc2 ), rnc2 - 1 )
          s = min( k + 1 + floor( drand( seed ) * ( rnc2 - ( k + 1 ) + 1  ) ), rnc2 )
          if ( .not.flipallowed ) then
             op1 = 1
          else
             op1 = min( 1 + floor( drand( seed ) * min( nr, 2 ) ), min( nr, 2 ) )
          end if
          op2 = min( 1 + floor( drand( seed ) * 2 ), 2 )
          if ( debug ) write(*,*) 'nr = ',nr,' j = ',j,' k = ',k,' s = ',s,' flip = ',flip(op1),' plus = ',plus(op2)
       else
          if ( debug ) write(*,*) 'nr = ',nr,' j = ',j,' k = 1',' flip = ',flip(1),' plus = ',plus(1)
       end if
          
       tmp = sum(inst%dem(sol%route(j:j+nr-1,i)))
       aggdem = sol%aggdem(i) - tmp
       aggdem2 = sol%aggdem(i2) + tmp
       
       if ( aggdem .gt. inst%maxcap .or. aggdem2 .gt. inst%maxcap ) then
          aggdem  = - 1.0d0
          aggdem2 = - 1.0d0
          length  = - 1.0d0
          length2 = - 1.0d0
          
       else
          call RemoveClients(stat,inst,sol,ptype,recomp,i,j,j+nr-1,rnc,length,route,pathind,pathcoo)
          if ( debug ) write(*,*) 'route = ',route(1:rnc),' length = ',length,' aggdem = ',aggdem

          if ( rnc2 .gt. 0 ) then
             call InsertClientsPlus(stat,inst,sol,ptype,recomp,i,j,nr,flip(op1),plus(op2),i2,k,s,rnc2,length2,route2,pathind2, &
                  pathcoo2)
             if ( debug ) write(*,*) 'route2 = ',route2(1:rnc2),' length2 = ',length2,' aggdem2 = ',aggdem2
          else
             call InsertClients(stat,inst,sol,ptype,recomp,i,j,j+nr-1,.false.,i2,1,rnc2,length2,route2,pathind2,pathcoo2)
             if ( debug ) write(*,*) 'route2 = ',route2(1:rnc2),' length2 = ',length2,' aggdem2 = ',aggdem2
          end if
                       
          if ( length .gt. inst%maxlen .or. length2 .gt. inst%maxlen ) then
             length  = - 1.0d0
             length2 = - 1.0d0
          end if
       end if
                
    else ! if ( otype .eq. 'B' .or. otype .eq. 'F' ) then
       do nr = nrini,nrend
          length = sol%length(i)
          length2 = sol%length(i2)
          
          do j = 1,sol%nclients(i) + 1 - nr
             rnctrial = sol%nclients(i)
             
             tmp = sum(inst%dem(sol%route(j:j+nr-1,i)))
             adtrial = sol%aggdem(i) - tmp
             adtrial2 = sol%aggdem(i2) + tmp
             
             if ( debug ) write(*,*) 'nr = ',nr,' j = ',j,' adtrial = ',adtrial,' adtrial2 = ',adtrial2
             
             if ( adtrial .le. inst%maxcap .and. adtrial2 .le. inst%maxcap ) then

                call RemoveClients(stat,inst,sol,ptype,recomp,i,j,j+nr-1,rnctrial,ltrial,rtrial,pindtrial,pcootrial)
                
                if ( debug ) write(*,*) 'rtrial = ',rtrial(1:rnctrial),' ltrial = ',ltrial    

                if ( sol%nclients(i2) .gt. 0 ) then
                   do k = 0,sol%nclients(i2) - 1                
                      do s = k + 1,sol%nclients(i2)
                         if ( .not. flipallowed ) then
                            op1max = 1
                         else
                            op1max = 2
                         end if
                         
                         do op1 = 1,min( nr, op1max )
                            if ( nr .eq. 1 .and. rnc2 .eq. 1 ) then
                               op2max = 1
                            else
                               op2max = 2
                            end if
                            
                            do op2 = 1,op2max
                               rnctrial2 = sol%nclients(i2)                    
                               
                               call InsertClientsPlus(stat,inst,sol,ptype,recomp,i,j,nr,flip(op1),plus(op2),i2,k,s,rnctrial2, &
                                    ltrial2,rtrial2,pindtrial2,pcootrial2)
                               
                               if ( debug ) write(*,*) 'k = ',k,' s = ',s,' rtrial2 = ',rtrial2(1:rnctrial2),' ltrial2 = ',ltrial2     
                               
                               if ( ltrial .le. inst%maxlen .and. ltrial2 .le. inst%maxlen ) then
                                  if ( ltrial + ltrial2 .lt. length + length2 ) then
                                     length = ltrial
                                     aggdem = adtrial
                                     rnc  = rnctrial
                                     route(1:rnctrial) = rtrial(1:rnctrial)
                                     pathind(1:rnctrial) = pindtrial(1:rnctrial)
                                     pathcoo(1:2,1:rnctrial) = pcootrial(1:2,1:rnctrial)
                                     
                                     length2 = ltrial2
                                     aggdem2 = adtrial2
                                     rnc2 = rnctrial2
                                     route2(1:rnctrial2) = rtrial2(1:rnctrial2)
                                     pathind2(1:rnctrial2) = pindtrial2(1:rnctrial2)
                                     pathcoo2(1:2,1:rnctrial2) = pcootrial2(1:2,1:rnctrial2)
                                     if ( otype .eq. 'F' ) return
                                  end if
                               end if
                            end do
                         end do
                      end do
                   end do

                else
                   rnctrial2 = 0
                   !length2 = 0.0d0
                   
                   call InsertClients(stat,inst,sol,ptype,recomp,i,j,j+nr-1,.false.,i2,1,rnctrial2,ltrial2,rtrial2, &
                        pindtrial2,pcootrial2)
                   
                   if ( debug ) write(*,*) 'rtrial2 = ',rtrial2(1:rnctrial2),' ltrial2 = ',ltrial2
                   
                   if ( ltrial .le. inst%maxlen .and. ltrial2 .le. inst%maxlen ) then
                      if ( ltrial + ltrial2 .lt. length + length2 ) then
                         length = ltrial
                         aggdem = adtrial
                         rnc  = rnctrial
                         route(1:rnctrial) = rtrial(1:rnctrial)
                         pathind(1:rnctrial) = pindtrial(1:rnctrial)
                         pathcoo(1:2,1:rnctrial) = pcootrial(1:2,1:rnctrial)
                         
                         length2 = ltrial2
                         aggdem2 = adtrial2
                         rnc2 = rnctrial2
                         route2(1:rnctrial2) = rtrial2(1:rnctrial2)
                         pathind2(1:rnctrial2) = pindtrial2(1:rnctrial2)
                         pathcoo2(1:2,1:rnctrial2) = pcootrial2(1:2,1:rnctrial2)
                         if ( otype .eq. 'F' ) return
                      end if
                   end if
                end if
                
             end if
          end do
       end do
    end if
    
  contains
    
    subroutine InsertClientsPlus(stat,inst,sol,ptype,recomp,isource,jini,nr,flip,plus,idest,kini,kend,rnc,length,route, &
         pathind,pathcoo)
      
      ! Insert after position kini of route idest clients from jini
      ! to jini+nr-1 of route isource.
      
      ! If flip = true, then the segment to be inserted in reversed.
      
      ! If plus = false, client kend of route idest remains where it
      ! is. Otherwise, client kend of route idest is relocated (in
      ! its own route) inbetween kini and the new segment that came
      ! from route isource.
      
      implicit none
      
      ! SCALAR ARGUMENTS
      logical, intent(in) :: flip,plus,recomp
      integer, intent(in) :: isource,idest,jini,kini,kend,nr
      integer, intent(inout) :: rnc
      character, intent(in) :: ptype
      real(kind=8), intent(out) :: length
      type(instance_type), intent(in) :: inst
      type(solution_type), intent(in) :: sol
      type(statistics_type), intent(inout) :: stat
      
      ! ARRAY ARGUMENTS
      integer, intent(out) :: pathind(inst%nclients),route(inst%nclients)
      real(kind=8), intent(out) :: pathcoo(2,inst%nclients)
      
      ! LOCAL SCALARS
      integer :: ell,rncorig
      
      rncorig = rnc
      rnc = rnc + nr
      
      if ( .not. flip ) then
         if ( .not. plus ) then
            route(1:rnc) = (/ sol%route(1:kini,idest), sol%route(jini:jini+nr-1,isource), sol%route(kend,idest), &
                 sol%route(kini+1:kend-1,idest), sol%route(kend+1:rncorig,idest) /)
         else
            route(1:rnc) = (/ sol%route(1:kini,idest), sol%route(kend,idest), sol%route(jini:jini+nr-1,isource), &
                 sol%route(kini+1:kend-1,idest), sol%route(kend+1:rncorig,idest) /)
         end if
      else
         if ( .not. plus ) then
            route(1:rnc) = (/ sol%route(1:kini,idest), sol%route(jini+nr-1:jini:-1,isource), sol%route(kend,idest), &
                 sol%route(kini+1:kend-1,idest), sol%route(kend+1:rncorig,idest) /)
         else
            route(1:rnc) = (/ sol%route(1:kini,idest), sol%route(kend,idest), sol%route(jini+nr-1:jini:-1,isource), &
                 sol%route(kini+1:kend-1,idest), sol%route(kend+1:rncorig,idest) /)
         end if
      end if
      
      if ( recomp ) then 
         if ( ptype .eq. 'G' ) then
            call DijkstraRoute(inst,rnc,route,pathind,length)
            do ell = 1,rnc
               pathcoo(1:2,ell) = inst%p(1:2,pathind(ell),route(ell))
            end do
         else ! if ( type .eq. 'N' ) then
            if ( .not. flip ) then
               if ( .not. plus ) then
                  pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini,idest), sol%pathcoo(1:2,jini:jini+nr-1,isource), &
                       sol%pathcoo(1:2,kend,idest), sol%pathcoo(1:2,kini+1:kend-1,idest), &
                       sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )
               else
                  pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini,idest), sol%pathcoo(1:2,kend,idest), &
                       sol%pathcoo(1:2,jini:jini+nr-1,isource), sol%pathcoo(1:2,kini+1:kend-1,idest), &
                       sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )
               end if
            else
               if ( .not. plus ) then
                  pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini,idest), sol%pathcoo(1:2,jini+nr-1:jini:-1,isource), &
                       sol%pathcoo(1:2,kend,idest), sol%pathcoo(1:2,kini+1:kend-1,idest), &
                       sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )
               else
                  pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini,idest), sol%pathcoo(1:2,kend,idest), &
                       sol%pathcoo(1:2,jini+nr-1:jini:-1,isource), sol%pathcoo(1:2,kini+1:kend-1,idest), &
                       sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )                   
               end if
            end if
            call BCD(stat,inst,rnc,route,pathcoo,length)
            pathind(1:rnc) = - 1
         end if
         
      else ! if ( .not. recomp ) then
         if ( .not. flip ) then
            if( .not. plus ) then
               pathind(1:rnc) = (/ sol%pathind(1:kini,idest), sol%pathind(jini:jini+nr-1,isource),sol%pathind(kend,idest), &
                    sol%pathind(kini+1:kend-1,idest), sol%pathind(kend+1:rncorig,idest) /)
               pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini,idest), sol%pathcoo(1:2,jini:jini+nr-1,isource), &
                    sol%pathcoo(1:2,kend,idest), sol%pathcoo(1:2,kini+1:kend-1,idest), &
                    sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )
            else
               pathind(1:rnc) = (/ sol%pathind(1:kini,idest), sol%pathind(kend,idest), sol%pathind(jini:jini+nr-1,isource), &
                    sol%pathind(kini+1:kend-1,idest), sol%pathind(kend+1:rncorig,idest) /)
               pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini,idest), sol%pathcoo(1:2,kend,idest), &
                    sol%pathcoo(1:2,jini:jini+nr-1,isource), sol%pathcoo(1:2,kini+1:kend-1,idest), &
                    sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )                
            end if
         else
            if( .not. plus ) then
               pathind(1:rnc) = (/ sol%pathind(1:kini,idest), sol%pathind(jini+nr-1:jini:-1,isource),sol%pathind(kend,idest), &
                    sol%pathind(kini+1:kend-1,idest), sol%pathind(kend+1:rncorig,idest) /)
               pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini,idest), sol%pathcoo(1:2,jini+nr-1:jini:-1,isource), &
                    sol%pathcoo(1:2,kend,idest), sol%pathcoo(1:2,kini+1:kend-1,idest), &
                    sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )
            else
               pathind(1:rnc) = (/ sol%pathind(1:kini,idest), sol%pathind(kend,idest), sol%pathind(jini+nr-1:jini:-1,isource), &
                    sol%pathind(kini+1:kend-1,idest), sol%pathind(kend+1:rncorig,idest) /)
               pathcoo(1:2,1:rnc) = reshape( (/ sol%pathcoo(1:2,1:kini,idest), sol%pathcoo(1:2,kend,idest), &
                    sol%pathcoo(1:2,jini+nr-1:jini:-1,isource), sol%pathcoo(1:2,kini+1:kend-1,idest), &
                    sol%pathcoo(1:2,kend+1:rncorig,idest) /), (/ 2, rnc /) )
            end if
         end if
         
         if ( ptype .eq. 'G' ) then
            length = inst%dist(pathind(1),route(1),1,0)
            do ell = 2,rnc
               length = length + inst%dist(pathind(ell),route(ell),pathind(ell-1),route(ell-1))
            end do
            length = length + inst%dist(1,0,pathind(rnc),route(rnc))
            
         else ! if ( ptype .eq. 'N' ) then
            length = norm2( pathcoo(1:2,1) - inst%p(1:2,1,0) )
            do ell = 2,rnc
               length = length + norm2( pathcoo(1:2,ell) - pathcoo(1:2,ell-1) )
            end do
            length = length + norm2( inst%p(1:2,1,0) - pathcoo(1:2,rnc) )
         end if
      end if
      
    end subroutine InsertClientsPlus
    
  end subroutine MoveReExcPlus
  
  ! *****************************************************************
  ! MOVEMENT
  ! *****************************************************************
  
  subroutine Shake(stat,inst,sol,ptype,recomp,nrend,flipallowed,seed,singleroute,feasible,moveind)
    
    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: flipallowed,recomp,singleroute
    logical, intent(out) :: feasible
    character, intent(in) :: ptype
    integer, intent(in) :: moveind,nrend
    real(kind=8), intent(inout) :: seed
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(inout) :: sol
    type(statistics_type), intent(inout) :: stat
    
    ! PARAMETER
    integer, parameter :: maxtrials = 100, ncombsizemax = 100    
    
    ! LOCAL SCALARS
    integer :: fmax,fmax2,i,itemp,itrial,i2,i2temp,j,ncombsize,ner,nr,nr2,rnc,rnc2
    real(kind=8) :: aggdem,aggdem2,length,length2
    
    ! LOCAL ARRAYS
    integer :: pathind(inst%nclients),route(inst%nclients),pathind2(inst%nclients),route2(inst%nclients), &
         routesID(inst%nclients),combsize(2,ncombsizemax)
    real(kind=8) :: pathcoo(2,inst%nclients),pathcoo2(2,inst%nclients)
    
    ner = count( sol%nclients(1:sol%nroutes) .gt. 0 )
    routesID(1:ner) = pack( (/ (j, j=1,sol%nroutes) /), sol%nclients(1:sol%nroutes) .gt. 0 )
    if ( debug ) write(*,*) 'ner = ',ner,' routesID = ',routesID(1:ner)
    
    itrial = 0
    
100 continue
    
    itrial = itrial + 1
    
    itemp = min( 1 + floor( drand( seed ) * ner ), ner )
    i = routesID(itemp)           
    rnc = sol%nclients(i)
    
    if ( singleroute ) then
       if ( .not. ( 1 .le. moveind .and. moveind .le. 3 ) ) then
          write(*,*) 'In Shake, moveind must be an integer between 1 and 3 when singleroute = TRUE!'
          stop
       end if
       
       if ( moveind .eq. 1 ) then
          if ( rnc .ge. 2 ) then
             call Move2Opt(stat,inst,sol,'S',ptype,recomp,seed,i,rnc,length,route,pathind,pathcoo)
          else
             length = - 1.0d0
          end if
          
       else if ( moveind .eq. 2 ) then
          if ( rnc .ge. 3 ) then
             call Move3Opt(stat,inst,sol,'S',ptype,recomp,seed,i,rnc,length,route,pathind,pathcoo)
          else
             length = - 1.0d0
          end if
          
       else if ( moveind .eq. 3 ) then
          call MoveOrOpt(stat,inst,sol,'S',ptype,recomp,seed,i,rnc,length,route,pathind,pathcoo,1,nrend,flipallowed) ! nrend=4 (with nrini=1, nrend>=nrini makes sense)
       end if
       
       feasible = length .ge. 0.0d0
       if ( debug ) write(*,*) 'feasible = ',feasible
       
    else
       if( ner .lt. 2 .and. ( moveind .eq. 5 .or. moveind .eq. 6 .or. moveind .eq. 9 ) ) then
          feasible = .false.
          if ( debug ) write(*,*) 'There is only one route in the current solution. Shakes of type 5, 6 and 9 are not possible.'
          return
       end if
       
       if ( .not. ( 4 .le. moveind .and. moveind .le. 9 ) ) then
          write(*,*) 'In Shake, moveind must be an integer between 4 and 9 when singleroute = FALSE!'
          stop         
       end if

       if ( moveind .eq. 5 .or. moveind .eq. 6 .or. moveind .eq. 9 .or. ner .eq. sol%nroutes ) then
          i2temp = min( 1 + floor( drand( seed ) * ( ner - 1 ) ), ner - 1 )
          if ( i2temp .ge. itemp ) i2temp = i2temp + 1
          i2 = routesID(i2temp) 
       else
          i2temp = min( 1 + floor( drand( seed ) * ner ), ner )
          if ( i2temp .eq. itemp ) then
             i2 = minloc( sol%nclients(1:sol%nroutes), dim=1 )
          else
             i2 = routesID(i2temp) 
          end if
       end if
       
       rnc2 = sol%nclients(i2)
          
       if ( moveind .eq. 4 ) then
          call Move10Exc(stat,inst,sol,'S',ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
               i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2)
          
       else if ( moveind .eq. 5 ) then
          call Move11Exc(stat,inst,sol,'S',ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
               i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2)
          
       else if ( moveind .eq. 6 ) then
          call Move12Exc(stat,inst,sol,'S',ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
               i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2,flipallowed)
          
       else if ( moveind .eq. 7 ) then
          call MoveReExcPlus(stat,inst,sol,'S',ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
               i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2,2,nrend,flipallowed) ! nrend=4 (with nrini=2, nrend>=nrini makes sense)
          
       else if ( moveind .eq. 8 ) then
          call Move2OptStar(stat,inst,sol,'S',ptype,recomp,seed,i,rnc,length,aggdem,route,pathind,pathcoo, &
               i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2)
          
       else if ( moveind .eq. 9 ) then
          ncombsize = 0
          do nr = 1,nrend
             do nr2 = 1,nrend          
                if ( .not. ( ( nr .eq. 1 .and. nr2 .eq. 1 ) .or. ( nr .eq. 1 .and. nr2 .eq. 2 ) .or. &
                     ( nr .eq. 2 .and. nr2 .eq. 1) ) ) then                   
                   ncombsize = ncombsize + 1
                   if ( ncombsize .gt. ncombsizemax ) then
                      write(*,*) 'In Neighbor: increase ncombsizemax and re-run.'
                      stop
                   end if
                   
                   combsize(1,ncombsize) = nr
                   combsize(2,ncombsize) = nr2
                end if
             end do
          end do

          if ( debug ) write(*,*) 'ncombsize = ',ncombsize,' combsize = ',combsize(1:2,1:ncombsize)
          
          if ( .not. flipallowed ) then
             fmax  = 0
             fmax2 = 0
          else
             fmax  = 1
             fmax2 = 1
          end if
          
          call MoveCrossExc(stat,inst,sol,'S',ptype,recomp,seed,ncombsize,combsize,0,fmax,i,rnc,length,aggdem,route, &
               pathind,pathcoo,0,fmax2,i2,rnc2,length2,aggdem2,route2,pathind2,pathcoo2) ! nrend=3 (with nrini=2, nrend>=nrini makes sense) 
       end if
       
       feasible = aggdem .ge. 0.0d0 .and. aggdem2 .ge. 0.0d0 .and. length .ge. 0.0d0 .and. length2 .ge. 0.0d0
       if ( debug ) write(*,*) 'feasible = ',feasible
    end if
    
    if ( .not. feasible ) then
       if ( itrial .lt. maxtrials ) then
          go to 100
       else
          if ( debug ) write(*,*) 'Shake did not find a feasible solution with moveind = ',moveind
       end if
       
    else       
       ! Apply the movement
       if ( singleroute ) then
          sol%length(i) = length
          sol%route(1:rnc,i) = route(1:rnc)
          sol%pathind(1:rnc,i) = pathind(1:rnc)
          sol%pathcoo(1:2,1:rnc,i) = pathcoo(1:2,1:rnc)
       else
          sol%nclients(i) = rnc
          sol%length(i) = length
          sol%aggdem(i) = aggdem
          sol%route(1:rnc,i) = route(1:rnc)
          sol%pathind(1:rnc,i) = pathind(1:rnc)
          sol%pathcoo(1:2,1:rnc,i) = pathcoo(1:2,1:rnc)
          
          sol%nclients(i2) = rnc2
          sol%length(i2) = length2
          sol%aggdem(i2) = aggdem2
          sol%route(1:rnc2,i2) = route2(1:rnc2)
          sol%pathind(1:rnc2,i2) = pathind2(1:rnc2)
          sol%pathcoo(1:2,1:rnc2,i2) = pathcoo2(1:2,1:rnc2)             
       end if
    end if
    
  end subroutine Shake
  
  !*************************************************************
  !*************************************************************

  subroutine Neighbor(stat,inst,sol,otype,ptype,recomp,nrend,flipallowed,seed,singleroute,improve,moveind)
    
    implicit none
       
    ! SCALAR ARGUMENTS
    logical, intent(in) :: flipallowed,recomp,singleroute
    logical, intent(out) :: improve
    character, intent(in) :: otype,ptype
    integer, intent(in) :: moveind,nrend
    real(kind=8), intent(inout) :: seed
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(inout) :: sol
    type(statistics_type), intent(inout) :: stat

    ! PARAMETERS
    integer, parameter :: ncombsizemax = 100
    
    ! LOCAL SCALARS
    logical :: emptyroute
    integer :: fmax,fmax2,i,i2,i2ini,ncombsize,nr,nr2,rnc,rnc2,ibest,ibest2,rncbest,rncbest2
    real(kind=8) :: aggdem,aggdem2,adtrial,adtrial2,delta,ltrial,ltrial2,length,length2
    
    ! LOCAL ARRAYS
    integer :: pathind(inst%nclients),pathind2(inst%nclients),route(inst%nclients),route2(inst%nclients), &
         pindtrial(inst%nclients),pindtrial2(inst%nclients),rtrial(inst%nclients),rtrial2(inst%nclients), &
         combsize(2,ncombsizemax)
    real(kind=8) :: pathcoo(2,inst%nclients),pathcoo2(2,inst%nclients), pcootrial(2,inst%nclients), &
         pcootrial2(2,inst%nclients)
   
    if ( singleroute ) then
       if ( .not. ( 1 .le. moveind .and. moveind .le. 3 ) ) then
          write(*,*) 'In Neighbor, moveind must be an integer between 1 and 3 when singleroute = TRUE!'
          stop
       end if
       
       i = 1
       improve = .false.
       do while ( i .le. sol%nroutes .and. .not. ( otype .eq. 'F' .and. improve ) )
          rnc = sol%nclients(i)
          if ( rnc .gt. 0 ) then             
             length = sol%length(i)
             
             if ( moveind .eq. 1 ) then
                call Move2Opt(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,ltrial,rtrial,pindtrial,pcootrial)
                
             else if ( moveind .eq. 2 ) then
                call Move3Opt(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,ltrial,rtrial,pindtrial,pcootrial)
                
             else if ( moveind .eq. 3 ) then
                call MoveOrOpt(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,ltrial,rtrial,pindtrial,pcootrial,1,nrend,flipallowed)          
             end if
             
             if ( ltrial .lt. length * ( 1.0d0 - sqrt( epsilon(1.0d0) ) ) ) then
                improve = .true.
                sol%length(i) = ltrial
                sol%route(1:rnc,i) = rtrial(1:rnc)
                sol%pathind(1:rnc,i) = pindtrial(1:rnc)
                sol%pathcoo(1:2,1:rnc,i) = pcootrial(1:2,1:rnc)
                if ( debug ) write(*,*) 'A better neighbor was found with moveind = ',moveind
             else
                if ( debug ) write(*,*) 'No better neighbor exists with moveind = ',moveind
             end if
          end if
          i = i + 1
       end do
       
    else     
       if ( .not. ( 4 .le. moveind .and. moveind .le. 9 ) ) then
          write(*,*) 'In Neighbor, moveind must be an integer between 4 and 9 when singleroute = FALSE!'
          stop
       end if
       
       improve = .false.
       delta = huge(1.0d0)
       outer_loop: do i = 1,sol%nroutes
          rnc = sol%nclients(i)
          if ( rnc .gt. 0 ) then
             emptyroute = .false.                
             if ( moveind .ne. 5 ) then
                i2ini = 1
             else
                i2ini = i + 1
             end if
             inner_loop: do i2 = i2ini,sol%nroutes
                if ( i2 .ne. i ) then
                   rnc2 = sol%nclients(i2)
                   if ( rnc2 .gt. 0 .or. ( ( moveind .eq. 4 .or. moveind .eq. 7 ) .and. .not. emptyroute ) ) then
                      if ( rnc2 .eq. 0 ) emptyroute = .true.       
                      
                      if ( moveind .eq. 4 ) then
                         call Move10Exc(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,ltrial,adtrial,rtrial,pindtrial,pcootrial, &
                              i2,rnc2,ltrial2,adtrial2,rtrial2,pindtrial2,pcootrial2)
                         
                      else if ( moveind .eq. 5 ) then
                         call Move11Exc(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,ltrial,adtrial,rtrial,pindtrial,pcootrial, &
                              i2,rnc2,ltrial2,adtrial2,rtrial2,pindtrial2,pcootrial2)
                         
                      else if ( moveind .eq. 6 ) then
                         call Move12Exc(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,ltrial,adtrial,rtrial,pindtrial,pcootrial, &
                              i2,rnc2,ltrial2,adtrial2,rtrial2,pindtrial2,pcootrial2,flipallowed)
                         
                      else if ( moveind .eq. 7 ) then
                         call MoveReExcPlus(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,ltrial,adtrial,rtrial,pindtrial,pcootrial, &
                              i2,rnc2,ltrial2,adtrial2,rtrial2,pindtrial2,pcootrial2,2,nrend,flipallowed)   
                         
                      else if ( moveind .eq. 8 ) then
                         call Move2OptStar(stat,inst,sol,otype,ptype,recomp,seed,i,rnc,ltrial,adtrial,rtrial,pindtrial,pcootrial, &
                              i2,rnc2,ltrial2,adtrial2,rtrial2,pindtrial2,pcootrial2)
                         
                      else if ( moveind .eq. 9 ) then
                         ncombsize = 0
                         do nr = 1,nrend
                            do nr2 = 1,nrend
                               if ( .not. ( ( nr .eq. 1 .and. nr2 .eq. 1 ) .or. ( nr .eq. 1 .and. nr2 .eq. 2 ) .or. &
                                    ( nr .eq. 2 .and. nr2 .eq. 1) ) ) then
                                  ncombsize = ncombsize + 1
                                  if ( ncombsize .gt. ncombsizemax ) then
                                     write(*,*) 'In Neighbor: increase ncombsizemax and re-run.'
                                     stop
                                  end if
                                  
                                  combsize(1,ncombsize) = nr
                                  combsize(2,ncombsize) = nr2
                               end if
                            end do
                         end do
                         
                         if ( debug ) write(*,*) 'ncombsize = ',ncombsize,' combsize = ',combsize(1:2,1:ncombsize)
                         
                         if ( .not. flipallowed ) then
                            fmax  = 0
                            fmax2 = 0
                         else
                            fmax  = 1
                            fmax2 = 1
                         end if
                         
                         call MoveCrossExc(stat,inst,sol,otype,ptype,recomp,seed,ncombsize,combsize,0,fmax,i,rnc,ltrial,adtrial, &
                              rtrial,pindtrial,pcootrial,0,fmax2,i2,rnc2,ltrial2,adtrial2,rtrial2,pindtrial2,pcootrial2)
                      end if
                      
                      if ( ltrial + ltrial2 .lt. ( sol%length(i) + sol%length(i2) ) * ( 1.0d0 - sqrt( epsilon(1.0d0) ) ) ) then
                         if ( ltrial + ltrial2 .lt. delta ) then
                            improve = .true.
                            delta = ltrial + ltrial2
                            
                            ibest = i
                            rncbest = rnc
                            length = ltrial
                            aggdem = adtrial
                            route(1:rncbest) = rtrial(1:rncbest)
                            pathind(1:rncbest) = pindtrial(1:rncbest)
                            pathcoo(1:2,1:rncbest) = pcootrial(1:2,1:rncbest)
                            
                            ibest2 = i2
                            rncbest2 = rnc2
                            length2 = ltrial2
                            aggdem2 = adtrial2
                            route2(1:rncbest2) = rtrial2(1:rncbest2)
                            pathind2(1:rncbest2) = pindtrial2(1:rncbest2)
                            pathcoo2(1:2,1:rncbest2) = pcootrial2(1:2,1:rncbest2)                         
                            if ( otype .eq. 'F' ) exit outer_loop
                         end if
                      end if
                      
                   end if
                end if
             end do inner_loop
          end if
       end do outer_loop
       
       if ( improve ) then
          sol%nclients(ibest) = rncbest
          sol%length(ibest) = length
          sol%aggdem(ibest) = aggdem
          sol%route(1:rncbest,ibest) = route(1:rncbest)
          sol%pathind(1:rncbest,ibest) = pathind(1:rncbest)
          sol%pathcoo(1:2,1:rncbest,ibest) = pathcoo(1:2,1:rncbest)
          
          sol%nclients(ibest2) = rncbest2
          sol%length(ibest2) = length2
          sol%aggdem(ibest2) = aggdem2
          sol%route(1:rncbest2,ibest2) = route2(1:rncbest2)
          sol%pathind(1:rncbest2,ibest2) = pathind2(1:rncbest2)
          sol%pathcoo(1:2,1:rncbest2,ibest2) = pathcoo2(1:2,1:rncbest2)
          
          if ( debug ) write(*,*) 'A better neighbor was found with moveind = ',moveind
       else
          if ( debug ) write(*,*) 'No better neighbor exists with moveind = ',moveind
       end if
    end if
    
  end subroutine Neighbor
  
end module neighborhoods
