module Dijkstra
  
  use instance, only: instance_type
  
  implicit none
  
  logical, parameter :: debug = .false.
  
  public :: DijkstraRoute
  
contains
  
  ! *****************************************************************
  ! *****************************************************************
  
  subroutine DijkstraRoute(inst,rnc,route,pathind,length)
    
    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: rnc
    real(kind=8), intent(out) :: length
    type(instance_type), intent(in) :: inst
    
    ! ARRAY ARGUMENTS
    integer, intent(in) :: route(rnc)
    integer, intent(out) :: pathind(rnc)
    
    ! LOCAL SCALARS
    integer :: i,j,k,last
    
    ! LOCAL ARRAYS
    integer :: prev(rnc+1,inst%npmax),routewd(0:rnc+1)
    real(kind=8) :: dist(0:rnc+1,inst%npmax)

    routewd(0:rnc+1) = (/ 0, route(1:rnc), 0 /)
   
    dist(1:rnc+1,1:inst%npmax) = huge(1.0d0)
    dist(0,1:inst%npmax) = 0.0d0
    
    do i = 1,rnc+1
       do k = 1,inst%np(routewd(i))
          ! k-th point of i-th client of the route
          do j = 1,inst%np(routewd(i-1))
             if ( dist(i-1,j) + inst%dist(j,routewd(i-1),k,routewd(i)) .lt. dist(i,k) ) then
                dist(i,k) = dist(i-1,j) + inst%dist(j,routewd(i-1),k,routewd(i))              
                prev(i,k) = j
             end if
          end do
       end do
    end do
    
    last = 1
    do i = rnc,1,-1
       last = prev(i+1,last)
       pathind(i) = last
    end do
    
    length = dist(rnc+1,1)
    
    if ( debug ) then
       write(*,*) 'Route: ',route(1:rnc)
       write(*,*) 'Path: ',pathind(1:rnc)
       write(*,*) 'Minimum distance: ',length
    end if
    
  end subroutine DijkstraRoute
  
end module Dijkstra
