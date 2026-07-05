module solution

  use instance, only: instance_type
  
  implicit none

  logical, parameter :: debug = .false.
  
  ! nroutes is the number of routes
  ! nclients(r) is the number of clients of route r
  ! route(1:nclients(r),r) is the list of clients of route r
  ! pathind(1:nclients(r),r) index of the client's point that is visited; -1 if the point isn't a client point
  ! pathcoo(1:2,1:nclients(r),r) Cartesian coordinate of the point associated with each client of the route
  ! clientroute(c) indicates to which route client c belongs
  ! length(r) saves the length of route r
  ! aggdem(r) saves the aggregated demand of route r
  
  type solution_type
     integer :: nroutes
     integer, allocatable :: nclients(:),route(:,:),pathind(:,:),clientroute(:)
     real(kind=8), allocatable :: length(:),aggdem(:),pathcoo(:,:,:)
  end type solution_type

  public :: solution_type,CreateSolution,DeleteSolution,CopySolution,drawsol
  
contains
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine CreateSolution(inst,sol)

    implicit none
    
    ! SCALAR ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(out) :: sol

    ! LOCAL SCALARS
    integer :: allocerr
    
    allocate(sol%nclients(inst%nclients), &
         sol%route(inst%nclients,inst%nclients), &
         sol%pathind(inst%nclients,inst%nclients), &
         sol%pathcoo(2,inst%nclients,inst%nclients), &
         sol%clientroute(inst%nclients), &
         sol%length(inst%nclients), &
         sol%aggdem(inst%nclients), &
         stat=allocerr)

    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error.'
       stop
    end if

  end subroutine CreateSolution
    
  ! *****************************************************************
  ! *****************************************************************
  
  subroutine DeleteSolution(sol)
    
    implicit none
    
    ! SCALAR ARGUMENTS
    type(solution_type), intent(inout) :: sol
    
    ! LOCAL SCALARS
    integer :: allocerr
    
    deallocate(sol%nclients,sol%route,sol%pathind,sol%pathcoo,sol%clientroute,sol%length,sol%aggdem,stat=allocerr)
    
    if ( allocerr .ne. 0 ) then
       write(*,*) 'Deallocation error.'
       stop
    end if
    
  end subroutine DeleteSolution
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine CopySolution(inst,sol,copy)

    implicit none
    
    ! SCALAR ARGUMENTS
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
    type(solution_type), intent(out) :: copy
    
    ! LOCAL SCALARS
    integer :: allocerr,i

    allocate(copy%nclients(inst%nclients), &
         copy%route(inst%nclients,inst%nclients), &
         copy%pathind(inst%nclients,inst%nclients), &
         copy%pathcoo(2,inst%nclients,inst%nclients), &
         copy%clientroute(inst%nclients), &
         copy%length(inst%nclients), &
         copy%aggdem(inst%nclients), &
         stat=allocerr)

    if ( allocerr .ne. 0 ) then
       write(*,*) 'Allocation error.'
       stop
    end if
    
    copy%nroutes = sol%nroutes
    
    do i = 1,copy%nroutes
       copy%nclients(i) = sol%nclients(i)
       copy%length(i) = sol%length(i)
       copy%aggdem(i) = sol%aggdem(i)
       copy%route(1:copy%nclients(i),i) = sol%route(1:sol%nclients(i),i)
       copy%pathind(1:copy%nclients(i),i) = sol%pathind(1:sol%nclients(i),i)
       copy%pathcoo(1:2,1:copy%nclients(i),i) = sol%pathcoo(1:2,1:sol%nclients(i),i)
    end do

    copy%clientroute(1:inst%nclients) = sol%clientroute(1:inst%nclients)
    
  end subroutine CopySolution
    
  ! **************************************************************
  ! **************************************************************
  
  subroutine drawsol(inst,sol,filename,soltype)
    
    ! SCALAR ARGUMENTS
    character, intent(in) :: soltype
    type(instance_type), intent(in) :: inst
    type(solution_type), intent(in) :: sol
        
    ! ARRAY ARGUMENTS
    character(len=80), intent(in) :: filename
    
    ! LOCAL SCALARS
    integer :: j,k,rnc
    real(kind=8) :: xymax

    xymax = 1.0d0
    do k = 1,inst%nclients
       do j = 1,inst%np(k)
          xymax = max( xymax, maxval( abs( inst%p(1:2,j,k) ) ) )
       end do
    end do

    open(unit=10,file=trim(filename))

    write(10,10)

    ! draw the "filled" polygons.
    if ( soltype .eq. 'N' ) then
       do k = 1,inst%nclients
          write(10,43)
          do j = 1,inst%np(k) - 1
             write(10,44) inst%p(1:2,j,k) / xymax
          end do
          write(10,45) inst%p(1:2,inst%np(k),k) / xymax
       end do
    end if
    
    ! drawing the routes
    if ( soltype .eq. 'G' ) then
       do k = 1,sol%nroutes
          rnc = sol%nclients(k)
          
          if ( rnc .gt. 0 ) then
             write(10,20) inst%p(1:2,1,0) / xymax,inst%p(1:2,sol%pathind(1,k),sol%route(1,k)) / xymax
             do j = 1,rnc - 1
                write(10,20) inst%p(1:2,sol%pathind(j,k),sol%route(j,k)) / xymax, &
                     inst%p(1:2,sol%pathind(j+1,k),sol%route(j+1,k)) / xymax
             end do
             write(10,20) inst%p(1:2,sol%pathind(rnc,k),sol%route(rnc,k)) / xymax,inst%p(1:2,1,0) / xymax
          end if
          
       end do
    end if

    if ( soltype .eq. 'N' ) then
       do k = 1,sol%nroutes
          rnc = sol%nclients(k)
          
          if ( rnc .gt. 0 ) then          
             write(10,20) inst%p(1:2,1,0) / xymax,sol%pathcoo(1:2,1,k) / xymax
             do j = 1,rnc - 1
                write(10,20) sol%pathcoo(1:2,j,k) / xymax,sol%pathcoo(1:2,j+1,k) / xymax
             end do
             write(10,20) sol%pathcoo(1:2,rnc,k) / xymax,inst%p(1:2,1,0) / xymax       
          end if
          
       end do
    end if

    write(10,30) inst%p(1:2,1,0) / xymax ! draw the depot
    
    ! draw the edges and vertices of the polygons.
    do k = 1,inst%nclients
       if ( soltype .eq. 'N' .or. soltype .eq. 'G' ) then
          do j = 1,inst%np(k) - 1
             write(10,40) inst%p(1:2,j,k) / xymax,inst%p(1:2,j+1,k) / xymax
          end do
          write(10,40) inst%p(1:2,inst%np(k),k) / xymax,inst%p(1:2,1,k) / xymax
       end if
       
       do j = 1,inst%np(k)
          write(10,32) inst%p(1:2,j,k) / xymax
       end do
    end do
    
    ! color the point to be visited in each polygon/client with a different color.
    do k = 1,sol%nroutes
       rnc = sol%nclients(k)
       if ( rnc .gt.0 ) then
          
          do j = 1,rnc
             if ( soltype .eq. 'G' ) then 
                write(10,33) inst%p(1:2,sol%pathind(j,k),sol%route(j,k)) / xymax
             end if
             if ( soltype .eq. 'N' ) then 
                write(10,33) sol%pathcoo(1:2,j,k) / xymax
             end if
          end do
          
       end if
    end do
    
    write(10,50)

10  format('prologues := 3;',/, &
           'outputtemplate := "%j.mps";',/, &
           'beginfig(1);',/, &
           'u = 7.0cm;',/, & 
           'defaultscale := 1.0;') 
20  format('draw( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 0.05pt withcolor blue;')
30  format('drawdot(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 0.5pt withcolor green;')
32  format('drawdot(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 0.25pt withcolor black;')
33  format('drawdot(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 0.25pt withcolor red;')
40  format('draw( ',f20.10,'u,',f20.10,'u)--(',f20.10,'u,',f20.10,'u) withpen pencircle scaled 0.1pt withcolor (0.75, 0.75, 0.75);')
43  format('fill')
44  format(' (', f20.10, 'u,', f20.10, 'u) -- ')
45  format(' (', f20.10, 'u,', f20.10, 'u) -- cycle withcolor (0.95, 0.95, 0.95);')

50  format('endfig;',/,'end;')

    close(10)
    
  end subroutine drawsol
  
end module solution
