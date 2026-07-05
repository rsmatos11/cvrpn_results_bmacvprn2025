module statistics

  implicit none

  integer, parameter :: hlenmax = 10000
  
  type statistics_type
     integer :: bcdncalls,bcdtotit,bcdmaxit,bcdstp(0:1),ncycles,ncsol,nneighs,nnsol,nshakes,nssol
     real(kind=8) :: tsol

     integer :: hlen
     real(kind=8) :: vtime(hlenmax),vsolval(hlenmax)
  end type statistics_type

  public :: statistics_type,CreateStatistics,CleanStatistics,hlenmax

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine CreateStatistics(stat)

    implicit none

    ! SCALAR ARGUMENTS
    type(statistics_type), intent(inout) :: stat

    stat%hlen = 0
    
    ! BCD statistics
    stat%bcdncalls = 0
    stat%bcdtotit = 0
    stat%bcdmaxit = 0
    stat%bcdstp(0:1) = 0

    ! GVNS statistics
    stat%ncycles = 0
    stat%nshakes = 0
    stat%nneighs = 0
    
    stat%tsol    = 0.0d0
    stat%ncsol   = 0
    stat%nssol   = 0
    stat%nnsol   = 0
    
  end subroutine CreateStatistics
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine CleanStatistics(stat)

    implicit none

    ! SCALAR ARGUMENTS
    type(statistics_type), intent(inout) :: stat

    stat%hlen = 0
    
    ! BCD statistics
    stat%bcdncalls = 0
    stat%bcdtotit = 0
    stat%bcdmaxit = 0
    stat%bcdstp(0:1) = 0

    ! GVNS statistics
    stat%ncycles = 0
    stat%nshakes = 0
    stat%nneighs = 0
    stat%tsol    = 0.0d0
    stat%ncsol   = 0
    stat%nssol   = 0
    stat%nnsol   = 0
    
  end subroutine CleanStatistics
  
  ! *****************************************************************
  ! *****************************************************************

  ! subroutine DeleteStatistics(stat)

  !   implicit none

  !   ! SCALAR ARGUMENTS
  !   type(statistics_type), intent(inout) :: stat

  ! end subroutine DeleteStatistics
  
end module statistics
