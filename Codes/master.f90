program master

  use instance, only: instance_type,CreateInstance,DeleteInstance
  use solution, only: solution_type,CreateSolution,DeleteSolution,CopySolution,drawsol
  use statistics, only: statistics_type,CreateStatistics,CleanStatistics
  use dijkstra, only: DijkstraRoute
  use constructiveheuristcs, only: CW,NN
  use blockcoordinatedescent, only: BCD
  use neighborhoods, only: Neighbor,Shake
  use metaheuristics, only: GVNS

  implicit none

  ! EXTERNAL FUNCTIONS
  real(kind=8), external :: drand 

  ! PARAMETERS
  logical, parameter :: debug = .false.
  real(kind=8), parameter :: seedbase = 123456.0d0,macheps12 = sqrt( epsilon( 1.0d0 ) )

  ! LOCAL SCALARS
  logical :: recomp,flip
  integer :: bestrun,i,k,ncycles,nneighs,norder,nrun,nshakes,qmax,rmax
  real(kind=8) :: bestgvns,fbest,fcw,fgvns,finish,fnn,seed,start,sumgvns,tcw,tlim,tnn,tsol
  type(instance_type) :: inst
  type(solution_type) :: sol,solbest,solcw,solcvrpnbest,solnn
  type(statistics_type) :: stat

  ! LOCAL ARRAYS
  character(len=80) :: filename
  character(len=256) :: arg_filename
  character(len=2) :: otype,ptype
  character(len=2) :: inisol

  ! ------------------------------------------------------------
  ! Read command line arguments
  ! ------------------------------------------------------------
  if ( command_argument_count() .ne. 1 ) then
     write(*,*) 'ERROR: Invalid number of command-line arguments.'
     write(*,*) 'Usage: ./master <instance_file.dat>'
     stop
  endif

  call get_command_argument(1,arg_filename)
  filename = trim(arg_filename)

  if (debug) write(*,*) 'Instance: ',trim(filename)

  !----------------------------------------------------
  ptype  = 'N'        
  recomp = .true.
  seed   = seedbase
  !-----------------------------------------------------
  !                   GVNS Parameters                  !
  !----------------------------------------------------- 
  otype  = 'F'         ! Type of operation: 'B' = Best Improvement, 'F' = First Improvement
  rmax   = 5           ! Maximum number of removed arcs or elements in the method
  norder = 4           ! Order of neighborhoods when applying the GVNS algorithm
  qmax   = 3           ! Maximum value for Q parameter (used in the method, e.g., shaking size)
  flip   = .false.     ! Boolean flag to activate the "flip" mechanism in the algorithm
  tlim   = 10.0d0      ! Time limit in seconds for solving each instance
  !------------------------------------------------------
  nrun   = 10          ! Number of independent GVNS algorithm executions
  !------------------------------------------------------
  
  call CreateInstance(inst,filename)  
  call CreateSolution(inst,solbest)
  call CreateSolution(inst,solcw)
  call CreateSolution(inst,solnn)
  call CreateSolution(inst,sol)
  call CreateSolution(inst,solcvrpnbest)
  call CreateStatistics(stat)

  write(*,*) 'Instance: ',trim(filename)
  write(*,*) 'Number of clients = ',inst%nclients
  write(*,*) 'Maximum vehicle capacity: ',int(inst%maxcap)
  !write(*,*) 'scalefactor = ',inst%scalefactor
  write(*,*)

  write(*,*) '--------------------------------------------------------------------------------------------'
  write(*,*) 'Solving the CVRPN via General Variable Neighborhood Search (GVNS) framework:'
  write(*,*) 'Initial solution: Best candidate from CW or NN constructive heuristics.'
  write(*,*) 'Number of independent runs = ',nrun
  write(*,*) 'Time limit per run (in seconds) = ',int(tlim)  
  write(*,*) '--------------------------------------------------------------------------------------------'
  write(*,*)
 
  !-------------------------------------------------------------------------------
  ! INITIAL SOLUTION: Best candidate between CW and NN
  !-------------------------------------------------------------------------------
  write(*,*) '----------------------------------------------------'
  write(*,*) '          CONSTRUCTIVE HEURISTICS PHASE             '
  write(*,*) '----------------------------------------------------'
  write(*,*)

  write(*,*) '------------------------------'
  ! Execute Clarke and Wright (CW) Savings heuristic
  call cpu_time(start)
  call CW(stat,inst,solcw,seed,ptype,-1.0d0)
  call cpu_time(finish)
  fcw = inst%scalefactor * sum( solcw%length(1:solcw%nroutes) )
  tcw = finish - start

  ! Execute Nearest Neighbor (NN) heuristic
  call cpu_time(start)
  call NN(stat,inst,solnn,seed,ptype,1,0.0d0)
  call cpu_time(finish)
  fnn = inst%scalefactor * sum( solnn%length(1:solnn%nroutes) )
  tnn = finish - start
  
  ! Display objective function values f for both constructive methods
  write(*,*) 'f(CW) = ',fcw
  write(*,*) 'f(NN) = ',fnn

  ! Select the best initial incumbent solution
  if ( fcw .le. fnn ) then
     call CopySolution(inst,solcw,solbest)
     inisol = 'CW'
  else
     call CopySolution(inst,solnn,solbest)
     inisol = 'NN'
  end if

  fbest = inst%scalefactor * sum( solbest%length(1:solbest%nroutes) )
  write(*,*) '------------------------------'
  write(*,*)

  !-------------------------------------------------------------------------------
  ! INITIAL SOLUTION DETAILED REPORT
  !-------------------------------------------------------------------------------
  write(*,*) '================  Initial Solution Report  ================'
  write(*,*) 'Heuristic: ',inisol
  write(*,*) 'Sum of all routes lengths = ',fbest
  write(*,*) 'Number of routes = ',count( solbest%nclients(1:solbest%nroutes) .gt. 0 )

  do i = 1,solbest%nroutes
     if ( solbest%nclients(i) .gt. 0 ) then
        write(*,*)
        write(*,*) 'RouteID = ',count( solbest%nclients(1:i) .gt. 0 )
        write(*,*) 'Route number of clients = ',solbest%nclients(i)
        write(*,*) 'Route length = ',inst%scalefactor * solbest%length(i)
        write(*,*) 'Route demand = ',int( solbest%aggdem(i) )
        write(*,*) 'Sequence = ',solbest%route(1:solbest%nclients(i),i)
        write(*,*) 'Visit points = ',inst%scalefactor * solbest%pathcoo(1:2,1:solbest%nclients(i),i)
     end if
  end do
  write(*,*)

  if ( inisol .eq. 'CW' ) then
     write(*,*) 'Constructive heuristic CPU Time (in seconds) = ',tcw
  else
     write(*,*) 'Constructive heuristic CPU Time (in seconds) = ',tnn
  end if

  write(*,*)

  seed = seedbase
  bestgvns = huge(1.0d0)
  sumgvns = 0.0d0
  !-------------------------------------------------------------------------------
  ! GVNS METAHEURISTIC PHASE
  !-------------------------------------------------------------------------------
  write(*,*) '----------------------------------------------------'
  write(*,*) '             GVNS METAHEURISTIC PHASE               '
  write(*,*) '----------------------------------------------------'
  write(*,*)

  do k = 1,nrun
     seed = k * seedbase    

     call CopySolution(inst,solbest,sol)   
     call CleanStatistics(stat)

     ! Execute the GVNS metaheuristic
     call GVNS(stat,inst,sol,seed,otype,ptype,recomp,rmax,norder,qmax,flip,tlim)

     ! Compute the objective function value obtained in this run
     fgvns = inst%scalefactor * sum( sol%length(1:sol%nroutes) )
     sumgvns = sumgvns + fgvns

     ! Results of the current run
     write(*,*) 'Run = ',k,' | Runtime (s) = ',stat%tsol,' | f(GVNS) = ',fgvns

     ! Save the best GVNS run
     if ( fgvns .lt. bestgvns - macheps12 * max( 1.0d0, abs( bestgvns ) ) ) then
        call CopySolution(inst,sol,solcvrpnbest)
        bestgvns = fgvns
        bestrun  = k
        ncycles  = stat%ncsol
        nshakes  = stat%nssol
        nneighs  = stat%nnsol
        tsol     = stat%tsol
     end if
  end do

  write(*,*)
  write(*,*) '----------------------------------------------------'
  write(*,*) '             GVNS EXECUTION CONCLUDED               '
  write(*,*) '----------------------------------------------------'
  write(*,*)

  !-------------------------------------------------------------------------------
  ! FINAL REPORT: BEST SOLUTION
  !-------------------------------------------------------------------------------
  write(*,*) '==========================================================='
  write(*,*) '            FINAL EXPERIMENTAL RESULTS (BEST RUN)          '
  write(*,*) '===========================================================' 
  write(*,*) 'Best run = ',bestrun
  write(*,*) 'Best f(GVNS) = ',bestgvns
  write(*,*) 'Mean f(GVNS) = ',sumgvns / nrun
  write(*,*) 'Number of routes = ',count( solcvrpnbest%nclients(1:solcvrpnbest%nroutes) .gt. 0 )
  ! Algorithmic metrics recorded at the moment the best solution was found
  write(*,*) 'GVNS cycles = ',ncycles
  write(*,*) 'Shaking steps = ',nshakes
  write(*,*) 'Neighborhood moves = ',nneighs
  write(*,*) 'CPU Time of best run (s) = ',tsol
  write(*,*) '==========================================================='
  write(*,*)
  write(*,*) 'DETAILED ROUTE CONFIGURATION (BEST RUN):'

  ! Print route details ONLY for the best run found
  do i = 1, solcvrpnbest%nroutes
     if ( solcvrpnbest%nclients(i) .gt. 0 ) then
        write(*,*)
        ! Dynamically compute the sequential Route ID bypassing empty slots
        write(*,*) 'RouteID = ',count( solcvrpnbest%nclients(1:i) .gt. 0 )
        write(*,*) 'Route number of clients = ',solcvrpnbest%nclients(i)
        write(*,*) 'Route length = ',inst%scalefactor * solcvrpnbest%length(i)
        write(*,*) 'Route demand = ',int( solcvrpnbest%aggdem(i) )
        write(*,*) 'Sequence = ',solcvrpnbest%route(1:solcvrpnbest%nclients(i),i)
        write(*,*) 'Visit points = ',inst%scalefactor * solcvrpnbest%pathcoo(1:2,1:solcvrpnbest%nclients(i),i)
     end if
  end do
  write(*,*)
  write(*,*) '==========================================================='
  write(*,*)

  call DeleteInstance(inst)
  call DeleteSolution(solbest)
  call DeleteSolution(solcw)  
  call DeleteSolution(solnn)
  call DeleteSolution(sol)
  call DeleteSolution(solcvrpnbest)

end program master
