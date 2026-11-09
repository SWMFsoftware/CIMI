!******************************************************************************
!
!                           integration__main.f90
!
! Integration test driver for cimi.
! Derived from cimi_main.f90 but stripped of MPI, timing, restart logic, with
! simple main time loop. Runs a single test point.
!
!******************************************************************************

program integration_test
  use ModMpi
  use ModCimiGrid,       ONLY: iProc, nProc, iComm
  use ModCimi,           ONLY: init_mod_cimi, Time, IsStandalone
  use ModReadParam,      ONLY: read_file,read_init
  use CIMI_SetParams,    ONLY: CIMI_set_parameters
  use ModCimiPlanet,     ONLY: init_cimi_planet_const
  use CON_planet,        ONLY: init_planet_const, set_planet_defaults
  use ModCimiTrace,      ONLY: init_mod_field_trace
  use ModCimiMethods,    ONLY: cimi_init
  use CIMI_DiffusionTests, ONLY: diffusion_test_run
 
  implicit none
 
  integer :: iError, iStep
  real :: DtAdvance

  !---------------------------------------------------------------------------

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************
  call MPI_INIT(iError)
  iComm = MPI_COMM_WORLD
  call MPI_COMM_RANK(iComm,iProc,iError)
  call MPI_COMM_SIZE(iComm,nProc,iError)
  
  !****************************************************************************
  ! Read the input file
  !****************************************************************************
  IsStandAlone = .true.

  ! Initial setup for the model
  call read_file('PARAM.in', iComm)
  call read_init('  ', iSessionIn=1, iLineIn=0)
  call CIMI_set_parameters('READ')

  ! Initialize the planetary constant library and set Earth as default
  call init_planet_const
  call set_planet_defaults

  !****************************************************************************
  ! Initialize the model
  !****************************************************************************
  call init_mod_cimi
  call init_mod_field_trace
  call init_cimi_planet_const

  ! Start Timing
  call timing_active(.true.)
  call timing_step(0)
  call timing_start('CIMI')

  ! General model initialization
  call timing_start('cimi_init')
  call cimi_init
  call timing_stop('cimi_init')
  
  write(*,*) 'Initialization complete. Starting diffusion integration test...'

  !****************************************************************************
  ! Single Test Point Execution
  !****************************************************************************
  
  ! Hardcode timestep for single step integration test
  DtAdvance = 5.0
  
  ! Run the specific diffusion test routine
  do iStep=0,25
    if (mod(Time,3600.0).eq.0) then
        write(*,'(a,i4,i8,f10.2)') 'IM: In Time Loop iProc, iStep, Time = ',&
                                                     iProc, iStep, Time
    endif
    call diffusion_test_run(DtAdvance)
  end do

  ! Finalize timing commands
  call timing_stop('CIMI')

  write(*,'(a)') 'Finished DiffusionTest run, (reporting timings)'
  write(*,'(a)') '-----------------------------------------------'
  call timing_report

  ! Finalize MPI
  call MPI_FINALIZE(iError)

end program integration_test

