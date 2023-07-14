program main_iono
  use ModGridIono,     only: init_grid_iono,iProc, nProc, iComm
  use GIMME_conductance , only: init_conductance, set_conductance, calc_coeficients
  use ModMagInput    , only: init_mag_input, test_input
  use ModGimmePlots  , only: gimme_plots
  use ModMpi
  use CON_planet, ONLY: init_planet_const, set_planet_defaults,is_planet_init
  
  character(len=100) :: NamePlanet='EARTH'
  integer :: errcode
  logical :: IsPlanetSet
  !-----------------------------------------------------------------------------
   !**************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !***************************************************************************

  call MPI_INIT(errcode)
  iComm = MPI_COMM_WORLD
  
  call MPI_COMM_RANK(iComm,iProc,errcode)
  call MPI_COMM_SIZE(iComm,nProc,errcode)
  
  !\         
  ! Initialize the planetary constant library and set Earth         
  ! as the default planet.
  !/
  call init_planet_const
  
  if (NamePlanet == 'EARTH') then
     call set_planet_defaults
     IsPlanetSet = .true.
  else
     write(*,*) NamePlanet
     IsPlanetSet = is_planet_init(NamePlanet)
  endif
  
  if (.not.IsPlanetSet) then
     call CON_stop('Planet not set. Stopping GIMME')
  endif

  
  
  call init_grid_iono
  call init_conductance
  call set_conductance
  call calc_coeficients
  call init_mag_input
  call test_input

  call gimme_potential
  
  call gimme_plots
end program main_iono

!============================================================================ 
! The following subroutines are here so that we can use SWMF library routines 
! Also some features available in SWMF mode only require empty subroutines 
! for compilation of the stand alone code.
!============================================================================
subroutine CON_stop(StringError)
  use ModGridIono, ONLY : iProc,iComm,TimeSimulation
  use ModMpi
  implicit none
  character (len=*), intent(in) :: StringError
  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------

  write(*,*)'Stopping execution! me=',iProc,' at time=',TimeSimulation,&
       ' with msg:'
  write(*,*)StringError
  call MPI_abort(iComm, nError, iError)
  stop
end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest = .false.; DoTestMe = .false.

end subroutine CON_set_do_test

subroutine CON_io_unit_new(iUnit)

  use ModIoUnit, ONLY: io_unit_new
  implicit none
  integer, intent(out) :: iUnit

  iUnit = io_unit_new()

end subroutine CON_io_unit_new
