!******************************************************************************
!
!                               main_cimi_sami.f90
!  Created on September 14, 2015 by Alex Glocer, Code 673, NASA GSFC.
!
!  Contact: Alex Glocer   at alex.glocer-1@nasa.gov,   301-286-9475.
!  Contact: Joe Huba at  huba@nrl.navy.mil
!  Contact: Mei-Ching Fok at mei-ching.h.fok@nasa.gov, 301-286-1083.
!
!******************************************************************************

program cimi_sami
  use ModCrcmGrid,    ONLY: iProcSAMI=>iProc,nProcSAMI=>nProc,iCommSAMI=>iComm
  use ModCrcmGrid,    ONLY: iProcCIMI=>iProc,nProcCIMI=>nProc,iCommCIMI=>iComm
  use ModCRCM,        ONLY: IsStandalone
  use ModMpi
  use ModCrcm,        ONLY: init_mod_crcm, Time
  use ModFieldTrace,  ONLY: init_mod_field_trace
  use ModImTime,      ONLY: TimeMax
  use ModCrcmRestart, ONLY: DtSaveRestart,crcm_write_restart
  use ModReadParam
  use ModPrerunField, ONLY: UsePrerun, read_prerun, read_prerun_IE, DtRead
  use ModIeCrcm,      ONLY: UseWeimer
  use CON_planet, ONLY: init_planet_const, set_planet_defaults
!  use ModPrerunField, ONLY: UsePrerun, read_prerun, read_prerun_IE
 
  implicit none
 
  !set some mpi parameters
  integer,parameter :: nProcCIMItmp = 2, nProcSAMItmp = 9
  integer :: iProcsCIMI_I(nProcCIMItmp),iProcsSAMI_I(nProcSAMItmp)
  integer :: iCommGlobal, nProcGloabl

  integer :: iError, iStep, iRank
  real    :: DtAdvance
  real    :: DtRestart = 300.0 ! currently set at 300s, should be read in later
  real    :: DtMax = 60.0      ! maximum timestep
  !---------------------------------------------------------------------------

  !****************************************************************************
  ! Initiallize MPI and get number of processors and rank of given processor
  !****************************************************************************
  
  !---------------------------------------------------------------------------
  call MPI_INIT(iError)
  iCommGlobal = MPI_COMM_WORLD
  
  call MPI_COMM_RANK(iCommGlobal,iProcGlobal,iError)
  call MPI_COMM_SIZE(iCommGlobal,nProcGlobal,iError)
  nProcCIMI = nProcCIMItmp
  nProcSAMI = nProcSAMItmp

  if (nProc /= nProcCIMI+nProcSAMI) &
       call con_stop('nProc not equal to nProcCIMI+nProcSAMI')

  ! create CIMI and SAMI proc lists and MPI groups
  do iRank = 0,nProcSAMI-1
     iProcsSAMI_I(iRank) = iRank
  enddo
  call mpi_group_incl(MPI_GROUP_WORLD,nProcSAMI,iProcsSAMI_I,SAMI_GROUP,iError)

  do iRank = nProcSAMI,nProcGlobal-1
     iProcsCIMI_I(iRank) = iRank
  enddo
  call mpi_group_incl(MPI_GROUP_WORLD,nProcCIMI,iProcsCIMI_I,CIMI_GROUP,iError)

  ! create CIMI and SAMI communicators from the proc groups
   call mpi_comm_create(MPI_COMM_WORLD,SAMI_GROUP,SAMI_COMM_WORLD,iError)
   call mpi_comm_create(MPI_COMM_WORLD,CIMI_GROUP,CIMI_COMM_WORLD,iError)
   
   iCommSAMI = SAMI_COMM_WORLD
   iCommCIMI = CIMI_COMM_WORLD

   ! assign sami and cimi rank
   call MPI_COMM_RANK(iCommSAMI,iProcSAMI,iError)
   call MPI_COMM_RANK(iCommCIMI,iProcCIMI,iError)

  !****************************************************************************
  ! Read the input file
  !****************************************************************************
  IsStandAlone=.true.

  ! Initial setup for the rbe model
  call read_file('PARAM.in',iCommCIMI)
  call read_init('  ',iSessionIn=1,iLineIn=0)
  call CRCM_set_parameters('READ')

  !if (usePrerun) call read_prerun(t)
  !if (usePrerun .and. iConvect==2) call read_prerun_IE(t)

  !\
  ! Initialize the planetary constant library and set Earth
  ! as the default planet.
  !/
  call init_planet_const
  call set_planet_defaults
  !****************************************************************************
  ! Initialize the model
  !****************************************************************************
  call init_mod_crcm
  call init_mod_field_trace
  
  ! Start Timing
  call timing_active(.true.)
  call timing_step(0)
  call timing_start('CRCM')

  !read initial prerun field
  if (UsePrerun) call read_prerun(Time)
  if (UsePrerun .and. .not.UseWeimer) call read_prerun_IE(Time)
  
  ! init model
  call timing_start('crcm_init')
  call crcm_init
  call timing_stop('crcm_init')
  
  if (iProc == 0)call timing_report_total
  call timing_reset('#all',3)

  !****************************************************************************
  ! start Timestepping
  !****************************************************************************
  iStep = 0
  TIMELOOP:do
     !report progress on proc 0
     if (iProc==0)write(*,*) 'In Time Loop iStep,Time = ', iStep,Time
     ! If Time exceeds max time then stop advancing
     if (Time >= TimeMax) exit TIMELOOP
     
     ! Set time to advance the model to either the time to reach TimeMax or 
     ! the next restart time
     DtAdvance = min(TimeMax - Time, DtMax)
  
     ! If DtAdvance is too small then just stop advancing
     if (DtAdvance < 1.0e-6) then
        Time = TimeMax
        exit TIMELOOP
     endif
     
     ! Call crcm_run to advance the Timestep
     call timing_step(iStep)
     call timing_start('crcm_run')     
     call crcm_run(DtAdvance)
     call timing_stop('crcm_run')
     
     ! Save restart at DtSaveRestart or TimeMax
     if (floor((Time+1.0e-5)/DtSaveRestart) /= &
          floor((Time+1.0e-5-DtAdvance)/DtSaveRestart)) then
        call crcm_write_restart
     endif
     
     ! Read new prerun file if using prerun fields and it is time to read
     if (floor((Time+1.0e-5)/DtRead) /= &
          floor((Time+1.0e-5-DtAdvance)/DtRead)) then 
        if (UsePrerun) call read_prerun(Time)
        if (UsePrerun .and. .not.UseWeimer) call read_prerun_IE(Time)     
     endif
     ! Advance the time iteration step
     iStep=iStep+1
  end do TIMELOOP

  ! Save restart at TimeMax
  call crcm_write_restart

  ! Finalize timing commands
  call timing_stop('CRCM')

  if (iProc == 0) then
     write(*,'(a)') 'Finished CRCM run, (reporting timings)'
     write(*,'(a)') '--------------------------------------'
     call timing_report
  endif
  
  ! Finalize MPI
  call MPI_FINALIZE(iError)

end program crcm
!============================================================================
subroutine CON_stop(StringError)
  use ModCrcmGrid,    ONLY: iProc,nProc
  use ModCrcm,        ONLY: Time
  use ModMpi
  implicit none
  character (len=*), intent(in) :: StringError
  
  ! Local variables:
  integer :: iError,nError
  !----------------------------------------------------------------------------
  
  write(*,*)'Stopping execution! me=',iProc,' at time=',Time,&
       ' with msg:'
  write(*,*)StringError
  call MPI_abort(MPI_COMM_WORLD, nError, iError)
  stop
end subroutine CON_stop
!============================================================================
subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe
  DoTest   = .false.
  DoTestMe = .false.
end subroutine CON_set_do_test
