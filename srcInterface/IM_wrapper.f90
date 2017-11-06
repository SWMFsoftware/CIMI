!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module IM_wrapper

  ! Wrapper for CIMI Internal Magnetosphere (IM) component

  implicit none

  private ! except

  public:: IM_set_param
  public:: IM_init_session
  public:: IM_run
  public:: IM_save_restart
  public:: IM_finalize

  ! Coupling with IE
  public:: IM_get_for_ie
  public:: IM_put_from_ie_mpi
  public:: IM_put_from_ie
  public:: IM_put_from_ie_complete

  ! Coupling with GM
  public:: IM_get_for_gm
  public:: IM_put_from_gm
  public:: IM_put_from_gm_line
  public:: IM_put_from_gm_crcm
  public:: IM_put_sat_from_gm

contains
  !===========================================================================
  subroutine IM_set_param(CompInfo, TypeAction)

    !USES:
    use CON_comp_info
    use ModUtilities
    use ModReadParam
    use CON_coupler, ONLY: Couple_CC, GM_, IM_, IE_
    use ModGmCimi,  ONLY: UseGm
    use ModIeCimi,  ONLY: UseIE
    use ModCimiGrid, ONLY: iProc,nProc,iComm
    implicit none

    character (len=*), intent(in)     :: TypeAction ! which action to perform
    type(CompInfoType), intent(inout) :: CompInfo   ! component information

    character (len=*), parameter :: NameSub='IM_set_param'
    integer :: iError
    character (len=100) :: NameCommand
    logical             :: UseStrict=.true.

    !------------------------------------------------------------------------
    UseGm = Couple_CC(GM_, IM_) % DoThis
    UseIE = Couple_CC(IE_, IM_) % DoThis

    !if(iProc>=0)then
    !   call IM_write_prefix;  
    !   write(iUnitOut,*) NameSub,' TypeAction= ',TypeAction, &
    !        '    iProc=',iProc
    !end if
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,                                     &
            Use=.true.,                                       &
            NameVersion='CIMI, M. Fok and A. Glocer', &
            Version=1.0)
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
       !if(nProc>1)call CON_stop('IM_ERROR this version can run on 1 PE only!')
    case('STDOUT')
       !iUnitOut=STDOUT_
       !StringPrefix='RB:'
    case('FILEOUT')
       !call get(CompInfo,iUnitOut=iUnitOut)
       !StringPrefix=''
    case('READ')
       call CIMI_set_parameters('READ')
    case('CHECK')
       ! We should check and correct parameters here
       !if(iProc==0)then
       !   call IM_write_prefix;  write(iUnitOut,*)&
       !        NameSub,': CHECK iSession =',i_session_read()
       !end if
    case('GRID')
       call IM_set_grid
    case default
       call CON_stop(NameSub//' IM_ERROR: invalid TypeAction='//TypeAction)
    end select

  end subroutine IM_set_param
  !============================================================================
  subroutine IM_set_grid

    use ModCimiGrid,      ONLY: np, nt, xlat, phi
    use ModNumConst,      ONLY: cDegToRad
    use CON_coupler,      ONLY: set_grid_descriptor, is_proc, IM_

    implicit none

    character (len=*), parameter :: NameSub='IM_set_grid'
    integer :: iSize,jSize
    real, dimension (:,:), allocatable :: gridLat,gridLT
    real :: Radius_I(1)
    logical :: IsInitialized=.false.
    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest)write(*,*)'IM_set_grid_descriptor called, IsInitialized=',&
         IsInitialized
    if(IsInitialized) return
    IsInitialized=.true.

    ! radial size of the ionosphere in meters
    Radius_I(1) = (6375.0+120.0)*1000.0 

    if(DoTest)then
       write(*,*)NameSub,' ir,ip=',np,nt
       write(*,*)NameSub,' size(xlat)=',size(xlat),' size(phi)=',size(phi)
       write(*,*)NameSub,' xlat      =',xlat
       write(*,*)NameSub,' phi       =',phi
       write(*,*)NameSub,' rIono     =',Radius_I(1)
    end if

    ! RB grid size in generalized coordinates
    call set_grid_descriptor( IM_,                 & ! component index
         nDim=2,                                   & ! dimensionality
         nRootBlock_D=(/1,1/),                     & ! single block
         nCell_D=(/np, nt/),                       & ! size of cell based grid
         XyzMin_D=(/0.5,    0.5/),                 & ! min gen.coords for cells
         XyzMax_D=(/np+0.5, nt+0.5/),              & ! max gen.coords for cells
         TypeCoord='SMG',                          & ! solar magnetic coord
         Coord1_I=(90.0-xlat(1:np))*cDegToRad,     & ! colatitude in radians
         Coord2_I=phi,                             & ! longitude in radians
         Coord3_I=Radius_I,                        & ! radial size in meters
         IsPeriodic_D=(/.false.,.true./),          & ! periodic in longitude
         nVar = 7,                                 & ! number of "fluid" vars
         NameVar = 'p rho ppar Hpp Opp Hprho Oprho') ! names of "fluid" vars

  end subroutine IM_set_grid

  !============================================================================

  subroutine IM_init_session(iSession, TimeSimulation)

    use ModCimiGrid,      ONLY: np, xlat
    use ModCimi,          ONLY: init_mod_cimi, Time
    use ModFieldTrace,    ONLY: init_mod_field_trace
    use ModImTime
    use ModTimeConvert,   ONLY: time_real_to_int
    use CON_physics,      ONLY: get_time

    implicit none

    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IM_init_session'
    !------------------------------------------------------------------------
    ! GM info needed before initialization just set up latitude/longitude grid
    call init_mod_cimi
    call init_mod_field_trace
    call cimi_init

    call get_time(tStartOut = StartTime)
    Time = TimeSimulation
    CurrentTime = StartTime + Time  
    call time_real_to_int(StartTime, iStartTime_I)

  end subroutine IM_init_session

  !============================================================================

  subroutine IM_run(TimeSimulation,TimeSimulationLimit)

    use CON_time,   ONLY: DoTimeAccurate
    use ModCimi,    ONLY: Time, dt, dtmax,iProc

    implicit none

    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded
    real, intent(inout):: TimeSimulation   ! current time of component

    Logical, save :: IsInitiallized = .false.
    character(len=*), parameter :: NameSub='IM_run'

    !------------------------------------------------------------------------

    if (.not. IsInitiallized) then
       call cimi_init
       IsInitiallized = .true.
    endif

    if( .not. DoTimeAccurate)then
       ! steady state mode
       dt = dtmax
       call cimi_run(2*dt)
    else
       ! time accurate mode
       dt = min(dtmax, 0.5*(TimeSimulationLimit - TimeSimulation))
       call cimi_run(TimeSimulationLimit - TimeSimulation)
       ! return time at the end of the time step to CON
       TimeSimulation   = TimeSimulationLimit
    end if

  end subroutine IM_run
  !===========================================================================

  subroutine IM_finalize(TimeSimulation)

    !USES:
    implicit none

    real,     intent(in) :: TimeSimulation   ! seconds from start time
    character(len=*), parameter :: NameSub='IM_finalize'

    !-------------------------------------------------------------------------

    !call IM_write_prefix; write(iUnitOut,*) &
    !     NameSub,' at TimeSimulation=',TimeSimulation

  end subroutine IM_finalize
  !===========================================================================

  subroutine IM_save_restart(TimeSimulation)
    use ModCimiRestart, ONLY: cimi_write_restart
    implicit none

    real,     intent(in) :: TimeSimulation   ! seconds from start time
    character(len=*), parameter :: NameSub='IM_save_restart'

    !-------------------------------------------------------------------------
    call cimi_write_restart

  end subroutine IM_save_restart
  !===========================================================================

  subroutine IM_put_from_gm_crcm(Buffer_IIV,iSizeIn,jSizeIn,nVarIn,&
       BufferLine_VI,nVarLine,nPointLine,NameVar,tSimulation)
    use ModGmCIMI
    use ModCimiGrid,  ONLY: nLat => np, nLon => nt, iProc, nProc
    use ModCimiPlanet,ONLY: rEarth => re_m
    use ModTsyInput,  ONLY: xnswa,vswa,bxw,byw,bzw,nsw,iyear,iday,UseSmooth
    use ModPrerunField,ONLY: DoWritePrerun, save_prerun
    use ModIoUnit, ONLY: UnitTmp_
    !  use ModPrerunField,ONLY: DoWritePrerun, save_prerun
    implicit none

    integer, intent(in) :: iSizeIn, jSizeIn, nVarIn
    real,    intent(in) :: Buffer_IIV(iSizeIn,jSizeIn,nVarIn)
    integer, intent(in) :: nVarLine, nPointLine
    real,    intent(in) :: BufferLine_VI(nVarLine, nPointLine)

    character (len=*),intent(in) :: NameVar
    real, intent(in) :: tSimulation

    real, parameter :: noValue=-99999.
    real :: SwDensMax, SwVelMax, SwDensMin, SwVelMin
    integer :: n,iLat,iLon
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='IM_put_from_gm_crcm'
    logical,save :: IsFirstCall = .true.
    character(len=100) :: NameOut
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    DoMultiFluidGMCoupling = .false.
    DoAnisoPressureGMCoupling = .false.

    if(NameVar == 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p:Hprho:Oprho:Hpp:Opp')then
       DoMultiFluidGMCoupling = .true.
    elseif(NameVar == 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p:ppar')then
       DoAnisoPressureGMCoupling = .true.
    elseif(NameVar /= 'x:y:bmin:I_I:S_I:R_I:B_I:rho:p')then
       call CON_stop(NameSub//' invalid NameVar='//NameVar)
    end if

    if(nVarLine /= nVar) then
       write(*,*)'nVarLine=',nVarLine
       call CON_stop(NameSub//' invalid nVarLine (should be 4)')
    end if
    DoneGmCoupling=.true.
    if(DoTestMe)then
       write(*,*)NameSub,' iSizeIn,jSizeIn,nVarIn=',&
            iSizeIn,jSizeIn,nVarIn
       write(*,*)NameSub,' nVarLine,nPointLine=',nVarLine,nPointLine
       !write(*,*)NameSub,' Buffer_IIV(21,1,:)=',Buffer_IIV(21,1,:)
       write(*,*)NameSub,' BufferLine_VI(:,1) =',BufferLine_VI(:,1)
       write(*,*)NameSub,' BufferLine_VI(:,2) =',BufferLine_VI(:,2)
       write(*,*)NameSub,' IMF: Density  = ',Buffer_IIV(1,1,6)
       write(*,*)NameSub,' IMF: Velocity = ',Buffer_IIV(2,1,6)
       write(*,*)NameSub,' IMF: Bx       = ',Buffer_IIV(5,1,6)
       write(*,*)NameSub,' IMF: By       = ',Buffer_IIV(6,1,6)
       write(*,*)NameSub,' IMF: Bz       = ',Buffer_IIV(7,1,6)
    end if

    if (allocated(StateLine_VI)) then
       deallocate(StateLine_VI,StateBmin_IIV)
    endif

    if (.not.allocated(StateLine_VI)) then
       allocate(StateLine_VI(nVarLine,nPointLine),&
            StateBmin_IIV(iSizeIn,jSizeIn,nVarIn))
    endif

    StateLine_VI      = BufferLine_VI
    StateBmin_IIV(:,:,:) = Buffer_IIV(:,:,:)

    ! convert unit of locations X and Y
    StateBmin_IIV(:,:,1:2) = StateBmin_IIV(:,:,1:2)/rEarth ! m --> Earth Radii

    nPoint    = nPointLine
    nVarBmin = nVarIn
    !Convert Units
    StateLine_VI(2,:) = StateLine_VI(2,:) / rEarth ! m --> Earth Radii
    StateLine_VI(3,:) = StateLine_VI(3,:) / rEarth ! m --> Earth Radii

    ! for anisopressure coupling 
    if(DoTest .and. DoAnisoPressureGMCoupling)then
       write(NameOut,"(a,f6.1)") 'StateBmin_t_',tSimulation
       open(UnitTmp_,FILE=NameOut)
       write(UnitTmp_,"(a)") &
            'IM_put_from_gm_crcm, StateBmin_IIV, last index 1:5 and 7 '
       write(UnitTmp_,"(a)") &
            'Xeq Yeq Beq rho p ppar'
       do iLat =iSizeIn,1,-1
          do iLon =1,jSizeIn
             write(UnitTmp_,"(100es18.10)")StateBmin_IIV(iLat,iLon,1:5), &
                  StateBmin_IIV(iLat,iLon,7)
          enddo
       enddo
       close(UnitTmp_)
    end if

    !  if (nProc == 1) then
    !     write(NameOut,"(a,i2.2)") 'StateBmin1_0_iproc_',iProc
    !     open(UnitTmp_,FILE=NameOut)
    !     do iLat =1,nLat
    !        write(UnitTmp_,"(100es18.10)")StateBmin_IIV(iLat,1:nLon,1)
    !     enddo
    !     close(UnitTmp_)
    !     write(*,*)'!!! iProc, StateBmin_IIV(21,25)',iProc, StateBmin_IIV(21,25,1)
    !     write(*,*) '!!! iSizeIn,jSizeIn',iSizeIn,jSizeIn
    !     call con_stop('')
    !  else
    !     write(NameOut,"(a,i2.2)") 'StateBmin1_iproc_',iProc
    !     open(UnitTmp_,FILE=NameOut)
    !     do iLat =1,nLat
    !        write(UnitTmp_,"(100es18.10)")StateBmin_IIV(iLat,1:nLon,1)
    !     enddo
    !     close(UnitTmp_)
    !     write(*,*)'!!! iProc, StateBmin_IIV(21,25)',iProc, StateBmin_IIV(21,25,1)
    !     write(*,*) '!!! iSizeIn,jSizeIn',iSizeIn,jSizeIn
    !     call con_stop('')
    !  endif


    !Solar wind values
    if(IsFirstCall .or. (.not. UseSmooth)) then
       xnswa(1) = Buffer_IIV(1,1,6)*1.0e-6                   !m^-3 -->/cc
       vswa (1) = sqrt(sum(Buffer_IIV(2:4,1,6)**2.0))*1.0e-3 !m/s-->km/s
    else
       ! Update Solar wind value, but do not let them change 
       ! more than 5 percent per update
       SwDensMax = 1.05*xnswa(1)
       SwDensMin = 0.95*xnswa(1)
       SwVelMax  = 1.05*vswa(1)
       SwVelMin  = 0.95*vswa(1)
       xnswa(1) = min(SwDensMax,Buffer_IIV(1,1,6)*1.0e-6)
       xnswa(1) = max(SwDensMin,xnswa(1))
       vswa(1)  = min(SwVelMax,sqrt(sum(Buffer_IIV(2:4,1,6)**2.0))*1.0e-3)
       vswa(1)  = max(SwVelMin,vswa(1))
    endif
    bxw(1) = Buffer_IIV(5,1,6)*1.0e9      !T --> nT
    byw(1) = Buffer_IIV(6,1,6)*1.0e9      !T --> nT
    bzw(1) = Buffer_IIV(7,1,6)*1.0e9      !T --> nT

    nsw = 1

    iyear=2002
    iday=1


    ! create an index array on the first call
    if (IsFirstCall) then
       n = 0
       do iLon = 1, nLon
          do iLat = 1, nLat
             n = n+1
             iLineIndex_II(iLon,iLat) = n
          end do
       end do
       IsFirstCall = .false.
    endif

    if (iProc==0 .and. DoWritePrerun) call save_prerun(tSimulation)
  end subroutine IM_put_from_gm_crcm

  !============================================================================
  subroutine IM_put_from_gm_line(nRadiusIn, nLonIn, Map_DSII, &
       nVarLineIn, nPointLineIn, BufferLine_VI, NameVar)

    integer, intent(in) :: nRadiusIn, nLonIn
    real,    intent(in) :: Map_DSII(3,2,nRadiusIn,nLonIn)
    integer, intent(in) :: nVarLineIn, nPointLineIn
    real,    intent(in) :: BufferLine_VI(nVarLineIn,nPointLineIn)
    character(len=*), intent(in) :: NameVar

    character (len=*),parameter :: NameSub='IM_put_from_gm_line'

    call CON_stop(NameSub//' should not be called for IM/CIMI')

  end subroutine IM_put_from_gm_line
  !============================================================================
  subroutine IM_put_from_gm(Buffer_IIV,iSizeIn,jSizeIn,nVarIn,NameVar)

    integer, intent(in) :: iSizeIn,jSizeIn,nVarIn
    real, dimension(iSizeIn,jSizeIn,nVarIn), intent(in) :: Buffer_IIV
    character (len=*),intent(in)       :: NameVar

    character (len=*),parameter :: NameSub = 'IM_put_from_gm'

    call CON_stop(NameSub//' should not be called for IM/CIMI')

  end subroutine IM_put_from_gm
  !============================================================================
  subroutine IM_put_from_ie_mpi(nTheta, nPhi, Potential_II)

    integer, intent(in):: nTheta, nPhi
    real,    intent(in):: Potential_II(nTheta, nPhi)

    character(len=*), parameter   :: NameSub='IM_put_from_ie_mpi'

    call CON_stop(NameSub//' cannot be used by CIMI!')

  end subroutine IM_put_from_ie_mpi

  !============================================================================
  subroutine IM_put_sat_from_gm(nSats, Buffer_I, Buffer_III)

    ! Puts satellite locations and names from GM into IM.
    use ModImSat, ONLY: nImSats, DoWriteSats, ReadRestartSat, &
         NameSat_I, SatLoc_3I
    use ModNumConst,   ONLY: cDegToRad
    use ModCimiGrid,   ONLY: iProc
    use ModIoUnit,     ONLY: UnitTmp_

    implicit none
    character (len=*),parameter :: NameSub='IM_put_sat_from_gm'

    ! Arguments
    integer, intent(in)            :: nSats
    real, intent(in)               :: Buffer_III(4,2,nSats)
    character(len=100), intent(in) :: Buffer_I(nSats)

    ! Internal variables
    integer :: iError, iSat, l1, l2, iRow
    !--------------------------------------------------------------------------
    ! Activate satellite writing in RCM
    DoWriteSats = .true.
    nImSats = nSats

    ! Check allocation of sat tracing variables
    if(allocated(SatLoc_3I)) deallocate(SatLoc_3I)
    if(allocated(NameSat_I)) deallocate(NameSat_I)

    allocate(SatLoc_3I(4,2,nImSats), stat=iError)
    allocate(NameSat_I(nImSats),     stat=iError)

    if (iProc == 0 .and. ReadRestartSat) then
       open(UnitTmp_,file="IM/restartIN/restart.sat",&
            status="old", form="unformatted")

       read(UnitTmp_) nImSats
       
       do iSat=1,nImSats
          
          read(UnitTmp_) NameSat_I(iSat)
          
       end do

       do iSat=1,nImSats
          
          do iRow=1,2
             
             read(UnitTmp_) SatLoc_3I(1:4,iRow,iSat)
             
          end do
          
       end do
       
       close(UnitTmp_)

       ReadRestartSat = .false.
       
    else
       
       ! Assign incoming values, remove path and extension from name.
       SatLoc_3I = Buffer_III
       SatLoc_3I(4,2,:)=SatLoc_3I(4,2,:)
       do iSat=1, nSats
          l1 = index(Buffer_I(iSat), '/', back=.true.) + 1
          l2 = index(Buffer_I(iSat), '.') - 1
          if (l1-1<=0) l1=1
          if (l2+1<=0) l2=len_trim(Buffer_I(iSat))
          NameSat_I(iSat) = Buffer_I(iSat)(l1:l2)
       end do

    endif
 
  end subroutine IM_put_sat_from_gm

  !============================================================================

  subroutine IM_get_for_gm(Buffer_IIV,iSizeIn,jSizeIn,nVar,NameVar)

    !use CON_time, ONLY : get_time
    use ModCimiGrid,  ONLY: iSize=>np, jSize=>nt, iProc
    use ModCimi,      ONLY: Pressure_IC, PressurePar_IC, Bmin_C, Time, Pmin
    use ModGmCimi,    ONLY: Den_IC, iLatMin, DoMultiFluidGMCoupling, &
         DoAnisoPressureGMCoupling
    use ModFieldTrace,ONLY: iba
    use ModCimiPlanet,ONLY: nspec,amu_I
    use ModNumConst,  ONLY: cRadToDeg
    use ModConst,     ONLY: cProtonMass
    use ModIoUnit, ONLY: UnitTmp_

    implicit none
    character (len=*),parameter :: NameSub='IM_get_for_gm'

    integer, intent(in)                                :: iSizeIn,jSizeIn,nVar
    real, dimension(iSizeIn,jSizeIn,nVar), intent(out) :: Buffer_IIV
    character (len=*),intent(in)                       :: NameVar

    !LOCAL VARIABLES:
    real :: tSimulation
    integer, parameter :: pres_=1, dens_=2, parpres_=3, bmin_=4,&
         Hpres_=3, Opres_=4, Hdens_=5, Odens_=6

    integer :: i,j,k
    logical :: DoTest, DoTestMe
    character(len=100) :: NameOut
    !--------------------------------------------------------------------------

    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if (DoTestMe) &
         write(*,*)NameSub,' starting with iSizeIn,jSizeIn,nVar,NameVar=',&
         iSizeIn,jSizeIn,nVar,NameVar

    DoMultiFluidGMCoupling = .false.
    DoAnisoPressureGMCoupling = .false.

    if(NameVar == 'p:rho:Hpp:Opp:Hprho:Oprho')then
       DoMultiFluidGMCoupling = .true.
    elseif(NameVar == 'p:rho:ppar:bmin')then
       DoAnisoPressureGMCoupling = .true.
    elseif(NameVar /= 'p:rho')then
       call CON_stop(NameSub//' invalid NameVar='//NameVar)
    end if

    if(iSizeIn /= iSize .or. jSizeIn /= jSize)then
       write(*,*)NameSub//' incorrect buffer size=',iSizeIn,jSizeIn
       call CON_stop(NameSub//' SWMF_ERROR')
    end if

    Buffer_IIV = 0.

    !Fill pressure and density
    do i=1,iSize; do j=1,jSize
       if( i<iLatMin .or.  i > iba(j) ) then
          Buffer_IIV(i,j,pres_) = -1.
          Buffer_IIV(i,j,dens_) = -1.
          if(DoAnisoPressureGMCoupling)then
             Buffer_IIV(i,j,parpres_) = -1.
             Buffer_IIV(i,j,bmin_) = -1.
          end if
          if(DoMultiFluidGMCoupling)then
             !  Multifluid case
             Buffer_IIV(i,j,Hpres_) = -1.
             Buffer_IIV(i,j,Opres_) = -1.
             Buffer_IIV(i,j,Hdens_) = -1.
             Buffer_IIV(i,j,Odens_) = -1.
          end if
       else
          ! make sure pressure passed to GM is not lower than Pmin [nPa]
          ! to avoid too low GM pressure 
          Buffer_IIV(i,j,pres_) = max(sum(Pressure_IC(:,i,j)), Pmin)*1e-9
          Buffer_IIV(i,j,dens_) = &
               sum(Den_IC (1:nspec-1,i,j)*cProtonMass*amu_I(1:nspec-1))
          if(DoAnisoPressureGMCoupling)then
             Buffer_IIV(i,j,parpres_) = &
                  max(sum(PressurePar_IC(:,i,j)), Pmin)*1e-9
             ! fill minimum B  
             Buffer_IIV(i,j,bmin_) = Bmin_C(i,j)
          end if
          if(DoMultiFluidGMCoupling)then
             Buffer_IIV(i,j,Hpres_) = max(Pressure_IC(1,i,j), Pmin)*1e-9
             Buffer_IIV(i,j,Hdens_) = &
                  Den_IC (1,i,j)*cProtonMass*amu_I(1)
             Buffer_IIV(i,j,Opres_) = max(Pressure_IC(2,i,j), Pmin)*1e-9
             Buffer_IIV(i,j,Odens_) = &
                  Den_IC (2,i,j)*cProtonMass*amu_I(2)

          end if
       end if

       ! Only a not-a-number can be less than zero and larger than one
       if(  .not. Buffer_IIV(i,j,pres_) > 0 .and. &
            .not. Buffer_IIV(i,j,pres_) < 1) then
          write(*,*)NameSub,': ERROR IN PRESSURE'
          write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,pres_)
          call CON_stop(NameSub//': Not a number found in IM pressure!')
       end if
       if(  .not. Buffer_IIV(i,j,dens_) > 0 .and. &
            .not. Buffer_IIV(i,j,dens_) < 1) then
          write(*,*)NameSub,': ERROR IN DENSITY'
          write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,dens_)
          call CON_stop(NameSub//': Not a number found in IM density!')
       end if

       !multi-fluid
       if(DoMultiFluidGMCoupling)then
          if(  .not. Buffer_IIV(i,j,Hpres_) > 0 .and. &
               .not. Buffer_IIV(i,j,Hpres_) < 1) then
             write(*,*)NameSub,': ERROR IN PRESSURE'
             write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,Hpres_)
             call CON_stop(NameSub//': Not a number found in IM Hp pressure!')
          end if
          if(  .not. Buffer_IIV(i,j,Hdens_) > 0 .and. &
               .not. Buffer_IIV(i,j,Hdens_) < 1) then
             write(*,*)NameSub,': ERROR IN DENSITY'
             write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,Hdens_)
             call CON_stop(NameSub//': Not a number found in IM Hp density!')
          end if
          if(  .not. Buffer_IIV(i,j,Opres_) > 0 .and. &
               .not. Buffer_IIV(i,j,Opres_) < 1) then
             write(*,*)NameSub,': ERROR IN PRESSURE'
             write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,Opres_)
             call CON_stop(NameSub//': Not a number found in IM Op pressure!')
          end if
          if(  .not. Buffer_IIV(i,j,Odens_) > 0 .and. &
               .not. Buffer_IIV(i,j,Odens_) < 1) then
             write(*,*)NameSub,': ERROR IN DENSITY'
             call CON_stop(NameSub//': Not a number found in IM Op density!')
          end if
       end if

       !aniso pressure
       if(DoAnisoPressureGMCoupling)then
          if(  .not. Buffer_IIV(i,j,parpres_) > 0 .and. &
               .not. Buffer_IIV(i,j,parpres_) < 1) then
             write(*,*)NameSub,': ERROR IN PRESSURE'
             write(*,*)NameSub,': i,j,Buffer =',i,j,Buffer_IIV(i,j,parpres_)
             call CON_stop(NameSub// &
                  ': Not a number found in IM parallel pressure !')
          end if
       endif
    end do; end do

    !for anisopressure IM coupling
    if(DoTest .and. iProc == 0) then
       write(NameOut,"(a,f6.1,a)") 'IM_Buffer_IIV_t', Time, '.out'
       open(UnitTmp_,FILE=NameOut)
       write(UnitTmp_,"(a)") 'IM_get_for_GM, Buffer_IIV:'
       !    write(UnitTmp_,"(a)") ' '
       !    write(UnitTmp_,"(a)") ' '
       !    write(UnitTmp_,"(a)") ' '
       write(UnitTmp_,"(a)") 'bmin rho p ppar'
       do i = iSizeIn, 1, -1                                                                       
          do j =1,jSizeIn 
             write(UnitTmp_,"(4es18.10)") &
                  Buffer_IIV(i,j,bmin_), Buffer_IIV(i,j,dens_), &
                  Buffer_IIV(i,j,pres_), Buffer_IIV(i,j,parpres_)
          end do
       end do
       close(UnitTmp_)
    end if

    ! Units of rcm_mass_density are kg/m3
    !  where(Buffer_IIV(:,:,dens_) > 0.0) &
    !       Buffer_IIV(:,:,dens_) = Buffer_IIV(:,:,dens_)

    !  if(DoMultiFluidGMCoupling)then
!!!! ADD MultiFluid Coupling Later !!!
    !     ! MultiFluid                                                                  
    !     where(Buffer_IIV(:,:,Hpres_) > 0.0) &
    !          Buffer_IIV(:,:,Hpres_) = Buffer_IIV(:,:,Hpres_) * 1.67E-35
    !     where(Buffer_IIV(:,:,Opres_) > 0.0) &
    !          Buffer_IIV(:,:,Opres_) = Buffer_IIV(:,:,Opres_) * 1.67E-35
    !     ! Units of rcm_mass_density are kg/m3                                                
    !     where(Buffer_IIV(:,:,Hdens_) > 0.0) &
    !          Buffer_IIV(:,:,Hdens_) = Buffer_IIV(:,:,Hdens_) / 6.37E+15
    !     where(Buffer_IIV(:,:,Odens_) > 0.0) &
    !          Buffer_IIV(:,:,Odens_) = Buffer_IIV(:,:,Odens_) / 6.37E+15
    !  endif

    if(DoTestMe)write(*,*) NameSub,' finished'

  end subroutine IM_get_for_gm

  !============================================================================

  subroutine IM_put_from_ie(nPoint,iPointStart,Index,Weight,DoAdd,Buff_V,nVar)

    use ModCimiGrid,  ONLY: np, nt
    use CON_router,   ONLY: IndexPtrType, WeightPtrType
    use ModIeCimi,    ONLY: Pot!, birk_mhd, iSize, jSize, sigmaH_mhd,sigmaP_mhd
    
    implicit none
    character(len=*), parameter   :: NameSub='IM_put_from_ie'
    integer,intent(in)            :: nPoint, iPointStart, nVar
    real, intent(in)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight
    logical,intent(in)            :: DoAdd
    integer :: iBlock,i,j
    !--------------------------------------------------------------------------
    if(nPoint>1)then
       write(*,*)NameSub,': nPoint,iPointStart,Weight=',&
            nPoint,iPointStart,Weight % Weight_I
       call CON_stop(NameSub//': should be called with 1 point')
    end if
    if(DoAdd)then
       write(*,*)NameSub,': nPoint,iPointStart,Weight=',&
            nPoint,iPointStart,Weight % Weight_I
       write(*,*)NameSub,': WARNING DoAdd is true'
    end if

    i = Index % iCB_II(1,iPointStart)
    j = Index % iCB_II(2,iPointStart)

    if(i<1.or.i>np.or.j<1.or.j>nt)then
       write(*,*)'i,j,DoAdd=',i,j,DoAdd
       call CON_stop('IM_put_from_ie: index out of range')
    end if

    if(DoAdd)then
       Pot(i,j)        = pot(i,j)        + Buff_V(1)
       !     birk_mhd(i,j) = birk_mhd(i,j) + Buff_V(2)
       !     sigmaH_mhd(i,j) = sigmaH_mhd(i,j) + Buff_V(3)
       !     sigmaP_mhd(i,j) = sigmaP_mhd(i,j) + Buff_V(4)
    else
       pot(i,j)        = Buff_V(1)
       !     birk_mhd(i,j) = Buff_V(2)
       !     sigmaH_mhd(i,j) = Buff_V(3)
       !     sigmaP_mhd(i,j) = Buff_V(4)
    end if
    
    
    
  end subroutine IM_put_from_ie

  !============================================================================

  subroutine IM_put_from_ie_complete
    use ModCimiGrid,  ONLY: np,nt,iComm, iProc
    use ModIeCimi,    ONLY: Pot
    use ModPrerunField,ONLY: DoWritePrerun, save_prerun_IE
    use ModCimi,    ONLY: Time

    use ModMpi
    implicit none
    integer:: iError,nsize
    !--------------------------------------------------------------------------
    ! Bcast  potential to all procs
    nsize=np*nt
    call MPI_bcast(pot,nsize,MPI_REAL,0,iComm,iError)
    
    if (iProc == 0 .and. DoWritePrerun) call save_prerun_IE(Time)
    return
  end subroutine IM_put_from_ie_complete

  !============================================================================

  subroutine IM_get_for_ie(nPoint,iPointStart,Index,Weight,Buff_V,nVar)

    ! Provide current for IE
    ! The value should be interpolated from nPoints with
    ! indexes stored in Index and weights stored in Weight
    ! The variables should be put into Buff_V

    use ModCimi,      ONLY: FAC_C, nLat=>np, nLon=>nt
    use CON_router,   ONLY: IndexPtrType, WeightPtrType
    implicit none
    character(len=*), parameter :: NameSub='IM_get_for_ie'

    integer,intent(in)            :: nPoint, iPointStart, nVar
    real,intent(out)              :: Buff_V(nVar)
    type(IndexPtrType),intent(in) :: Index
    type(WeightPtrType),intent(in):: Weight

    integer :: iLat, iLon, iBlock, iPoint
    real    :: w

    Buff_V = 0.0

    do iPoint = iPointStart, iPointStart + nPoint - 1

       iLat   = Index % iCB_II(1,iPoint)
       iLon   = Index % iCB_II(2,iPoint)
       iBlock = Index % iCB_II(3,iPoint)
       w      = Weight % Weight_I(iPoint)

       if(iBlock/=1)then
          write(*,*)NameSub,': iPoint,Index % iCB_II=',&
               iPoint,Index%iCB_II(:,iPoint)
          call CON_stop(NameSub//&
               ' SWMF_ERROR iBlock should be 1=North in IM-IE coupling')
       end if

       if (iLon > nLon) iLon = mod(iLon,nLon)
       if(iLat<1 .or. iLat>nLat*2 .or. iLon<0 .or. iLon>nLon)then
          write(*,*)'iLat,iLon=',iLat, nLat, iLon, nLon
          call CON_stop(NameSub//' SWMF_ERROR index out of range')
       end if

       ! Only worry about the northern hemisphere....  IE can fix the southern hemisphere.

       if (iLat <= nLat .and. iLon <= nLon) then
          Buff_V(1) = Buff_V(1) - w * FAC_C(iLat,iLon)/2.0 ! / 1.0e6
          ! Fill with -1 for now. In the future we will determine these values
          Buff_V(2) = -1.0
          Buff_V(3) = -1.0

          !        Buff_V(1) = Buff_V(1) - w * birk(iLat,iLon)/2 / 1.0e6
          !        Buff_V(2) = Buff_V(2) + w * eflux(iLat,iLon,1)
          !        Buff_V(3) = Buff_V(3) + w * eavg(iLat,iLon,1) / 1000.0
       endif

       !     if (iLat > nLat .and. iLon <= nLon) &
       !          Buff_V(1) = Buff_V(1) + w * birk(2*nLat-iLat+1,iLon)/2 / 1.0e6

       !     if (iLat > IONO_nTheta .and. iLon <= IONO_nPsi) &
       !          Buff_V(1) = Buff_V(1) + w * IONO_SOUTH_RCM_JR(2*IONO_nTheta-iLat+1,iLon)

    end do

  end subroutine IM_get_for_ie

end module IM_wrapper
