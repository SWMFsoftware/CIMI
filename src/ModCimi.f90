Module ModCimi
  use ModCimiGrid
  use ModCimiPlanet,ONLY: nspec
  implicit none

  SAVE

  real    :: dt=1., dtmax=1. ! typical time step of cimi
  real    :: Time = 0.0
  logical :: UseMcLimiter = .false., UseStrongDiff = .false., UseDecay =.false.
  logical :: UseFLC=.false., useElectronDecay=.false.
  real    :: BetaLimiter = 1.5
  integer :: minLon, maxLon, nTotalCell

  real    :: Pmin = 1e-6

  real    :: DecayTimescale = 36000. ! Seconds
  real, allocatable:: SDtime(:,:,:,:), f2(:,:,:,:,:)
  real, allocatable:: phot(:,:,:), Ppar_IC(:,:,:), &
       Pressure_IC(:,:,:), PressurePar_IC(:,:,:), PressureOut_IC(:,:,:)
  real, allocatable:: FAC_C(:,:)
  real, allocatable:: Bmin_C(:,:)

  !Variable for E gain and loss
  real, allocatable :: driftin(:), driftout(:), &
       rbsumLocal(:), rbsumGlobal(:), rcsumLocal(:), rcsumGlobal(:)
  real, allocatable :: &
       eChangeOperator_VICI(:,:,:,:,:), pChangeOperator_VICI(:,:,:,:,:)
  real, allocatable :: eChangeLocal(:,:), eChangeGlobal(:,:)
  real, allocatable :: eTimeAccumult_ICI(:,:,:,:), pTimeAccumult_ICI(:,:,:,:)
  real, allocatable :: energy(:,:), Ebound(:,:), delE(:,:)
  real, allocatable :: &
       preF(:,:,:,:), preP(:,:,:,:), Eje1(:,:,:) ! precipitation output
  real, allocatable :: precipEnergyLoss(:,:,:,:), &
                       precipNumberLoss(:,:,:,:)
  integer, parameter :: &
       nOperator = 10, OpDrift_ = 1, OpBfield_ = 2, OpChargeEx_ = 3, &
       OpWaves_ = 4, OpStrongDiff_ = 5, OpDecay_ = 6, &
       OpLossCone_ = 7, OpFLC_ = 8, OpLossCone0_ = 9, OpFLC0_ = 10
! Note order and number of operators has been changed Waves are added
! and OpLossCone0_=9 is previous in time OpLossCone, same for OpFLC0_=10

! Variables brought over from Module BoundaryCheck in
! ModCimiBoundary. neng can now have a variable size making it
! necessary that these variables be allocatable.

  real, allocatable :: &
       vdr_q1(:,:,:), vdr_q3(:,:,:), vgyr_q1(:,:,:), vgyr_q3(:,:,:)
  real, allocatable :: &
       eng_q1(:,:,:), eng_q3(:,:,:), vexb(:,:,:), dif_q1(:,:,:)
  real, allocatable :: &
       Part_phot(:,:,:,:), dif_q3(:,:,:)

  ! Values brought over from cimi_run, need to be allocated to reduce
  ! stack size
  real, allocatable :: flux(:,:,:,:,:), psd(:,:,:,:,:), &
        vlEa(:,:,:,:,:), vpEa(:,:,:,:,:)
  real, allocatable :: achar(:,:,:,:,:)
  real, allocatable :: vl(:,:,:,:,:), vp(:,:,:,:,:), &
        fb(:,:,:,:)

  ! Global output array for proc 0 only
  real, allocatable :: FACGlobal_C(:,:), PressureGlobal_IC(:,:,:), PressureParGlobal_IC(:,:,:), &
    BminGlobal_C(:,:), photGlobal(:,:,:), PparGlobal_IC(:,:,:), DenGlobal_IC(:,:,:), &
    prePGlobal(:,:,:,:), preFGlobal(:,:,:,:), Eje1Global(:,:,:), &
    vdr_q3Global(:,:,:), eng_q3Global(:,:,:), part_photGlobal(:,:,:,:), vexbGlobal(:,:,:), dif_q3Global(:,:,:), &
    fluxGlobal(:,:,:,:,:), psdGlobal(:,:,:,:,:), vlEaGlobal(:,:,:,:,:), vpEaGlobal(:,:,:,:,:)

! in CIMI: eChangeOperator=xle(ns,ir,ip,je+2) Here we create one array
!   for all operators (plus one dimension)
!   eTimeAccumulatv=esum(ns,ir,ip,je+2) No need for additional
!   dimention because it's total energy

  logical :: IsStandAlone = .false.
  logical :: DoCalcPrecip = .true.
  logical :: IsStrictDrift = .false.
  real :: DtCalcPrecip = 10.
contains

  subroutine init_mod_cimi

    if(allocated(f2)) RETURN

    if (iProc < mod(nt,nProc))then
      nLonPar=(nt+nProc-1)/nProc
    else
      ! nLonPar is the number of points per proccess
      nLonPar=nt/nProc
    endif
    minLon = -nGhostCell+1
    maxLon = nLonPar+nGhostCell
    nTotalCell = maxLon - minLon

    allocate(                         &
         phi(minLon:maxLon), xmlt(minLon:maxLon), &
         f2(np,minLon:maxLon,nm,nk,nspec),     &
         SDtime(np,minLon:maxLon,nm,nk),         &
         phot(np,minLon:maxLon,nspec),           &
         Ppar_IC(np,minLon:maxLon,nspec),        &
         Pressure_IC(np,minLon:maxLon,nspec),    &
         PressurePar_IC(np,minLon:maxLon,nspec), &
         FAC_C(np,minLon:maxLon),                &
         Bmin_C(np,minLon:maxLon),              &
         d4Element_C(np,nm,nk,nspec))
    f2 = -912392.0
    ! minimum B field along each field line
    ! passed from GM to IM, now as an output of IM

    ! now allocate arrays for energy tracking
    allocate(eChangeOperator_VICI(np,minLon:maxLon,neng+2,nOperator,nspec), &
         driftin(nspec), driftout(nspec), &
         rbsumLocal(nspec),rbsumGlobal(nspec), &
         rcsumLocal(nspec),rcsumGlobal(nspec))

    allocate(pChangeOperator_VICI(np,minLon:maxLon,neng+2,nOperator,nspec), &
         eChangeLocal(nOperator,nspec),eChangeGlobal(nOperator,nspec), &
         eTimeAccumult_ICI(np,minLon:maxLon,neng+2,nspec), &
         pTimeAccumult_ICI(np,minLon:maxLon,neng+2,nspec))

    allocate(preP(np,minLon:maxLon,neng+2,nspec), preF(np,minLon:maxLon,neng+2,nspec), &
         Eje1(np,minLon:maxLon,nspec), precipEnergyLoss(np,minLon:maxLon,neng+2,nspec), &
         precipNumberLoss(np,minLon:maxLon,neng+2,nspec))
    preP = 0.0
    preF = 0.0
    Eje1 = 0.0
    allocate(flux(np,minLon:maxLon,neng,npit,nspec), psd(np,minLon:maxLon,nm,nk,nspec), &
        vlEa(np,minLon:maxLon,neng,npit,nspec), vpEa(np,minLon:maxLon,neng,npit,nspec))
    allocate(vl(0:np,minLon:maxLon,nm,nk,nspec), vp(0:np,minLon:maxLon,nm,nk,nspec), &
        fb(minLon:maxLon,nm,nk,nspec), achar(0:np,minLon:maxLon,nm,nk,nspec))
    vl = 0.0
    vp = 0.0

    allocate(energy(neng,nspec), Ebound(neng,nspec), delE(neng,nspec))

    allocate( vdr_q1( np, minLon:maxLon, nspec ), vdr_q3( np, minLon:maxLon, nspec ), &
         vgyr_q1( np, minLon:maxLon, nspec ), vgyr_q3( np, minLon:maxLon, nspec ) )
    allocate( eng_q1( np, minLon:maxLon, nspec ), eng_q3( np, minLon:maxLon, nspec ), &
         vexb( np, minLon:maxLon, nspec ), dif_q1( np, minLon:maxLon, nspec ) )
    allocate( Part_phot( np, minLon:maxLon, neng, nspec ), dif_q3( np, minLon:maxLon, nspec ) )

    if(iProc == 0) then
      allocate( FACGlobal_C(np,nt), PressureGlobal_IC(np,nt,nspec), PressureParGlobal_IC(np,nt,nspec), &
           BminGlobal_C(np,nt), photGlobal(np,nt,nspec), PparGlobal_IC(np,nt,nspec), DenGlobal_IC(np,nt,nspec), &
           prePGlobal(np,nt,neng+2,nspec), preFGlobal(np,nt,neng+2,nspec), Eje1Global(np,nt,nspec), &
           vdr_q3Global(np,nt,nspec), eng_q3Global(np,nt,nspec), part_photGlobal(np,nt,neng+2,nspec), &
           vexbGlobal(np,nt,nspec), dif_q3Global(np,nt,nspec), &
           fluxGlobal(np,nt,neng,npit,nspec), psdGlobal(np,nt,nm,nk,nspec), &
           vlEaGlobal(np,nt,neng,npit,nspec), vpEaGlobal(np,nt,neng,npit,nspec))
    endif

  end subroutine init_mod_cimi

end Module ModCimi
