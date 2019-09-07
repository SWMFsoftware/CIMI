Module ModCimi
  use ModCimiGrid
  use ModCimiPlanet,ONLY: nspec
  implicit none

  SAVE

  real    :: dt=1., dtmax=1. ! typical time step of cimi
  real    :: Time = 0.0
  logical :: UseMcLimiter = .false., UseStrongDiff = .false., UseDecay =.false.
  real    :: BetaLimiter = 1.5
  real    :: Pmin = 1e-6
  real    :: DecayTimescale = 36000. ! Seconds
  real, allocatable:: SDtime(:,:,:,:), f2(:,:,:,:,:)
  real, allocatable:: phot(:,:,:), Ppar_IC(:,:,:), &
       Pressure_IC(:,:,:), PressurePar_IC(:,:,:)
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
       preF(:,:,:,:), preP(:,:,:,:), Eje1(:,:,:) ! presipitation output
  integer, parameter :: &
       nOperator = 8, OpDrift_ = 1, OpBfield_ = 2, OpChargeEx_ = 3, &
       OpWaves_ = 4, OpStrongDiff_ = 5, OpDecay_ = 6, &
       OpLossCone_ = 7, OpLossCone0_ = 8
! Note order and number of operators has been changed Waves are added
! and OpLossCone0_=8 is previous in time OpLossCone (needed for
! precipitation, see cimi_precip_calc subroutine in cimi.f90)

! Variables brought over from Module BoundaryCheck in
! ModCimiBoundary. neng can now have a variable size making it
! necessary that these variables be allocatable.
  
  real, allocatable :: &
       vdr_q1(:,:,:), vdr_q3(:,:,:), vgyr_q1(:,:,:), vgyr_q3(:,:,:)
  real, allocatable :: &
       eng_q1(:,:,:), eng_q3(:,:,:), vexb(:,:,:), dif_q1(:,:,:)
  real, allocatable :: &
       Part_phot(:,:,:,:), dif_q3(:,:,:)

! in CIMI: eChangeOperator=xle(ns,ir,ip,je+2) Here we create one array
!   for all operators (plus one dimension)
!   eTimeAccumulatv=esum(ns,ir,ip,je+2) No need for additional
!   dimention because it's total energy

  logical :: IsStandAlone = .false.
  logical :: DoCalcPrecip = .false.
  logical :: IsStrictDrift = .false.
  real :: DtCalcPrecip = 10.
  
contains

  subroutine init_mod_cimi

    if(allocated(f2)) RETURN

    allocate(                         &
         f2(nspec,np1,nt1,nm,nk),     &
         SDtime(np,nt,nm,nk),         &
         phot(nspec,np,nt),           &
         Ppar_IC(nspec,np,nt),        &
         Pressure_IC(nspec,np,nt),    &
         PressurePar_IC(nspec,np,nt), &
         FAC_C(np,nt),                &
         Bmin_C(np,nt))
    ! minimum B field along each field line
    ! passed from GM to IM, now as an output of IM

    ! now allocate arrays for energy tracking
    allocate(eChangeOperator_VICI(nspec,np,nt,neng+2,nOperator), &
         driftin(nspec), driftout(nspec), &
         rbsumLocal(nspec),rbsumGlobal(nspec), &
         rcsumLocal(nspec),rcsumGlobal(nspec))

    allocate(pChangeOperator_VICI(nspec,np,nt,neng+2,nOperator), &
         eChangeLocal(nspec,nOperator),eChangeGlobal(nspec,nOperator), &
         eTimeAccumult_ICI(nspec,np,nt,neng+2), &
         pTimeAccumult_ICI(nspec,np,nt,neng+2))
   
    allocate(preP(nspec,np,nt,neng+2), preF(nspec,np,nt,neng+2), &
         Eje1(nspec,np,nt))

    allocate(energy(nspec,neng), Ebound(nspec,neng), delE(nspec,neng))

    allocate( vdr_q1( nspec, np, nt ), vdr_q3( nspec, np, nt ), &
         vgyr_q1( nspec, np, nt ), vgyr_q3( nspec, np, nt ) )
    allocate( eng_q1( nspec, np, nt ), eng_q3( nspec, np, nt ), &
         vexb( nspec, np, nt ), dif_q1( nspec, np, nt ) )
    allocate( Part_phot( nspec, np, nt, neng ), dif_q3( nspec, np, nt ) )

  end subroutine init_mod_cimi

end Module ModCimi
