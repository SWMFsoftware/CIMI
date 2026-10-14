Module ModCimiGrid
 
  use ModCimiPlanet,   ONLY: nspec

  implicit none
  
  logical, parameter :: DoUseUniformLGrid = .false.
  logical :: DoVerboseLatGrid = .true.
  logical :: DoSingleGrid = .true.
  
  ! define dimensions of CIMI grids
  integer,parameter :: np1=1,nt1=1 
  integer,parameter :: nm=48,nk=29 ! dimension of CIMI magnetic moment and K
  integer,parameter :: nGhostCell=0

  integer,parameter :: np=1    ! dimension of the CIMI latitude grid
  integer,parameter :: nt=1    ! dimension of the CIMI local-time grid
  integer	    :: neng=15  ! dimension of the CIMI energy grid
  integer,parameter :: npit=18  ! dimension of the CIMI pitch-angle grid

  ! These have to be initialized so that IM_set_grid does not fail on non-IM PEs
  real:: xlat(np) = 0.0, phiGlobal(nt)=0.0


  real :: xlatr(np), dlat(np), sinAo(npit), xmltGlobal(nt), dphi
  real, allocatable :: phi(:), xmlt(:)
  logical :: DoReadParamGrid = .false.
  real :: LatitudeDeg=62.0, dLatitudeDeg=1.0,&
          LongitudeRad=0.0, dLongitudeRad=1.0
  integer :: minLon, maxLon, nTotalCell

  real :: varL(0:np+1)=0.,&
          dvarL=1.,&
          Lfactor(0:np+1) = 0.,&
          Lfactor1(0:np+1) = 0.,&
          xlatmin=9.8403398252184129,&
          xlatmax=72.4356255492731975,&
          varNpower=1.5
  logical :: DoDefineVarNpower = .false.

  real, allocatable :: d4Element_C(:,:,:,:) !4D element (dlat*dphi*dmm*dk)

  ! Defines the outer boundary maximum distance.
  real	  :: rb = 10.0
  
  ! Define MPI parameters affecting grid
  integer :: iProc, nProc, iComm, nLonPar, MinLonPar, MaxLonPar, &
     MinLonGlobal, MaxLonGlobal
  integer, allocatable :: nLonPar_P(:),nLonBefore_P(:)
  integer ::iProcLeft, iProcRight, iLonLeft, iLonRight
  integer :: iProcMidnight, iLonMidnight

  ! Define parameters for specifying control over Energy grid
  real :: MinIonEnergy=0.1, MaxIonEnergy=316.22777, ElectronGridRatio=10.
  logical :: UseLogEGrid = .true.

  ! Define parameters for specifying Energy grid corresponding to RBSP
  ! flux (MagEIS and REPT) instruments
  integer, parameter :: nRBSPEnergy = 30
  real, dimension(1:nRBSPEnergy) :: &
       energy_RBSP = (/  19.5,   32.7,    50.0,    71.8,    98.1, &
                        129.0,  146.0,   166.0,   207.0,   232.0, &
                        334.0,  450.0,   580.0,   722.0,   877.0, &
                       1040.0, 1200.0,  1740.0,  2300.0,  2520.0, &
                       2850.0, 3600.0,  3770.0,  4500.0,  5600.0, &
                       7150.0, 8800.0, 11650.0, 22550.0, 59450.0 /)
  logical :: UseRBSPGrid = .false.

contains

!-----------------------------------------------------------------------------
 subroutine init_lat
 ! This subroutine initializes and defines latitude grid
!-----------------------------------------------------------------------------
 use ModNumConst,    ONLY: cDegToRad
 implicit none

 integer i

 ! CIMI xlat grid for non-uniform grid
 xlat = LatitudeDeg
 dlat = 1.0
 xlatr(:)=xlat(:)*cDegToRad

 end subroutine init_lat
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
 subroutine init_mlt
 ! This subroutine initializes and defines MLT grid
!-----------------------------------------------------------------------------
 use ModNumConst,    ONLY: cDegToRad,cPi
 implicit none

 integer i

 ! CIMI xmlt grid
 dphi = dLongitudeRad
 phiGlobal=LongitudeRad
 phi = LongitudeRad

 xmltGlobal=mod(phiGlobal*12.0/cPi + 12.0,24.0)
 xmlt=mod(phi*12.0/cPi + 12.0,24.0)

 iLonMidnight = 1

 end subroutine init_mlt
!-----------------------------------------------------------------------------
     
end Module ModCimiGrid


