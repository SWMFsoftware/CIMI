Module ModGmCimi

  use ModCimiGrid,  ONLY: nLat => np, nLon => nt, MinLonPar, MaxLonPar
  use ModCimiPlanet,ONLY: nspec

  implicit none

  real, allocatable :: StateLine_VI(:,:),StateBmin_IIV(:,:,:)
  integer :: iLineIndex_II(nLon,1:nLat),nPoint, nVarBmin
  integer :: TotalRho_=-1, TotalP_=-1, TotalPpar_=-1, TotalPe_=-1
  integer, dimension(nspec-1) :: iBufferRho_I 
  integer, dimension(nspec-1) :: iBufferP_I
  integer, dimension(nspec-1) :: iBufferPpar_I   
  real :: rBodyGm=2.5

  ! index of conjugate lat and lon in StateBmin array
  integer :: ConjugateLat_=4,ConjugateLon_=5
  logical :: UseTotalRhoGm = .false.
  logical :: UseTotalPGm   = .false.
  logical :: UseTotalPparGm= .false.
  logical :: UseMultiRhoGm = .false.
  logical :: UseMultiPGm   = .false.
  logical :: UseMultiPparGm= .false.
  logical :: DoFeedbackPs  = .false.
  logical :: UsePeGm       = .false.
  !number of variables extacted along field
  integer,parameter :: nVar=4

  real :: KpGm=-1.0, AeGm=-1.0
  logical :: UseGmKp=.false.,UseGmAe=.false.
  
  real, allocatable :: Den_IC(:,:,:), Temp_IC(:,:,:), &
       Temppar_IC(:,:,:)
  integer :: iLatMin=22 !Minimum latitude in MHD boundary
  
  logical :: UseGm                  = .false.
  logical :: DoneGmCoupling         = .false.

 public :: init_cimi_gm
contains
 subroutine init_cimi_gm
   integer :: iSpecies
   !----------------------------------------------------------------------
   if(allocated(Den_IC)) RETURN

   allocate(Den_IC(nLat,MinLonPar:MaxLonPar,nspec), &
            Temp_IC(nLat,MinLonPar:MaxLonPar,nspec), &
            Temppar_IC(nLat,MinLonPar:MaxLonPar,nspec))
 end subroutine init_cimi_gm
  
end Module ModGmCimi
