Module ModGmCimi

  use ModCimiGrid,  ONLY: nLat => np, nLon => nt
  use ModCimiPlanet,ONLY: nspec

  implicit none

  real, allocatable :: StateLine_VI(:,:),StateBmin_IIV(:,:,:)
  integer :: iLineIndex_II(nLon,1:nLat),nPoint, nVarBmin
  integer, parameter :: AveDens_=4, AveP_=5, AvePpar_=7
  integer, dimension(nspec-1) :: AveDen_I 
  integer, dimension(nspec-1) :: AveP_I   

  integer,parameter :: nVar=4

  real :: KpGm=-1.0
  logical :: UseGmKp=.false.
  
  real :: Den_IC(nspec,nLat,nLon) = 0.0, Temp_IC(nspec,nLat,nLon) = 0.0, &
       Temppar_IC(nspec,nLat,nLon) = 0.0
  integer :: iLatMin=22 !Minimum latitude in MHD boundary
  
  logical :: UseGm                  = .false.
  logical :: DoneGmCoupling         = .false.
  logical :: DoMultiFluidGMCoupling = .false.
  logical :: DoAnisoPressureGMCoupling      = .false.

  public :: init_gm_cimi
contains
  subroutine init_gm_cimi
    integer :: iSpecies
    !----------------------------------------------------------------------
    do iSpecies =1,nspec-1
       AveDen_I(iSpecies) = 6+iSpecies
       AveP_I(iSpecies)=6+nspec-1+iSpecies
    enddo
  end subroutine init_gm_cimi
  
end Module ModGmCimi
