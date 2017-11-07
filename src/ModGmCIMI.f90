Module ModGmCimi

  use ModCimiGrid,  ONLY: nLat => np, nLon => nt
  use ModCimiPlanet,ONLY: nspec

  implicit none

  real, allocatable :: StateLine_VI(:,:),StateBmin_IIV(:,:,:)
  integer :: iLineIndex_II(nLon,1:nLat),nPoint, nVarBmin
  integer, parameter :: AveDens_=4, AveP_=5, AvePpar_=7, AveHpRho_=7, &
       AveOpRho_=8, AveHpP_=9, AveOpP_=10
  integer, parameter,dimension(nspec-1) :: AveDen_I = (/7,8/)
  integer, parameter,dimension(nspec-1) :: AveP_I   = (/9,10/)

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

end Module ModGmCimi
