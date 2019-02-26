Module ModCimiInitialize
  use ModCimiGrid ,ONLY: nm, nk, neng,npit1,np1
  use ModCimiPlanet,ONLY: nspec
  implicit none
  
  real :: xmm(nspec,0:nm+1),xk(0:nk+1),dphi,dmm(nspec,nm),dk(nk),&
          dmu(npit1),xjac(nspec,np1,nm)

  ! Need to initialize these variables to avoid compilation errors
  real :: dvarL = 1., varNpower = 3., varL(0:np1+1) = 0., &
       Lfactor(0:np1+1) = 0., Lfactor1(0:np1+1) = 0.
  
  logical :: IsEmptyInitial=.false., IsDataInitial=.false., &
       IsRBSPData=.false., IsGmInitial=.true., DoDefineVarNpower = .false.

end Module ModCimiInitialize
