Module ModCimiInitialize
  use ModCimiGrid ,ONLY: nm, nk, neng,npit1,np1
  use ModCimiPlanet,ONLY: nspec
  implicit none
  
  real :: xmm(nspec,0:nm+1),xk(0:nk+1),dphi,dmm(nspec,nm),dk(nk),&
          dmu(npit1),xjac(nspec,np1,nm),varL(0:np1+1),dvarL,varNpower,&
          Lfactor(0:np1+1),Lfactor1(0:np1+1)

  logical :: IsEmptyInitial=.false., IsDataInitial=.false., &
       IsRBSPData=.false., IsGmInitial=.true.,DoDefineVarNpower=.false.

end Module ModCimiInitialize
