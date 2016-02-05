Module ModCrcmInitialize
  use ModCrcmGrid
  use ModCrcmPlanet,ONLY: nspec
  implicit none
  
  real :: xmm(nspec,0:nm+1),xk(0:nk+1),dphi,dmm(nspec,nm),dk(nk),delE(nspec,neng),&
          dmu(npit1),xjac(nspec,np1,nm)

  logical :: IsEmptyInitial=.false., IsDataInitial=.false., &
       IsRBSPData=.false., IsGmInitial=.true.

end Module ModCrcmInitialize
