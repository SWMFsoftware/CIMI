Module ModCrcmInitialize
  use ModCrcmGrid
  use ModCrcmPlanet,ONLY: nspec
  implicit none
  
  real :: xmm(nspec,nm),xk(nk),dphi,dmm(nspec,nm),dk(nk),delE(nspec,neng),&
          dmu(npit1),xjac(nspec,np1,nm)

  logical :: IsEmptyInitial=.false., IsDataInitial=.false., &
       IsRBSPData=.false., IsGmInitial=.true.

end Module ModCrcmInitialize
