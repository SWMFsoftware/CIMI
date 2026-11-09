Module ModCimiInitialize
  use ModCimiGrid,	ONLY: nm, nk, neng,npit,np
  use ModCimiPlanet,	ONLY: nspec
  implicit none
  
  real		:: &
       xmm(0:nm+1,nspec), xk(0:nk+1), dmm(nm,nspec), dk(nk),&
       dmu(npit), xjac(np,nm,nspec)

  logical	:: &
       IsEmptyInitial = .false., IsDataInitial = .false., &
       IsRBSPData = .false., IsGmInitial = .true., &
       IsInitialElectrons = .false., DoLstarInitialization = .false.

  character(16) :: NameFinFile = 'quiet_x.fin'
  character(11)	:: FinFilePrefix = 'xxxxxxxx'

end Module ModCimiInitialize
