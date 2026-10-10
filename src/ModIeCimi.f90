module ModIeCimi
  use ModCimiGrid,ONLY: np, nt, MinLonPar, MaxLonPar
  implicit none 
  
  real,allocatable :: pot(:,:)
  real :: potGlobal(np,nt)=0.0
  logical :: UseIe =.true.
  logical :: UseWeimer =.false.
contains

  subroutine init_mod_iecimi

    if(allocated(pot)) RETURN
    
    allocate(pot(np,MinLonPar:MaxLonPar))
    pot = 0.0

  end subroutine init_mod_iecimi

end module ModIeCimi
