module ModIeCrcm
  use ModCrcmGrid,ONLY: np, nt
  implicit none 
  
  real :: pot(np,nt)
  logical :: UseIe =.true.
  logical :: UseWeimer =.false.
contains

  subroutine init_mod_iecrcm
    
    pot = 0.0

  end subroutine init_mod_iecrcm

end module ModIeCrcm
