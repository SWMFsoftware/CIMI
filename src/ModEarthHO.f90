Module ModCimiPlanet
  implicit none

  real :: re_m, dipmom
  real :: Hiono = 120 ! ionosphere altitude in km
  
  !Species Information: Number of Ion Species (Earth=H+,O+)
  integer,parameter :: nspec=3  ! number of species

  ! a0,a1,... are coef. of a polynomial which defines the 
  ! exponent of the charge exchange cross section
  real,parameter,dimension(nspec-1) :: a0_I = (/-18.767,-18.987/)
  real,parameter,dimension(nspec-1) :: a1_I = (/-0.11017,-0.10613/)
  real,parameter,dimension(nspec-1) :: a2_I = (/-3.8173e-2,-5.4841e-3/)
  real,parameter,dimension(nspec-1) :: a3_I = (/-0.1232,-1.6262e-2/)
  real,parameter,dimension(nspec-1) :: a4_I = (/-5.0488e-2,-7.0554e-3/)

  !set species Mass in amu
  real,parameter,dimension(nspec)   :: amu_I = (/1.0,16.0,5.4462e-4/)


  !set density and temp factor
  real, dimension(nspec) :: dFactor_I =(/0.85,0.15,1.0/)
  real, dimension(nspec) :: tFactor_I =(/1.0,1.0,0.128205/)

  !set plot parameters
  character(len=300), parameter :: & 
       NamePlotVar='x y P[nP] HpP[nP] OpP[nP] eP[nP] ' &
       //'Ppar[nP] HpPpar[nP] OpPpar[nP] ePpar[nP] ' &
       //'Phot[nP] HpPhot[nP] OpPhot[nP] ePhot[nP] ' &
       //'Pparhot[nP] HpPparhot[nP] OpPparhot[nP] ePparhot[nP] ' &
       //'N[/m3] HpN[/m3] OpN[/m3] eN[/m3] ' &
       //'Beq[T] Vol[m3/Wb] Pot[Volts] FAC[Amp/m2] Lstar[Re] Plas[/m3] ' &
       //'g rbody'

  integer, dimension(nspec+1) :: iPplot_I    =(/1,2,3,4/) 
  integer, dimension(nspec+1) :: iPparplot_I =(/5,6,7,8/)
  integer, dimension(nspec+1) :: iPhotplot_I =(/9,10,11,12/)
  integer, dimension(nspec+1) :: iPparhotplot_I =(/13,14,15,16/)
  integer, dimension(nspec+1) :: iNplot_I    =(/17,18,19,20/)
  integer, parameter          :: Beq_=21,Vol_=22,Pot_=23, FAC_=24, Lstar_=25, &
                                 Plas_=26, nVar=26

  !set Logplot parameters

!  integer, parameter :: nLogVars = 8
  character(len=500), parameter :: & 
       NamePlotVarLog=&
       'it t dst '// &
       'RbSumH RcSumH HpDrift HpBfield HpChargeEx HpWaves '// &
       'HpStrongDiff HpDecay HpLossCone HpDriftIn HpDriftOut '// &
       'RbSumO RcSumO OpDrift OpBfield OpChargeEx OpWaves '// &
       'OpStrongDiff OpDecay OpLossCone OpDriftIn OpDriftOut '// &
       'RbSume RcSume  eDrift  eBfield  eChargeEx  eWaves '// &
       ' eStrongDiff  eDecay  eLossCone  eDriftIn  eDriftOut'


end Module ModCimiPlanet
