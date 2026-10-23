module ModCimiBoundary

  use ModCimiGrid,    ONLY: MinLonPar,MaxLonPar,nt,np
  use ModCimiPlanet,  ONLY: nSpec, amu_I, dFactor_I, tFactor_I
  use ModUtilities,   ONLY: CON_stop

  implicit none
  
  SAVE

  private ! EXCEPT:

  public :: cimi_set_boundary_mhd
  public :: cimi_set_boundary_empirical
  public :: init_cimi_boundary
  public :: cimi_set_boundary_f2
  
  !when in standalone mode which boundary type do you use for the RC?
  logical,public :: UseBoundaryEbihara =.true.
  logical,public :: UseYoungEtAl =.false.
  logical,public :: CIMIboundary = .false.
  logical,public :: Outputboundary = .false.

  ! Dynamic boundary kappa based on Eyelade (2021)
  logical,public :: UseDynamicKappa  = .false.
  logical,public :: UseKappaTotalFit = .true.   ! temp to use for kappa fit
  logical,public :: DoRescaleKappaTemp = .true.  ! rescale T to conserve energy
  real, public :: KappaConst_I(nSpec) = -100, KappaMin_I(nSpec) = 2.0, &
    KappaMax_I(nSpec) = 30.0
  
  !boundary density in m^-3 and temperature in eV
  real, public, allocatable :: BoundaryDens_IC(:,:),BoundaryTemp_IC(:,:),&
   BoundaryTempPar_IC(:,:)

contains
  !============================================================================
  subroutine init_cimi_boundary
    use ModCimiGrid,    ONLY: MinLonPar,MaxLonPar
    use ModCimiPlanet,  ONLY: nSpec

    if(allocated(BoundaryDens_IC)) RETURN
    allocate(BoundaryDens_IC(MinLonPar:MaxLonPar,nSpec), &
             BoundaryTemp_IC(MinLonPar:MaxLonPar,nSpec), &
             BoundaryTempPar_IC(MinLonPar:MaxLonPar,nSpec))

    ! default kappa for electrons
    if(KappaConst_I(nSpec) == -100.) KappaConst_I(nSpec) = 3.0  

  end subroutine init_cimi_boundary
  !============================================================================
  subroutine cimi_set_boundary_mhd
    use ModGmCimi,      ONLY: Den_IC,Temp_IC,Temppar_IC,StateBmin_IIV,&
         TotalP_,TotalPpar_,TotalRho_, iBufferRho_I,iBufferP_I,&
         iBufferPpar_I, iLatMin,&
         UseTotalRhoGm, UseTotalPGm, UseTotalPparGm, UseMultiRhoGm , &
         UseMultiPGm, UseMultiPparGm, DoFeedbackPs, UsePeGm, TotalPe_
    use ModCimiGrid, ONLY: MinLonGlobal 
    
    use ModCimiTrace,  ONLY: irm, iba
    use ModConst, ONLY: cProtonMass
    integer :: iSpecies, iLon, iLat, jLat, iLonGlobal
    real, parameter :: cJouleToEv=6.2415e18  
    real :: FactorTotalDens
    integer :: ib1
    !---------------------------------------------------------------------------

    ! If using separate GM electron pressure, all values should just be passed
    ! without any temperature partitioning
    if (UsePeGm) tFactor_I = 1.0
    
    !set boundary density and temperature inside irm
    if (UseTotalRhoGm &
         .and. .not. UseMultiRhoGm &
         .and. .not. UseMultiPGm   &
         .and. .not. UseMultiPparGm) then
     ! When not Multifluid we get the total density from Rho_MHD as follows:
     ! Rho = (m1n1+m2n2+...) = n * sum(m_i*dFactor_i)
     ! where sum(dFactor_i)=1 (over ions) and n_i=dFactor_i*n 
     ! n_i = dFactor_i*Rho_MHD/(sum(m_i*dFactor_i))
     ! n_total = Rho_MHD/sum(m_i*dfactor_i)
     ! FactorTotalDens = sum(m_i*dfactor_i)
     FactorTotalDens = sum(dFactor_I(1:nSpec-1)*amu_I(1:nSpec-1))
     do iSpecies = 1, nSpec
        do iLon=MinLonPar,MaxLonPar
           iLonGlobal = MinLonGlobal + iLon - 1
           ! Set boundary index consistant with boundaryIM
           ib1=min(iba(iLon)+1,irm(iLon))
           do iLat=1,irm(iLon)
              ! Use boundary values for iLat>iLatMin
              jLat = max(iLat, iLatMin)

              !write(*,*)StateBmin_IIV(jLat,iLon,TotalRho_)
              Den_IC(iLat,iLon,iSpecies) = dFactor_I(iSpecies) * &
                   StateBmin_IIV(jLat,iLonGlobal,TotalRho_)&
                   /(cProtonMass * FactorTotalDens)
              Temp_IC(iLat,iLon,iSpecies) = tFactor_I(iSpecies) * &
                   StateBmin_IIV(jLat,iLonGlobal,TotalP_) &
                   * cProtonMass * FactorTotalDens &
                   / StateBmin_IIV(jLat,iLonGlobal,TotalRho_) &
                   * cJouleToEv
              if(UseTotalPparGm) then
                 Temppar_IC(iLat,iLon,iSpecies) = tFactor_I(iSpecies) * &
                      StateBmin_IIV(jLat,iLonGlobal,TotalPpar_) &
                      * cProtonMass * FactorTotalDens &
                      / StateBmin_IIV(jLat,iLonGlobal,TotalRho_) &
                      * cJouleToEv
              else
                 Temppar_IC(iLat,iLon,iSpecies) = Temp_IC(iLat,iLon,iSpecies)
              endif
           end do
        end do
     end do
  elseif (UseMultiRhoGm  &
       .and. UseTotalPGm &
       .and. .not. UseMultiPGm   &
       .and. .not. UseMultiPparGm) then
     ! Multispecies Case
     
     do iLon=MinLonPar,MaxLonPar
        iLonGlobal = MinLonGlobal + iLon - 1
        ! Set boundary index consistant with boundaryIM
        ib1=min(iba(iLon)+1,irm(iLon))
        do iLat=1,irm(iLon) 
           ! Use boundary values for iLat>iLatMin
           jLat = max(iLat, iLatMin)
           do iSpecies = 1, nSpec-1
              Den_IC(iLat,iLon,iSpecies) = &
                   StateBmin_IIV(jLat,iLonGlobal,iBufferRho_I(iSpecies))&
                   / (cProtonMass*amu_I(iSpecies))
           enddo
           ! set electron density
           Den_IC(iLat,iLon,nSpec) = sum(Den_IC(iLat,iLon,1:nSpec-1))
           
           ! set temperatures
           Temp_IC(iLat,iLon,1:nSpec) = tFactor_I(1:nSpec) * &
                StateBmin_IIV(jLat,iLonGlobal,TotalP_)  &
                /(Den_IC(iLat,iLon,nSpec)) &
                * cJouleToEv
           if(UseTotalPparGm) then
              Temppar_IC(iLat,iLon,1:nSpec) = tFactor_I(1:nSpec) * &
                   StateBmin_IIV(jLat,iLonGlobal,TotalPpar_) &
                   / (Den_IC(iLat,iLon,nSpec)) &
                   * cJouleToEv
           else
              Temppar_IC(iLat,iLon,1:nSpec) = Temp_IC(iLat,iLon,1:nSpec)
           endif
           
        end do
     end do


  elseif (UseMultiRhoGm .and. UseMultiPGm ) then
     !Multifluid Case

     do iLon=MinLonPar,MaxLonPar
        iLonGlobal = MinLonGlobal + iLon - 1
        ! Set boundary index consistant with boundaryIM
        ib1=min(iba(iLon)+1,irm(iLon))
        do iLat=1,irm(iLon) 
           ! Use boundary values for iLat>iLatMin
           jLat = max(iLat, iLatMin)

           do iSpecies = 1, nSpec-1
              Den_IC(iLat,iLon,iSpecies) = &
                   StateBmin_IIV(jLat,iLonGlobal,iBufferRho_I(iSpecies))&
                   / (cProtonMass*amu_I(iSpecies))
              
              !set temperatures
              Temp_IC(iLat,iLon,iSpecies) = &
                   StateBmin_IIV(jLat,iLonGlobal,iBufferP_I(iSpecies))  &
                   /(Den_IC(iLat,iLon,iSpecies)) &
                   * cJouleToEv
              if(UseMultiPparGm) then
                 Temppar_IC(iLat,iLon,iSpecies) = &
                      StateBmin_IIV(jLat,iLonGlobal,iBufferPpar_I(iSpecies)) &
                      / (Den_IC(iLat,iLon,iSpecies)) &
                      * cJouleToEv
              else
                 Temppar_IC(iLat,iLon,iSpecies) = Temp_IC(iLat,iLon,iSpecies)
              endif
           enddo

           ! set electron density
           Den_IC(iLat,iLon,nSpec) = sum(Den_IC(iLat,iLon,1:nSpec-1))
           
           ! set electron temperatures
           Temp_IC(iLat,iLon,nSpec) = &
                sum(StateBmin_IIV(jLat,iLonGlobal,iBufferP_I(1:nSpec-1)))  &
                /(Den_IC(iLat,iLon,nSpec)) &
                * cJouleToEv
           Temppar_IC(iLat,iLon,nSpec) = Temp_IC(iLat,iLon,nSpec)
        end do
     end do
  else
     call CON_stop('Unknown condition in ModCimiBoundary')
  endif
  
  !Overwrite e temperature if solving for Pe in GM
  if (UsePeGm) then
     do iLon=MinLonPar,MaxLonPar
     iLonGlobal = MinLonGlobal + iLon - 1
        ! Set boundary index consistant with boundaryIM
        ib1=min(iba(iLon)+1,irm(iLon))
        do iLat=1,irm(iLon)
           jLat = max(iLat, iLatMin)
           Temp_IC(iLat,iLon,nSpec) = StateBmin_IIV(jLat,iLonGlobal,TotalPe_)&
                / Den_IC(iLat,iLon,nSpec)* cJouleToEv
           Temppar_IC(iLat,iLon,nSpec) = Temp_IC(iLat,iLon,nSpec)
        end do
     end do
  endif
  
  ! Fill the boundary n, T, and T||
  do iSpecies = 1,nSpec
     BOUNDARY_FILL: do iLon=MinLonPar,MaxLonPar
        ! Fill boundary arrays
        ib1=min(iba(iLon)+1,irm(iLon))
        !write(*,*) iSpecies, ib1,iLon,Temp_IC (iSpecies,ib1,iLon) 
        BoundaryDens_IC(iLon,iSpecies) = &
             Den_IC (ib1,iLon,iSpecies)
        BoundaryTemp_IC(iLon,iSpecies) = &
             Temp_IC(ib1,iLon,iSpecies)
        BoundaryTempPar_IC(iLon,iSpecies) = &
             Temppar_IC(ib1,iLon,iSpecies)
        
     end do BOUNDARY_FILL
  enddo

  end subroutine cimi_set_boundary_mhd


  !============================================================================
  subroutine cimi_set_boundary_empirical
    use ModImTime,    ONLY: CurrentTime, TimeLagBoundary
    use ModCimiTrace,ONLY: iba, ro, xmlto, irm
    use ModGmCimi,      ONLY: Den_IC,Temp_IC
    use ModImIndices,   ONLY: UseKpApF107IndicesFile,get_im_indices_F107,&
         get_im_indices_Kp
    use ModIndicesInterfaces
    use ModNumConst,       ONLY: cPi
    
    
    real,parameter :: cCm3toM3 =1.0e6, cKevTOeV=1000.0
    real :: TimeBoundary, Bz, vsw, xnsw, Rtsy, phit
    integer :: iError, iLon, ib1, iSpecies
    
    real :: xn1(MinLonPar:MaxLonPar), xkt1(MinLonPar:MaxLonPar)
    !variable for young et al density split
    real :: O_H_ratio, F107, Kp
    !--------------------------------------------------------------------------
    
    if(tFactor_I(1) < 1) tFactor_I = tFactor_I / tFactor_I(1)
    ! If using Young et al and nSpec=3 then overwrite dFactor_I
    if(UseYoungEtAl .and. nSpec==3) then
       !Get Inputs
       if (UseKpApF107IndicesFile) then
          call get_im_indices_Kp(CurrentTime, Kp)
          call get_im_indices_F107(CurrentTime, F107)
       else
          call get_kp(CurrentTime, Kp, iError)
          call get_F107(CurrentTime, F107, iError)
       endif
       !Set O+:H+ ratio
       O_H_ratio = get_OtoH_young(F107, Kp)
       dFactor_I(1) = 1.0/(O_H_ratio+1.0)        !H+
       dFactor_I(2) = O_H_ratio/(O_H_ratio+1.0)  !O+

    endif
    
    If (UseBoundaryEbihara) then
       ! Use Ebihara-Ejiri-Borovsky or Tsyganenko-Mukai model for 
       ! boundary n and kT
       ! Set TimeBoundary to CurrentTime with a time lag
       TimeBoundary=CurrentTime - TimeLagBoundary
       call get_IMF_Bz(TimeBoundary, Bz, iError)
       call get_SW_V  (TimeBoundary, vsw, iError)
       call get_SW_N  (TimeBoundary, xnsw, iError)
       
       ! Density in m^-3, Ebihara & Ejiri 2000
       xn1(MinLonPar:MaxLonPar)=(0.025*xnsw+0.395)*cCm3toM3
       
       ! kT_ion in eV, Borovsky etal 1998
       xkt1(MinLonPar:MaxLonPar)=(0.019*vsw-3.65)*cKevTOeV
    else
       ! use Tsyganenko-Mukai PS model
       call get_IMF_Bz(CurrentTime, Bz, iError)
       call get_SW_V  (CurrentTime, vsw, iError)
       call get_SW_N  (CurrentTime, xnsw, iError)
       do iLon=MinLonPar,MaxLonPar
          ! Set boundary index consistant with boundaryIM
          ib1=min(iba(iLon)+1,irm(iLon))
          Rtsy=ro(ib1,iLon)
          if (Rtsy.lt.10.) Rtsy=10.   ! model only for r > 10
          phit=xmlto(ib1,iLon)*cPi/12.
          call get_tsy_plasma(Bz,vsw,xnsw,Rtsy,phit,xkt1(iLon),xn1(iLon))
       enddo
    endif

    ! Fill boundary arrays for all species
    do iSpecies = 1, nSpec
       BoundaryDens_IC(MinLonPar:MaxLonPar,iSpecies) = dFactor_I(iSpecies)*xn1(MinLonPar:MaxLonPar)
       BoundaryTemp_IC(MinLonPar:MaxLonPar,iSpecies) = tFactor_I(iSpecies)*xkt1(MinLonPar:MaxLonPar)
       BoundaryTempPar_IC(MinLonPar:MaxLonPar,iSpecies) = BoundaryTemp_IC(MinLonPar:MaxLonPar,iSpecies)
    end do
 
 end subroutine cimi_set_boundary_empirical

 !============================================================================
 real function get_OtoH_young(F107, Kp) result(O_H_ratio)
   real, intent(in) :: F107, Kp
   O_H_ratio=4.5e-2*exp(0.17*Kp+0.01*F107)  ! O+/H+ ratio (Young etal 1982)
 end function get_OtoH_young
 
 !============================================================================
 !-----------------------------------------------------------------------------
 subroutine get_tsy_plasma(Bz,Vsw,Nsw,r,phi,t_tsy,n_tsy)
   !---------------------------------------------------------------------------
   ! Subroutine of Tsyganenko-Mukai plasma sheet model.
   ! Reference: Tsyganenko, N. A., and T. Mukai, Tail plasma sheet models 
   !               derived from Geotail particle data, J. Geophys. Res., 
   !               108(A3), 1136,doi:10.1029/2002JA009707, 2003.
   use ModIoUnit, ONLY: UnitTmp_
   implicit none
   integer, parameter :: n_points_tsy=16
   
   ! variable for reading and storing data parameters
   real, save :: a_tsy(6,16)
   real    :: coef1(6), tmp
   integer :: i

   real Bz,Bn,Bs,V,Vsw,r,rho,phi,N,Nsw,t_tsy,n_tsy,temp1,temp2,temp3
   logical :: IsFirstCall = .true.
   !---------------------------------------------------------------------------
   
   ! Read data for Tsyganenko-Mukai plasmasheet model
   if (IsFirstCall) then
      open(UnitTmp_,file='TSYGANENKO_MUKAI.dat',status='old')
      do i=1,n_points_tsy
         read(UnitTmp_,*) tmp,coef1
         a_tsy(1:6,i)=coef1(1:6)
      enddo
      close(UnitTmp_)
      IsFirstCall = .false.
   endif

   if (Bz.ge.0.) then
      Bn=Bz/5.
      Bs=0.
   else
      Bn=0.
      Bs=-Bz/5.
   endif
   
   V=Vsw/500.
   rho=r/10.
   N=Nsw/10.
   
   ! temperature final fit (equation (4), page 5)
   temp1=a_tsy(1,9)*(V**a_tsy(1,15))+a_tsy(1,10)*Bn+a_tsy(1,11)*Bs
   temp1=temp1*(-1.)*(rho-1.)
   temp2=a_tsy(1,12)*(V**a_tsy(1,16))+a_tsy(1,13)*Bn+a_tsy(1,14)*Bs
   temp2=temp2*(-1.)*(rho-1.)
   temp3=a_tsy(1,5)*V+a_tsy(1,6)*Bn+a_tsy(1,7)*Bs+a_tsy(1,8)*exp(temp2)
   t_tsy=a_tsy(1,1)*V+a_tsy(1,2)*Bn+a_tsy(1,3)*Bs+a_tsy(1,4)*exp(temp1)+ &
        (temp3*sin(phi)*sin(phi))
   
   ! density final fit (equation (5), page 6)
   temp1=a_tsy(3,1)+a_tsy(3,2)*(N**a_tsy(3,10))+a_tsy(3,3)*Bn+a_tsy(3,4)*V*Bs
   temp2=a_tsy(3,5)*(N**a_tsy(3,11))+a_tsy(3,6)*Bn+a_tsy(3,7)*V*Bs
   n_tsy=temp1*(rho**a_tsy(3,8)) + temp2*(rho*a_tsy(3,9))*sin(phi)*sin(phi)
   
 end subroutine get_tsy_plasma
 !============================================================================
 real function dynamic_kappa(Beta, EnergyKeV, iSpecies, UseTotalFit)

  ! Empirical kappa = A*E**gamma with A,gamma = a*|log10(beta/beta0)|+b
  ! from Eyelade et al. (2021, ApJS 253:34, Table 5), THEMIS observations.
  ! Energy is in KeV

  real,    intent(in):: Beta, EnergyKeV
  integer, intent(in):: iSpecies
  logical, intent(in):: UseTotalFit

  ! (a, beta0, b) for A and gamma
  ! note that while in paper, Core energy fit performs better against data,
  ! our core energy needs to be set based on the kappa, so using it is circular
  ! in nature. Default is to use total energy fit. 
  real, parameter:: CoefIonACore_I(3) = [1.37, 0.06, 0.76]
  real, parameter:: CoefIonATotal_I(3) = [1.26, 0.07, 0.84]
  real, parameter:: CoefIonGCore_I(3) = [-0.53, 0.09, 1.28]
  real, parameter:: CoefIonGTotal_I(3) = [-0.40, 0.11, 1.10]
  real, parameter:: CoefEleACore_I(3) = [2.21, 0.15, 2.25]
  real, parameter:: CoefEleATotal_I(3) = [1.69, 0.15, 1.83]  
  real, parameter:: CoefEleGCore_I(3) = [-0.48, 0.09, 1.17]
  real, parameter:: CoefEleGTotal_I(3) = [-0.68, 0.10, 1.33]

  real:: CoefA_I(3), CoefG_I(3), BetaFit, EnergyFit, A, Gamma
  !----------------------------------------------------------------------------
  ! Pick right coefficients based on species and fit type
  if (iSpecies == nSpec) then
     EnergyFit = min(max(EnergyKeV, 0.25), 9.25)
     if (UseTotalFit) then
        CoefA_I = CoefEleATotal_I
        CoefG_I = CoefEleGTotal_I
     else
        CoefA_I = CoefEleACore_I
        CoefG_I = CoefEleGCore_I
     endif
  else
     EnergyFit = min(max(EnergyKeV, 0.2), 15.25)
     if (UseTotalFit) then
        CoefA_I = CoefIonATotal_I
        CoefG_I = CoefIonGTotal_I
     else
        CoefA_I = CoefIonACore_I
        CoefG_I = CoefIonGCore_I
     end if   
  end if
  
  ! Restrict beta to observed range
  ! Possibly should be restricted further on the > 1 side, since data doesn't
  ! seem great for it
  BetaFit = min(max(Beta, 1e-3), 1e2)

  ! Set fit parameters
  A     = CoefA_I(1)*abs(log10(BetaFit/CoefA_I(2))) + CoefA_I(3)
  Gamma = CoefG_I(1)*abs(log10(BetaFit/CoefG_I(2))) + CoefG_I(3)

  dynamic_kappa = A*EnergyFit**Gamma

 end function dynamic_kappa
 !============================================================================
 subroutine cimi_set_boundary_f2( neng,xjac,energy)

  !  *Relativistic version of boundary PSD.*
  ! Set up the boundary distribution for the CIMI. Distribution
  ! at the boundary is assumed to be Maxwellian for ions and a kappa
  ! (kappa=3) distribution for electrons. Boundary temperature and
  ! density are from MHD (or empirical models if Weimer/Tsyganenko fields are
  ! used)

  ! Input: nspec,np,nt,nm,nk,iba,irm,amu,xjac,Den_IC,Temp_IC,vel
  ! Output: fb

  use ModGmCimi, ONLY:  Temp_IC, Temppar_IC
  use ModCimi,        ONLY: MinLon,MaxLon, f2, fb
  use ModCimiTrace,  ONLY: sinA,ekev,iw2,pp, iba, irm, vel, bo
  use ModCimiUtil,     ONLY: ln_gamma
  use ModConst,      ONLY: cMu, cElectronCharge
  use ModCimiGrid, ONLY: nm, nk

  implicit none

  integer, intent(in):: neng
  real, intent(in):: xjac(np,nm,nspec), energy(neng,nspec)

  integer:: j,n,k,m,ib1,iLat
  real:: pi, xmass, chmass
  real:: fbb,fbb1,parE1,perE1,den1,temp32,y2,x2,Psq,&
       Pper2,Ppar2
  real:: erpp,Vsq,Vper2,Vpar2
  real:: BoundaryTempperp_IC(minLon:maxLon,nspec)
  real:: kappa,kappa_plus_one,kappa_minus_half,ln_gamma_diff,gamma_ratio
  real:: e_min, Beta, TempCoreFactor

  !real, external:: ln_gamma
  !----------------------------------------------------------------------------
  pi=acos(-1.)

  !   defining Pperp and Ppar for both cases, anisotropic and isotropic.
  !   Advantages: use the same form of Maxwellian/kappa dist
  !   since in the case of isotropic PSD

  do n=1,nspec
     do j=MinLonPar,MaxLonPar
        BoundaryTempperp_IC(j,n) = &
             (3*BoundaryTemp_IC(j,n) - BoundaryTempPar_IC(j,n))/2
     enddo
  enddo

  fb(MinLonPar:MaxLonPar,1:nm,1:nk,1:nspec)=0.
  do n=1,nspec
     e_min=0.
     if (n==nspec) e_min=0.5*energy(1,n)
     xmass=amu_I(n)*1.673e-27
     chmass=1.6e-19/xmass
     do j=MinLonPar,MaxLonPar
        den1 =BoundaryDens_IC(j,n)                     ! density in m-3
        TempCoreFactor=1.0
        ib1=iba(j)+1   ! difference bewteen cimi standalone where ib1=iba(j)
        if (ib1 > irm(j)) ib1=irm(j)

        ! Find Kappa
        kappa=KappaConst_I(n)
        if(UseDynamicKappa) then
          ! Find species plasma beta
          Beta = den1 * BoundaryTemp_IC(j,n) * cElectronCharge &
            / (bo(ib1,j)**2 / (2 * cMu))
          ! Use Eyelade (2021) fit for kappa as a function of beta + temp
          kappa = dynamic_kappa(Beta, BoundaryTemp_IC(j,n)/1000., &
                                    n, UseKappaTotalFit)
          kappa = min(max(kappa,KappaMin_I(n)),KappaMax_I(n))  ! enforce min/max
        end if
        if(DoRescaleKappaTemp .and. kappa > 0.) &
          TempCoreFactor = (kappa - 1.5)/kappa   ! rescale T to conserve energy
        parE1=TempCoreFactor*BoundaryTempPar_IC(j,n)/1000. ! characteristic in keV
        perE1= TempCoreFactor*BoundaryTempPerp_IC(j,n)/1000.! characteristic in keV
        
        temp32=sqrt(parE1)*perE1*1.6e-16**1.5
        if (kappa > 0.) then    !    kappa, for electrons only
           kappa_plus_one = kappa + 1.
           kappa_minus_half = kappa - 0.5
           ln_gamma_diff = ln_gamma(kappa_plus_one) - ln_gamma(kappa_minus_half)
           gamma_ratio = exp(ln_gamma_diff)
           fbb1=den1*gamma_ratio/temp32/(2.*pi*kappa*xmass)**1.5
        else
           fbb1=den1/temp32/(2.*pi*xmass)**1.5 ! Maxwellian
        endif

        !      do k=1,nm   ! magnetic moment
        do k=1,nk   ! K invariant
           y2=sinA(ib1,j,k)*sinA(ib1,j,k)
           x2=1.-y2
           do m=1,iw2(k,n)   ! magnetic moment
              fbb=0.
              if ((ekev(ib1,j,m,k,n) > e_min) .and. &
                   (BoundaryDens_IC(j,n) > 0.)) then

                 if (kappa > 0.) then   ! KAPPA FOR ELECTRONS
                    Psq=pp(ib1,j,m,k,n)*pp(ib1,j,m,k,n)
                    Ppar2=Psq*x2
                    Pper2=Psq*y2
                    erpp=(Ppar2/parE1+Pper2/perE1)/2./xmass/1.6e-16
                    fbb=fbb1/(1.+erpp/kappa)**kappa_plus_one

                 else ! MAXWELLIAN OTHERWISE:
                    Vsq=vel(ib1,j,m,k,n)*vel(ib1,j,m,k,n)
                    Vpar2=Vsq*x2
                    Vper2=Vsq*y2
                    erpp=(Vpar2/parE1+Vper2/perE1)/2./1000./chmass
                    fbb=0.
                    if (erpp < 500.) fbb=fbb1/exp(erpp)
                 endif
              endif      ! end of ekev(n,ib,j,m,k) > e_min.and.den1 > 0

              fb(j,m,k,n)=fbb*xjac(ib1,m,n)

              do iLat=ib1,np
                 f2(iLat,j,m,k,n)=fbb*xjac(ib1,m,n) ! f2 beyond rb
              enddo

              !        if ((n == nspec).and.(j == 1).and.(m == 10)) write(*,*) fbb,xjac(n,ib1,k), &
              !        ekev(ib1,j,k,m,n), vel(ib1,j,k,m,n),xmass,chmass, &
              !        Psq,Ppar2,Vsq,Vpar2, 'f,jac,e,vel,m,chm,Psq,Ppar,Vsq,Vpar'   !   NB Jan 20 2017
           enddo                ! end of m loop
        enddo                   ! end of k loop
     enddo                      ! end of j loop
  enddo                         ! end of n loop

end subroutine cimi_set_boundary_f2
!==============================================================================
 
end module ModCimiBoundary
