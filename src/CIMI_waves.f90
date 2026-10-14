!modules from standalone cimi as reference
!module constants
!  integer,parameter :: i_zero=0,i_one=1,m_one=-1
!  real,parameter :: pi=3.14159265358979
!  real,parameter :: re_m=6.3712e6	            ! earth's average radius (m)
!  real,parameter :: xmp=1.673e-27               ! mass of H+ in kg
!  real,parameter :: e_mass=9.11e-31             ! electron mass in kg
!  real,parameter :: echarge=1.6e-19	            ! electron charge
!  real,parameter :: EM_speed=2.998e8            ! speed of light (m/s)
!  real,parameter :: epsilon0=8.8542e-12         ! permittivity of free space
!end module constants
!
!
!module cWpower
!  use cimigrid_dim
!  real,parameter::densityP=4.e7   ! density to define plasmapause
!  real rppa(ip)                   ! plasmapause radial distance
!  real CHpower(ir,ip),HIpower(ir,ip),ompe(ir,ip),r_wave
!  real Cpower0,Hpower0,BLc0,BLh0
!end module cWpower

Module CIMI_waves
  !these data used to be taken from the constants module of the standalone cimi
  use ModConst, only: cEps,cLightSpeed, cElectronCharge, cElectronMass, cProtonMAss
  use ModNumConst,      ONLY: cPi
  use ModPlanetConst,	ONLY: Earth_, rPlanet_I
  use ModUtilities,     ONLY: CON_stop
  implicit none
  private !except

  !arrays for the wave power copy and paste directly from cWpower in
  !standalone code
  real,public, allocatable :: CHpower(:,:),HIpower(:,:),ompe(:,:),r_wave
  real,public  ::  Cpower0,Hpower0,BLc0,BLh0
  

  !copy and paste directly from waveDiffCoef in standalone code
  real,allocatable,public,dimension(:) :: cLshell,hLshell,&
                                          cOmpe,hOmpe,&
                                          ckeV,hkeV,&
                                          cPa,hPa
  !chorus diffusion coef
  real,public,allocatable,dimension(:,:,:,:) :: cDaa,cDEE,cDaE!cDpp,cDap,hDaa,hDpp,hDap
  !hiss diffusion coef
  real,public,allocatable,dimension(:,:,:,:) :: hDaa,hDEE,hDaE

  real,allocatable,public,dimension(:,:,:,:) :: E_Q2c0,E_Q2h0,&
                                                cDqq10,cDqq20,hDqq10,hDqq20
  
  !The type of wave model. In the standalone code this was ichor, ihiss
  character(len=100),public :: TypeHiss, TypeChorus

  !Name of input wave files
  character(len=100),public :: NameHissFile, NameChorusFile

  !flags that indicate if we are using waves and what kind
  logical, public :: UseWaves=.false., UseChorus, UseHiss
  logical, public :: UseDiagDiffusion=.true.

  ! dimension of LB chorus & hiss diff coef copied from cimigrid_dim module of
  ! standalone code
  ! dimension of LB chorus & hiss diff coef
  !integer,public :: ipc,iwc,iph,iwh
  !integer, public,parameter :: ipa=89
  !real,public :: cPA(ipa)
  integer,public :: icL,ihL,icfpefce,ihfpefce,icEnergy,ihEnergy,icPa,ihPa

  
  ! option to use KpIndex instead of AE. This was origional wave approach
  logical, public :: UseKpIndex

  public :: readDiffCoef
  public :: WavePower
contains

  !*****************************************************************************
  !                           readDiffCoef
  !  Routine reads chorus and hiss diffusion coeff.
  !  directly adapted from standalone code
  !*****************************************************************************
  subroutine readDiffCoef(ir,MinLon,MaxLon)

    implicit none
    integer, intent(in)::ir,MinLon,MaxLon
    integer :: m

    !allocate arrays
    allocate(CHpower(ir,MinLon:MaxLon),HIpower(ir,MinLon:MaxLon), &
         ompe(ir,MinLon:MaxLon))

    !set the filetype for hiss and chorus based on the file name
    if (UseHiss) then
       select case (NameHissFile)
       case('D_hiss_UCLA.dat')
          TypeHiss='UCLA'
       !case('D_hiss_Albert.dat')
       !   TypeHiss='Albert'
       case('D_hiss_BU.dat')
          TypeHiss='BU'
       case default
          call CON_stop('IM ERROR: Hiss wave file type not recognized.')
       end select
    endif

    if (UseChorus) then
       select case (NameChorusFile)
       case('D_LBchorus_QZ.dat','D_LBchorus_merge2.dat',&
            'D_LBchorus_merge3.dat')
          TypeChorus='UCLA'
       !case('D_LBchorus_merge.dat')
       !   TypeChorus='ChorusFileType2'
       case('D_LBchorus_BUMLAT010.dat',&
            'D_LBchorus_BUMLAT1020.dat',&
            'D_LBchorus_BUMLAT2030.dat')
          TypeChorus='BU'
       case default
          call CON_stop('IM ERROR: Chorus wave file type not recognized.')
       end select
    endif

    if (TypeHiss.ne.TypeChorus) then
        call CON_stop('IM ERROR: Wave file types must match.')
    endif

    
    ! pitch-angle grid
    !do m=1,ipa
    !   cPA(m)=real(m) !m*1.
    !enddo

    !read chorus data
    if (UseChorus) then
       select case(TypeChorus)
       case('UCLA')
         call read_chorus
       case('BU')
         if(UseDiagDiffusion) then
            call read_BU_waves_DQQ(NameChorusFile,&
                               Cpower0,icL,icfpefce,icEnergy,icPa,&
                               cLshell,cOmpe,ckeV,cPa,&
                               E_Q2c0,cDqq10,cDqq20)
         else
            call read_BU_waves_DaE(NameChorusFile,&
                               Cpower0,icL,icfpefce,icEnergy,icPa,&
                               cLshell,cOmpe,ckeV,cPa,&
                               cDaa,cDaE,cDEE)
         endif
       end select
    endif

    if (UseHiss) then
       select case (TypeHiss)
       case('UCLA')
         call read_hiss
       case('BU')
         if(UseDiagDiffusion) then
            call read_BU_waves_DQQ(NameHissFile,&
                               Hpower0,ihL,ihfpefce,ihEnergy,ihPa,&
                               hLshell,hOmpe,hkeV,hPa,&
                               E_Q2h0,hDqq10,hDqq20)
         else
            call read_BU_waves_DaE(NameHissFile,&
                            Hpower0,ihL,ihfpefce,ihEnergy,ihPa,&
                            hLshell,hOmpe,hkeV,hPa,&
                            hDaa,hDaE,hDEE)
         endif
       end select
    endif
  end subroutine readDiffCoef

  !***************************************************************************
  subroutine read_chorus
    use ModIoUnit, ONLY: UnitTmp_
    use ModCimiPlanet, only: xme => dipmom
    !index vars
    integer ::  j,k,m,n

    ! to hold record diffusion coef from file
    real:: pa,daa,dap,dpp,daE,dEE
    
    !reference L location for waves
    real :: Lc0 !
  
    character header*80
    !---------------------------------------------------------------------------
    
    open(unit=UnitTmp_,file='IM/'//NameChorusFile,status='old')
    !write(*,*) 'IM: Using Chorus wave file: ', NameChorusFile
    read(UnitTmp_,*) Cpower0,Lc0
    read(UnitTmp_,*) icfpefce,icEnergy,icPa
    icL = 1
    !note some additional variables will be needed when we go to p-k system
    !those are commented out for now
    allocate (cLshell(icL),&
              cOmpe(icfpefce),&
              ckeV(icEnergy),&
              cPa(icPa))!,ckeVlog(iwc))
    !allocate (wDaa(iwc,ipa),wDap(iwc,ipa),wDpp(iwc,ipa))
    !allocate (Daa1(iwc),Dap1(iwc),Dpp1(iwc))
    !allocate (cDaa(ipc,0:iw+1,ipa),cDap(ipc,0:iw+1,ipa),cDpp(ipc,0:iw+1,ipa))
    allocate (cDEE(icL,icfpefce,icEnergy,icPa),&
              cDaa(icL,icfpefce,icEnergy,icPa),&
              cDaE(icL,icfpefce,icEnergy,icPa))
    read(UnitTmp_,*) cOmpe
    read(UnitTmp_,*) ckeV
    !ckeVlog(:)=log10(ckeV(:))
    BLc0=xme/(Lc0*rPlanet_I(Earth_))**3       ! dipole B at Lc0
    cDaa(icL,:,:,:)=0.
    cDaE(icL,:,:,:)=0.
    cDEE(icL,:,:,:)=0.
    do j=1,icfpefce
       ! Read coefficients
       do k=1,icEnergy
          read(UnitTmp_,'(a80)') header
          read(UnitTmp_,'(a80)') header
          do m=1,icPa
             read(UnitTmp_,*) pa,daa,dap,dpp,daE,dEE
             cPa(m) = pa
             cDaa(icL,j,k,m) = daa    ! coeff in 1/sec
             cDaE(icL,j,k,m) = daE
             cDEE(icL,j,k,m) = dEE
          enddo
       enddo
    enddo
    close(UnitTmp_)
    
  end subroutine read_chorus
  
  !***************************************************************************
  subroutine read_hiss
    use ModIoUnit, ONLY: UnitTmp_
    use ModCimiPlanet, only: xme => dipmom
    !index vars
    integer ::  j,k,m,n

    ! to hold record diffusion coef from file
    real:: pa,daa,dap,dpp,daE,dEE
    
    !reference L location for waves
    real :: Lh0 , EEo
  
    character header*80
    !---------------------------------------------------------------------------
    
    open(unit=UnitTmp_,file='IM/'//NameHissFile,status='old')
    !write(*,*) 'IM: Using Hiss wave file: ', NameHissFile
    read(UnitTmp_,*) Hpower0,Lh0
    read(UnitTmp_,*) ihfpefce,ihEnergy,ihPa
    ihL = 1
    !note some additional variables will be needed when we go to p-k system
    !those are commented out for now
    allocate (hLshell(ihL),&
              hOmpe(ihfpefce),&
              hkeV(ihEnergy),&
              hPa(ihPa))!,ckeVlog(iwc))
    !allocate (wDaa(iwc,ipa),wDap(iwc,ipa),wDpp(iwc,ipa))
    !allocate (Daa1(iwc),Dap1(iwc),Dpp1(iwc))
    !allocate (cDaa(ipc,0:iw+1,ipa),cDap(ipc,0:iw+1,ipa),cDpp(ipc,0:iw+1,ipa))
    allocate (hDEE(ihL,ihfpefce,ihEnergy,ihPa),&
              hDaa(ihL,ihfpefce,ihEnergy,ihPa),&
              hDaE(ihL,ihfpefce,ihEnergy,ihPa))
    read(UnitTmp_,*) hOmpe
    read(UnitTmp_,*) hkeV
    !ckeVlog(:)=log10(ckeV(:))
    BLh0=xme/(Lh0*rPlanet_I(Earth_))**3       ! dipole B at Lc0
    hDaa(ihL,:,:,:)=0.
    hDaE(ihL,:,:,:)=0.
    hDEE(ihL,:,:,:)=0.
    do j=1,ihfpefce
       ! Read coefficients
       read(UnitTmp_,'(a80)') header
       read(UnitTmp_,'(a80)') header
       do k=1,ihEnergy
          do m=1,ihPa
             read(UnitTmp_,*) Lh0,EEo,pa,daa,dap,dpp,daE,dEE
             hPa(m) = pa
             hDaa(ihL,j,k,m) = daa    ! coeff in 1/sec
             hDaE(ihL,j,k,m) = daE
             hDEE(ihL,j,k,m) = dEE
          enddo
       enddo
    enddo
    close(UnitTmp_)
    
  end subroutine read_hiss

  !***************************************************************************
  subroutine read_BU_waves_DQQ(NameFile,&
                               power0,iL,ifpefce,iEnergy,iPa,&
                               Lshell,Ompe,Energy,Pa,&
                               EQ2,Dq1q1,Dq2q2)
    use ModIoUnit, ONLY: UnitTmp_

    ! Inputs
    character,intent(in):: NameFile
    ! Outputs
    real,intent(out):: power0
    integer,intent(out):: iL,ifpefce,iEnergy,iPa
    real, allocatable, intent(out), dimension(:) :: Lshell,Ompe,Energy,Pa
    real, allocatable, intent(out), dimension(:,:,:,:) :: EQ2,Dq1q1,Dq2q2
    ! index vars
    integer :: i, j, k, m

    ! dummy variables for reading lines
    real :: raw_daa, raw_dap, raw_dpp
    real :: eq2_val, raw_dq1q1, raw_dq2q2
    real :: blk_E, blk_fpe, blk_L

    character(len=100) :: header
    !---------------------------------------------------------------------------

    open(unit=UnitTmp_, file='IM/'//NameHissFile, status='old')

    ! Read nominal wave power
    read(UnitTmp_, *) power0, header

    ! Read grid sizes: nLshell, nfpefce, nEnergy/EQ2, nPa/alpha0
    read(UnitTmp_, *) iL, ifpefce, iEnergy, iPa, header

    ! Allocate 1D arrays for the grids
    allocate (Lshell(iL))
    allocate (Ompe(ifpefce), Energy(iEnergy), Pa(iPa))

    ! Allocate 4D arrays for Q-coordinate coefficients (Lshell, fpefce, EQ2, alpha0)
    allocate (EQ2(iL, ifpefce, iEnergy, iPa))
    allocate (Dq1q1(iL, ifpefce, iEnergy, iPa))
    allocate (Dq2q2(iL, ifpefce, iEnergy, iPa))

    Dq1q1(:,:,:,:) = 0.
    Dq2q2(:,:,:,:) = 0.

    ! Read 1D grid arrays (skip text headers)
    read(UnitTmp_, '(a80)') header  ! "L-shell"
    read(UnitTmp_, *) Lshell

    read(UnitTmp_, '(a80)') header  ! "fpe/fce"
    read(UnitTmp_, *) Ompe

    read(UnitTmp_, '(a80)') header  ! "energy (keV) / Q2 (keV)"
    read(UnitTmp_, *) Energy

    read(UnitTmp_, '(a80)') header  ! "pitch-angle (deg) / Q1 (deg)"
    read(UnitTmp_, *) Pa

    ! Loop over blocks (Lshell -> fpe/fce -> Energy/EQ2)
    do i = 1, iL
       do j = 1, ifpefce
          do k = 1, iEnergy

             ! Read block header
             read(UnitTmp_, *) blk_E, blk_fpe, blk_L

             ! Skip the column header line
             read(UnitTmp_, '(a80)') header

             do m = 1, iPa
                ! Column 5: Eq2(keV) -> eq2_val
                ! Column 6: Dq1q1    -> raw_dq1q1
                ! Column 7: Dq2q2/Q2^2 (sec^-1 -> raw_dq2q2
                read(UnitTmp_, *) &
                                  Pa(m), raw_daa, raw_dap, raw_dpp, &
                                  eq2_val, raw_dq1q1, raw_dq2q2
               ! Store the local EQ2 coordinate for this (L, fpe/fce, E, Pa)
                EQ2(i, j, k, m) = eq2_val
                ! Store Dq1q1, and Dq2q2 directly
                Dq1q1(i, j, k, m) = raw_dq1q1
                Dq2q2(i, j, k, m) = raw_dq2q2
             enddo
          enddo
       enddo
    enddo

    close(UnitTmp_)
  end subroutine read_BU_waves_DQQ

  !***************************************************************************
  subroutine read_BU_waves_DaE(NameFile,&
                           power0,&
                           iL,ifpefce,iEnergy,iPa,&
                           Lshell,Ompe,keV,Pa,&
                           Daa,DaE,DEE)
    use ModIoUnit, ONLY: UnitTmp_
    use ModConst,  ONLY: cRmeElectron,cKEV

    ! Input variables
    character(len=100),intent(in) :: NameFile
    ! Output variables
    real,intent(out) :: power0
    integer,intent(out) :: iL,ifpefce,iEnergy,iPa
    real,allocatable,dimension(:),intent(out) :: Lshell,Ompe,keV,Pa
    real,allocatable,dimension(:,:,:,:),intent(out) :: Daa,DaE,DEE
    ! index variables
    integer :: i, j, k, m 
    ! to hold record diffusion coef from file
    real :: dum_pa, raw_daa, raw_dap, raw_dpp
    real :: dum_Eq2, dum_Dq1, dum_Dq2
    real :: blk_E, blk_fpe, blk_L

    ! conversion variables
    real :: E_keV, cE_rest, gamma, v_c

    character(len=100) :: header
    !-------------------------------------------------------------------------

    ! Rest mass of electron in keV (used for momentum to energy conversion)
    cE_rest = cRmeElectron/cKEV

    open(unit=UnitTmp_, file='IM/'//NameFile, status='old')
    !write(*,*) 'IM: Using Hiss wave file: ', NameFile
    ! Read nominal wave power
    read(UnitTmp_, *) power0
    ! Read grid sizes: nLshell, nfpefce, nEnergy, nPa
    read(UnitTmp_, *) iL, ifpefce, iEnergy, iPa
    ! Allocate 1D arrays for the grids
    allocate (Lshell(iL),&
              Ompe(ifpefce),&
              keV(iEnergy),&
              Pa(iPa))
    ! Allocate 4D arrays for coefficients (Lshell, fpe/fce, Energy, PitchAngle)
    allocate (Daa(iL, ifpefce, iEnergy, iPa))
    allocate (DaE(iL, ifpefce, iEnergy, iPa))
    allocate (DEE(iL, ifpefce, iEnergy, iPa))
    Daa(:,:,:,:) = 0.
    DaE(:,:,:,:) = 0.
    DEE(:,:,:,:) = 0.
    ! Read 1D grid arrays (skip text headers)
    read(UnitTmp_, '(a80)') header  ! "L-shell"
    read(UnitTmp_, *) Lshell

    read(UnitTmp_, '(a80)') header  ! "fpe/fce"
    read(UnitTmp_, *) Ompe

    read(UnitTmp_, '(a80)') header  ! "energy (keV) / Q2 (keV)"
    read(UnitTmp_, *) keV

    read(UnitTmp_, '(a80)') header  ! "pitch-angle (deg) / Q1 (deg)"
    read(UnitTmp_, *) Pa

    ! Loop over blocks (Block order: Lshell -> fpe/fce -> Energy)
    do i = 1, iL
       do j = 1, ifpefce
          do k = 1, iEnergy
             ! Read block header
             read(UnitTmp_, *) blk_E, blk_fpe, blk_L
             ! Skip the column header line
             read(UnitTmp_, '(a80)') header
             ! Current energy in keV for conversion
             E_keV = keV(k)
             ! Relativistic conversion: velocity (v/c)
             gamma = (E_keV / cE_rest) + 1.0
             v_c = sqrt(gamma**2 - 1.0) / gamma
             do m = 1, iPa
                ! Read the 7 data columns
                read(UnitTmp_, *) Pa(m), raw_daa, raw_dap, raw_dpp, &
                                  dum_Eq2, dum_Dq1, dum_Dq2
                ! --- UNIT & PHYSICS CONVERSION ---
                ! <Daa> is already the correct units (Daa/p^2)
                Daa(i,j,k,m) = raw_daa
                ! <DaE> (1/sec) from Dap/E^2:
                ! D_aE = v * D_ap -> <DaE> rate = v_c * E_keV * (Dap / E^2)
                DaE(i,j,k,m) = v_c * E_keV * raw_dap
                ! <DEE> (1/sec) from Dpp/E^2:
                ! D_EE = v^2 * D_pp -> <DEE> rate = (v_c**2) * (Dpp / E^2)
                DEE(i,j,k,m) = (v_c**2) * raw_dpp
             enddo
          enddo
       enddo
    enddo

    close(UnitTmp_)

  end subroutine read_BU_waves_DaE
      
  
  !****************************************************************************
  subroutine WavePower(ir,t,density,AE,Kp,Bo,ro,xmlto,iba,MinLon,MaxLon,MinLonPar,MaxLonPar)
    ! Routine determines Chorus and hiss wave power (pT)^2 and fpe/fce as a
    ! function of location
    
    implicit none
    integer, intent(in)::ir,MinLonPar, MaxLonPar, MinLon, MaxLon
    integer,intent(in) :: iba(MinLon:MaxLon)
    real,   intent(in) :: t, AE, Kp,ro(ir,MinLon:MaxLon), density(ir,MinLonPar:MaxLonPar), Bo(ir,MinLon:MaxLon), &
         xmlto(ir,MinLon:MaxLon)
    real :: ro1,xmlt1,r_hiss, AE_chorus, AE_hiss
    integer :: iae,jae,i,j
    
    r_wave=2.0               ! lower boundary in RE of wave diffusion
    r_hiss=8.                ! upper boundary in RE of hiss
    
    
    if (UseKpIndex) then
       if (Kp.lt.2.) then
          iae=1
          jae=1
          AE_chorus = 0
          AE_hiss = 0
       endif
       if (Kp.ge.2..and.Kp.lt.4.) then
          iae=2
          jae=2
          AE_chorus = 100
          AE_hiss = 100
       endif
       if (Kp.ge.4.) then
          iae=3
          jae=3
          AE_chorus = 300
          AE_hiss = 500
       endif
    else
       ! Determine AE level in chorus wave power data provided by Meredith
       if (AE.lt.100.) iae=1
       if (AE.ge.100..and.AE.lt.300.) iae=2
       if (AE.ge.300.) iae=3
       
       ! Determine AE level in hiss wave power data provided by Meredith
       if (AE.lt.100.) jae=1
       if (AE.ge.100..and.AE.lt.500.) jae=2
       if (AE.ge.500.) jae=3
       AE_chorus = AE
       AE_hiss = AE
    endif
    
    ! Determine wave power (in pT^2), and ompe (fpe/fce)
    CHpower=0.
    HIpower=0.
    ompe=0.            ! initial values
    do j=MinLonPar,MaxLonPar
       do i=1,iba(j)
          if (density(i,j).eq.0.) then
             write(*,*) 'Error: density(i,j).eq.0, t,iba(j),i,j ',t,iba(j),i,j
             stop
          endif
          !write(*,*)'i,j,density(i,j),bo(i,j)',i,j,density(i,j),bo(i,j)
          ompe(i,j)=sqrt(density(i,j)*cElectronMass/cEps)/bo(i,j)    ! fpe/fce
          ro1=ro(i,j)
          xmlt1=xmlto(i,j)
          if (xmlt1.lt.0.) xmlt1=xmlt1+24.
          if (UseChorus .and.ro1.ge.r_wave) call ChorusBpower(ro1,xmlt1,iae, &
               CHpower(i,j), AE_chorus)
          if (UseHiss.and.ro1.ge.r_wave.and.ro1.le.r_hiss) &
               call HissBpower(ro1,xmlt1,jae,HIpower(i,j), AE_hiss)
       enddo
    enddo
    
  end subroutine WavePower
  
  
  !*****************************************************************************
  subroutine ChorusBpower(ro,xmlt,iae,Bpower, AE)
    ! Routine calculates the lower band Chorus power in (pT)^2 by Gaussian fits
    ! provided by Qiuhua.
    
    implicit none

    real,    intent(in) :: ro,xmlt, AE
    integer, intent(in) :: iae
    real   , intent(out):: Bpower
    
    real    :: xxbar,yybar,pi2,r12,sigx1,sigy1,qq,amp_fit, r_xy, x_bar, y_bar,&
               sig_x, sig_y
    real    :: Amp(3),xbar(3),sigx(3),ybar(3),sigy(3),rxy(3)
    !--------------------------------------------------------------------------
    data Amp/18951.,54237.,111104./  
    data xbar/7.46,7.82,5.20/        ! mlt of peak power
    data ybar/6.44,6.37,5.70/        ! ro of peak power
    data sigx/5.50,4.87,4.64/
    data sigy/1.30,1.68,1.77/
    data rxy/0.,-0.1,0./

    if (AE < 300) then
       amp_fit = 23245.615 * sqrt(0.249696 * AE + 55.73685) - 154594.057
       r_xy = 0.000005 * (AE ** 2) - (0.0015 * AE)
       x_bar = -0.0000556667 * (AE ** 2) + (0.00916667 * AE) + 7.46
       y_bar = 2.9966 + (6.44-2.9966) / (1 + (AE/520.915) ** 2.3479)
    else
       amp_fit = 111104.
       r_xy = 0.0
       x_bar = 5.20
       y_bar = 5.70
    end if
    if (AE < 1) then
       sig_x = 5.50734
       sig_y = 1.30006
    else
       sig_x = -0.14667 * log(AE) + 5.50734
       sig_y = 0.0824351 * log(AE) + 1.30006
    end if
    
    pi2=2.*cPi        
    r12=1./sqrt(1.-r_xy*r_xy)
    sigx1=sig_x
    sigy1=sig_y
    
    xxbar=abs(xmlt-x_bar)
    if (xxbar.gt.12.) xxbar=24.-xxbar
    yybar=abs(ro-y_bar)
    qq=xxbar**2/sigx1**2+yybar**2/sigy1**2-2.*r_xy*xxbar*yybar/sigx1/sigy1
    Bpower=amp_fit*r12/pi2/sigx1/sigy1/exp(0.5*r12*qq)

    !kludge
    !Bpower=0.0
    !Bpower=1.0e-40
    
  end subroutine ChorusBpower

  
  !*****************************************************************************
  subroutine HissBpower(ro,xmlt,iae,Bpower, AE)
    ! Routine calculates the CRRES hiss power in (pT)^2 by Gaussian fits
    ! provided by Qiuhua.
    
    implicit none
    
    real,    intent(in) :: ro,xmlt, AE
    integer, intent(in) :: iae
    real   , intent(out):: Bpower
    
    real xxbar,yybar,pi2,r12,sigx1,sigy1,qq, amp_fit,r_xy, x_bar, y_bar,&
               sig_x, sig_y
    real Amp(3),xbar(3),sigx(3),ybar(3),sigy(3),rxy(3)
    
    data Amp/37964.,71459.,130742./          
    data xbar/15.20,11.82,14.46/
    data ybar/3.0,3.25,3.55/
    data sigx/6.08,4.87,5.64/
    data sigy/0.81,1.31,1.32/
    data rxy/0.,0.,0./

    if (AE < 500) then
       amp_fit = 899.89849 * sqrt(30.45545 * AE + 497.35063) + 17895.0403
       r_xy = 0.0
       x_bar = 0.0000808 * (AE ** 2) - (0.04188 * AE) + 15.20
       y_bar = 0.00793319 * sqrt(9.60314 * AE) + 3.00148
       sig_x = 0.00002805 * (AE ** 2) - (0.014905 * AE) + 6.08
       sig_y = -0.00000995 * (AE ** 2) + (0.005995 * AE) + 0.81
    else
       amp_fit = 130742.
       r_xy = 0.0
       x_bar = 14.46
       sig_x = 5.64
       sig_y = 1.32
    end if
    

    pi2=2.*cPi          
    r12=1./sqrt(1.-r_xy*r_xy)
    sigx1=sig_x
    sigy1=sig_y
    
    xxbar=abs(xmlt-x_bar)
    if (xxbar.gt.12.) xxbar=24.-xxbar 
    yybar=abs(ro-y_bar)
    qq=xxbar**2/sigx1**2+yybar**2/sigy1**2-2.*r_xy*xxbar*yybar/sigx1/sigy1
    Bpower=amp_fit*r12/pi2/sigx1/sigy1/exp(0.5*r12*qq)


    !kludge
    !Bpower=0.0
    !Bpower=1.0e-40

  end subroutine HissBpower
  
  !note that other routines like diffuse_aa and diffuse_ee can be added later
  !from the standalone code
end Module CIMI_waves







  
