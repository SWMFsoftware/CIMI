Module DensityTemp
use ModCrcmGrid,ONLY: ir=>np, ip=>nt
implicit none
real,parameter ::  density(ir,ip)=1.e+6
EndModule


!!!!!!!!!!!!!!!!  Main wave module !!!!!!!!!!!!!1
Module ModWaveDiff

use ModCrcmGrid,ONLY: ir=>np, ip=>nt, iw=>nm , ik=>nk  
use ModCrcmPlanet,ONLY: nspec
use ModFieldTrace, ONLY: ekev,tya,y=>sinA,Bo,ro,xmlto,lintpIM,bm
use ModNumConst,       ONLY: pi => cPi


implicit none

logical            :: UseWaveDiffusion = .false.
logical            :: UseHiss  = .false.
logical            :: UseChorus  = .false.
logical            :: UseChorusUB = .false.
logical            :: testDiff_aa = .false.
logical            :: testDiff_EE = .false.      
real               :: DiffStartT = 0.
character(len=100)  :: HissWavesD = ' '
character(len=100)  :: ChorusWavesD= ' '
character(len=100)  :: ChorusUpperBandD = ' '


integer,parameter :: ipu=6,iwu=31,ipa=89    ! from module cimigrid_dim in CIMI: dimension of Qiuhua's UB chorus coef
integer ipc,iwc,iph,iwh            ! from module cimigrid_dim in CIMI: dimension of LB chorus & hiss diff coef
integer ichor,iUBC,ihiss           !   switches for three wave models    
!integer iba(ip),iw2(nspec,ik)

! from module WaveDiffCoef in cimi
real, allocatable, dimension(:) :: cOmpe,ckeV,hOmpe,hkeV
real, allocatable, dimension(:,:,:) :: cDEE,cDaa,cDaE,hDEE,hDaa,hDaE
real cPA(ipa),uOmpe(ipu),ukeV(iwu), &
     uDEE(ipu,iwu,ipa),uDaa(ipu,iwu,ipa),uDaE(ipu,iwu,ipa)

!from module cWpower in cimi
real CHpower(ir,ip),UBCpower(ir,ip),HIpower(ir,ip),ompe(ir,ip),r_wave
real Cpower0,Hpower0,BLc0,BLh0,BLu0

! from module constants
real,parameter :: re_m=6.375e6                ! earth's radius (m)
real,parameter :: e_mass=9.11e-31             ! electron mass in kg
real,parameter :: epsilon0=8.8542e-12         ! permittivity of free space


contains

  subroutine ReadDiffCoef(rc)

!  use constants  check later
!  use cimigrid_dim
!  use cread1
!  use waveDiffCoef
!  use cWpower

  use ModCrcmPlanet,  ONLY: xme=>dipmom
  use ModIoUnit, ONLY: UnitTmp_

  implicit none

 integer j,k,m  
 real Eo,c_cgs,pa,daa,dap,dpp,daE,dEE,Dnorm,cMeV,Lc0,Lh0
 real E_cgs,cE2,EEo,E2Eo,daan,dapn,dppn,daa0,dap0,dpp0
 real,parameter :: Lu0=6.
 character header*80
 real rc

! list of available wave diffusion coefficients models:
!  main lower band chorus model:  D_LBchorus_merge.dat
!  Q.Zheng lower band model:      D_LBchorus_QZ.dat
!  main hiss model:               D_hiss_UCLA.dat
!  additionall hiss model:        D_hiss_Albert.dat


 ihiss=0
 if (UseHiss) then
   if ( trim(HissWavesD) .eq.'D_hiss_UCLA.dat') ihiss=2
 else 
   ihiss=1
 endif

 ichor=0
 if (UseChorus) then 
   if ( trim(ChorusWavesD) .eq.'D_LBchorus_merge.dat') ichor=2
 else 
  ichor=1
 endif

 iUBC=0
 if (UseChorusUB) iUBC=1

 write(*,*) '*** WAVE MODEL FOR ELECTRONS IS ON ***'
 write(*,*) 'UseHiss:            ', UseHiss
 write(*,*) 'Lower (main) Chorus:', UseChorus
 write(*,*) 'Upper  Chorus:      ', UseChorusUB
 write(*,*) '***'

  Eo=511.                 ! electron rest energy in keV
  c_cgs=2.998e10          ! speed of light in cgs

! pitch-angle grid
  do m=1,ipa
     cPA(m)=m*1.
  enddo

! Ompe grid. Ratio of fpe/fpc and energy for upper-band chorus
  uOmpe(1:ipu)=(/1.,3.,6.,9.,12.,20./)
  ukeV(1:iwu)=(/0.1,0.1259,0.1585,0.1995,0.2512,0.3162,0.3981,0.5012, &
                0.631,0.7943,1.0,1.259,1.585,1.995,2.512,3.162, &
                3.981,5.012,6.31,7.943,10.0,12.59,15.85,19.95, &
                25.12,31.62,39.81,50.12,63.1,79.43,100.0/)


  if (UseChorus) then
! Read LB chorus ckeV(0.1keV-10MeV), Daa, DEE, DaE at L = 6.5
  write(*,*) 'Reading chorus data, ichor=',ichor
  write(*,*) 'Chorus Dcoef are from: ', trim(ChorusWavesD)
  open(unit=UnitTmp_,file='IM/'//trim(ChorusWavesD),status='old')
!  if (ichor.eq.1) open(unit=UnitTmp_,file='IM/D_LBchorus_QZ.dat',status='old')
!  if (ichor.eq.2) open(unit=UnitTmp_,file='IM/D_LBchorus_merge.dat',status='old')
  read(UnitTmp_,*) Cpower0,Lc0
  read(UnitTmp_,*) ipc,iwc
  allocate (cOmpe(ipc),ckeV(iwc))
  allocate (cDaa(ipc,iwc,ipa),cDaE(ipc,iwc,ipa),cDEE(ipc,iwc,ipa))
  read(UnitTmp_,*) cOmpe
  read(UnitTmp_,*) ckeV
  BLc0=xme/(Lc0*re_m)**3       ! dipole B at Lc0
  cDaa(:,:,:)=0.
  cDaE(:,:,:)=0.
  cDEE(:,:,:)=0.
  if (ichor.ge.1) then
     do j=1,ipc
        do k=1,iwc
           read(UnitTmp_,'(a80)') header
           read(UnitTmp_,'(a80)') header
           do m=1,ipa
              if (ichor.eq.1) read(UnitTmp_,*) pa,daa,dap,dpp,daE,dEE
              if (ichor.eq.2) read(UnitTmp_,*) pa,daa,daE,dEE
              cDaa(j,k,m)=daa    ! coeff in 1/sec
              cDaE(j,k,m)=daE
              cDEE(j,k,m)=dEE
           enddo
        enddo
     enddo
  endif
  close(UnitTmp_)
  endif ! if UseChorus

! Read UB chorus Daa, DEE, DaE at L = 6.0
  BLu0=xme/(Lu0*re_m)**3       ! dipole B at Lu0
  uDaa(:,:,:)=0.
  uDaE(:,:,:)=0.
  uDEE(:,:,:)=0.

  if (UseChorusUB) then
    write(*,*) 'Reading UBchorus, iUBC=',iUBC
    write(*,*) 'Upper band Dcoeff are from: ', trim(ChorusUpperBandD)
    open(unit=UnitTmp_,file='IM/'//trim(ChorusUpperBandD),status='old')
!    open(unit=UnitTmp_,file='IM/D_UBchorus.dat',status='old')
  if (iUBC.eq.1) then
      do j=1,6
         read(UnitTmp_,'(a80)') header
      enddo
      do j=1,ipu
         do k=1,iwu
            read(UnitTmp_,'(a80)') header
            do m=2,ipa
               read(UnitTmp_,*) pa,uDaa(j,k,m),uDaE(j,k,m),uDEE(j,k,m) ! D in 1/s
            enddo
         enddo
      enddo
      close(UnitTmp_)
  endif
  endif ! if UseChorusUB

! Read hiss Daa, DEE, DaE at L = 5.5


  if (UseHiss) then
  write(*,*) 'Reading hiss data, ihiss=',ihiss
  write(*,*) 'Hiss Dcoeff are from: ',trim(HissWavesD)
  open(unit=UnitTmp_,file='IM/'//trim(HissWavesD),status='old')
  !if (ihiss.eq.2) open(unit=UnitTmp_,file='IM/D_hiss_UCLA.dat',status='old')
  !if (ihiss.eq.1) open(unit=UnitTmp_,file='IM/D_hiss_Albert.dat',status='old')
  read(UnitTmp_,*) Hpower0,Lh0
  read(UnitTmp_,*) iph,iwh
  allocate (hOmpe(iph),hkeV(iwh))
  allocate (hDaa(iph,iwh,ipa),hDaE(iph,iwh,ipa),hDEE(iph,iwh,ipa))
  hDaa(:,:,:)=0.
  hDaE(:,:,:)=0.
  hDEE(:,:,:)=0.
  read(UnitTmp_,*) hOmpe
  read(UnitTmp_,*) hkeV
  BLh0=xme/(Lh0*re_m)**3       ! dipole B at Lh0
  do j=1,iph
     if (ihiss.eq.1) then
        do k=1,iwh
           read(UnitTmp_,'(56x,e10.4)') Dnorm
           read(UnitTmp_,'(a80)') header
           E_cgs=hkeV(k)*1.6e-9        ! energy in erg
           cE2=(c_cgs/E_cgs)**2
           EEo=hkeV(k)+Eo
           E2Eo=hkeV(k)+2.*Eo
           do m=1,ipa
              read(UnitTmp_,*) pa,daan,dapn,dppn,daa0,dap0,dpp0
              daa=Dnorm*(daan+daa0)    ! coeff in p^2/s
              dap=Dnorm*(dapn+dap0)    !
              dpp=Dnorm*(dppn+dpp0)    !
              hDaa(j,k,m)=daa*cE2*hkeV(k)/E2Eo           ! coeff in s-1
              hDaE(j,k,m)=dap*cE2*hkeV(k)/EEo            !
              hDEE(j,k,m)=dpp*cE2*hkeV(k)*E2Eo/EEo/EEo   !
           enddo
        enddo
     endif   ! endif (ihiss.eq.1)
     if (ihiss.eq.2) then
        read(UnitTmp_,'(a80)') header
        read(UnitTmp_,'(a80)') header
        do k=1,iwh
           do m=1,ipa
              read(UnitTmp_,*) Lh0,EEo,pa,hDaa(j,k,m),dpp,dap,hDaE(j,k,m),hDEE(j,k,m)
           enddo
        enddo
     endif   ! endif (ihiss.eq.2)
  enddo
  close(UnitTmp_)
  endif ! UseHiss

   end subroutine ReadDiffCoef



!*******************************************************************************
  subroutine WavePower(t,AE,iba)
!*******************************************************************************
! Routine determines Chorus and hiss wave power (pT)^2 and fpe/fce as a function
! of location
!
! Input: t,density,AE,Bo,ro,xmlto,iba
! Output: CHpower,UBCpower,HIpower,ompe,r_wave (through module cWpower)

 ! use constants
 ! use cWpower
  use ModCrcmGrid,   ONLY: MinLonPar,MaxLonPar
  use DensityTemp, ONLY: density
  implicit none
  integer iba(ip),iae,jae,i,j
  real t,AE,ro1,xmlt1
  real r_hiss

  r_wave=2.5               ! lower boundary in RE of wave diffusion
  r_hiss=8.                ! upper boundary in RE of hiss

! Determine AE level in chorus wave power data provided by Meredith
  if (AE.lt.100.) iae=1
  if (AE.ge.100..and.AE.lt.300.) iae=2
  if (AE.ge.300.) iae=3

! Determine AE level in hiss wave power data provided by Meredith
  if (AE.lt.100.) jae=1
  if (AE.ge.100..and.AE.lt.500.) jae=2
  if (AE.ge.500.) jae=3

! Determine wave power (in pT^2), and ompe (fpe/fce)
  CHpower=0.
  UBCpower=0.
  HIpower=0.
  ompe=0.            ! initial values
!  do j=1,ip
  do j=MinLonPar,MaxLonPar
     do i=1,iba(j)
        if (density(i,j).eq.0.) then
           write(*,*) 'Error: density(i,j).eq.0, t,iba(j),i,j ',t,iba(j),i,j
           call CON_STOP('CRCM dies in WavePower')
        endif
        ompe(i,j)=sqrt(density(i,j)*e_mass/epsilon0)/bo(i,j)    ! fpe/fce
        ro1=ro(i,j)
        xmlt1=xmlto(i,j)
        if (xmlt1.lt.0.) xmlt1=xmlt1+24.
        if (ichor.ge.1.and.ro1.ge.r_wave) call ChorusBpower(ro1,xmlt1,iae, &
                                                              CHpower(i,j))
        if (iUBC.eq.1.and.ro1.ge.r_wave) call UBCBpower(ro1,xmlt1,iae, &
                                                          UBCpower(i,j))
        if (ihiss.ge.1.and.ro1.ge.r_wave.and.ro1.le.r_hiss) &
             call HissBpower(ro1,xmlt1,jae,HIpower(i,j))
     enddo
  enddo

  end subroutine WavePower

  !*****************************************************************************
  subroutine ChorusBpower(ro,xmlt,iae,Bpower)
!*****************************************************************************'
! Routine calculates the lower band Chorus power in (pT)^2 by Gaussian fits
! provided by Qiuhua.

!  use constants
  implicit none

  integer iae
  real Bpower,ro,xmlt,xxbar,yybar,pi2,r12,sigx1,sigy1,qq
  real Amp(3),xbar(3),sigx(3),ybar(3),sigy(3),rxy(3)

  data Amp/18951.,54237.,111104./
  data xbar/7.46,7.82,5.20/
  data ybar/6.94,6.87,6.20/
  data sigx/5.50,4.87,4.64/
  data sigy/0.80,1.18,1.27/
  data rxy/0.,-0.1,0./

  pi2=2.*pi
  r12=1./sqrt(1.-rxy(iae)*rxy(iae))
  sigx1=sigx(iae)
  sigy1=sigy(iae)

  xxbar=abs(xmlt-xbar(iae))
  if (xxbar.gt.12.) xxbar=24.-xxbar
  yybar=abs(ro-ybar(iae))
  qq=xxbar**2/sigx1**2+yybar**2/sigy1**2-2.*rxy(iae)*xxbar*yybar/sigx1/sigy1
  Bpower=Amp(iae)*r12/pi2/sigx1/sigy1/exp(0.5*r12*qq)

  end subroutine ChorusBpower

  !*****************************************************************************
  subroutine UBCBpower(ro,xmlt,iae,Bpower)
!*****************************************************************************'
! Routine calculates the upper band Chorus power in (pT)^2 by Gaussian fits
! provided by Qiuhua.

!  use constants
  implicit none

  integer iae
  real Bpower,ro,xmlt,xxbar,yybar,pi2,r12,sigx1,sigy1,qq
  real Amp(3),xbar(3),sigx(3),ybar(3),sigy(3),rxy(3)

  data Amp/225.5,480.,1440./
  data xbar/9.,6.,4./
  data ybar/5.5,5.,5./
  data sigx/3.4,3.5,4./
  data sigy/0.75,0.75,0.75/
  data rxy/0.05,0.,0.1/

  pi2=2.*pi
  r12=1./sqrt(1.-rxy(iae)*rxy(iae))
  sigx1=sigx(iae)
  sigy1=sigy(iae)

  xxbar=abs(xmlt-xbar(iae))
  if (xxbar.gt.12.) xxbar=24.-xxbar
  yybar=abs(ro-ybar(iae))
  qq=xxbar**2/sigx1**2+yybar**2/sigy1**2-2.*rxy(iae)*xxbar*yybar/sigx1/sigy1
  Bpower=Amp(iae)*r12/pi2/sigx1/sigy1/exp(0.5*r12*qq)

  end subroutine UBCBpower

   !*****************************************************************************
  subroutine HissBpower(ro,xmlt,iae,Bpower)
!*****************************************************************************'
! Routine calculates the CRRES hiss power in (pT)^2 by Gaussian fits
! provided by Qiuhua.

!  use constants
  implicit none

  integer iae
  real Bpower,ro,xmlt,xxbar,yybar,pi2,r12,sigx1,sigy1,qq
  real Amp(3),xbar(3),sigx(3),ybar(3),sigy(3),rxy(3)

  data Amp/37964.,71459.,130742./
  data xbar/15.20,11.82,14.46/
  data ybar/3.0,3.25,3.55/
  data sigx/6.08,4.87,5.64/
  data sigy/0.81,1.31,1.32/
  data rxy/0.,0.,0./

  pi2=2.*pi
  r12=1./sqrt(1.-rxy(iae)*rxy(iae))
  sigx1=sigx(iae)
  sigy1=sigy(iae)

  xxbar=abs(xmlt-xbar(iae))
  if (xxbar.gt.12.) xxbar=24.-xxbar
  yybar=abs(ro-ybar(iae))
  qq=xxbar**2/sigx1**2+yybar**2/sigy1**2-2.*rxy(iae)*xxbar*yybar/sigx1/sigy1
  Bpower=Amp(iae)*r12/pi2/sigx1/sigy1/exp(0.5*r12*qq)

  end subroutine HissBpower


!!!!!!!!!!!!!!!!!!!   main  subroutine  for diffusion in pitch-angle !!!!!!!!!!!!

   subroutine diffuse_aa(f2,dt,xjac,iba,iw2)
      use DensityTemp
      use ModMPI
      use ModCrcmGrid,   ONLY: MinLonPar,MaxLonPar
      implicit none
      integer,parameter :: ie=40 ! diffusion subroutine has its own energy grid
      integer i,j,m,k,irun,n,nn,ier
      integer iba(ip),iw2(nspec,ik)
      real f2(nspec,ir,ip,iw,ik),xjac(nspec,ir,iw),dt  ! xjac has different dimension from CIMI !
           

      real ein_log(ie),einlog, &
           f1d(iw),e1d(iw),ein(ie),ao(0:ik+1),dao(ik),Gjac(0:ik+1),DD(0:ik+1),&
           um(ik),up(ik),a1d(ik),b1d(ik),c1d(ik),f0(0:ik+1),uDaoao, &
           f2d(ie,ik),df(ie,ik),df1(ie),fr(ik),ekevlog(iw,ik), &
           hDaoao,Cfactor,Upower0,Hfactor,Daoao,Ufactor, &
           ckeV_log(iwc),ukeV_log(iwu),hkeV_log(iwh)
      real u_max,u_max_log,ompe1,ro1,emin,emax,x0,x2,x,cDaoao
      real u_mx,u_mx_log,factor1,factor_1,DDm,DDp,ump_mx,xlam,alam,dpsd,ao_d
      real densityP,edmin,edmax


     u_max=5.              ! magic number for numerical method from M.C. Fok
     u_max_log=log10(u_max)
     Upower0=10000.        ! coeff based on UB chorus power of (100pT)^2
     densityP=2.e7         ! plasmaspheric density (m^-3) to define plasmapause;
!    needed to avoid applying chorus and hiss at the same time

     ckeV_log(:)=log10(ckeV(:))
     ukeV_log(:)=log10(ukeV(:))
     hkeV_log(:)=log10(hkeV(:))
     

! determine edmax and edmin
  if (iUBC.eq.1) edmax=ukeV(iwu)
  if (ihiss.ge.1) edmax=hkeV(iwh)
  if (ichor.ge.1) edmax=ckeV(iwc)
  if (ihiss.ge.1) edmin=hkeV(1)
  if (ichor.ge.1) edmin=ckeV(1)
  if (iUBC.eq.1) edmin=ukeV(1)

   n=nspec  ! use last index to access electrons; different from CIMI

!   do j=1,ip
    do j=MinLonPar,MaxLonPar
     do i=1,iba(j)
        ro1=ro(i,j)
        ompe1=ompe(i,j)                           ! fpe/fce
        Cfactor=(BLc0/bo(i,j))*CHpower(i,j)/Cpower0
        Ufactor=(BLu0/bo(i,j))*UBCpower(i,j)/Upower0
        Hfactor=(BLh0/bo(i,j))*HIpower(i,j)/Hpower0

!        ompe1=cOmpe(1)  ! for testing; this is because we do not have realistic plasmasphere now

        if (ompe1.ge.cOmpe(1).and.ompe1.le.cOmpe(ipc).and.ro1.ge.r_wave) then
            ! Set up the energy grid, ein
            emin=minval(ekev(n,i,j,1,1:ik))
            emax=maxval(ekev(n,i,j,iw,1:ik))
            ein(1)=emin
            ein(ie)=emax
            ein_log(1)=log10(emin)
            ein_log(ie)=log10(emax)
            x0=ein_log(1)
            x2=ein_log(ie)
            x=(x2-x0)/(ie-1)
            do k=2,ie-1
               ein_log(k)=x0+(k-1)*x
               ein(k)=10.**ein_log(k)
            enddo

            ! Map psd to ein grid, f2d
            f2d(:,:)=0.
            do m=1,ik
               do k=1,iw
                f1d(k)=-50.                      ! f1d is log(psd)
                if(f2(n,i,j,k,m).gt.0.)f1d(k)=log10(f2(n,i,j,k,m)/xjac(n,i,k))
                
                if (k.gt.iw2(n,m)) f1d(k)=f1d(iw2(n,m))
                ekevlog(k,m)=log10(ekev(n,i,j,k,m))
                e1d(k)=ekevlog(k,m)              ! e1d is log(ekev)
               enddo

               do k=1,ie
                  einlog=ein_log(k)
                  if (einlog.lt.e1d(1)) einlog=e1d(1)
                  if (einlog.le.e1d(iw)) then
                     call lintpIM(e1d,f1d,iw,einlog,x)
                     f2d(k,m)=10.**x      ! f2d is psd
                  endif
               enddo
            enddo


             ! calcuate ao, dao, and Gjac
            do m=0,ik+1
               ao(m)=asin(y(i,j,m))
               Gjac(m)=tya(i,j,m)*y(i,j,m)*sqrt(1.-y(i,j,m)*y(i,j,m))
            enddo
            do m=1,ik
               dao(m)=0.5*(ao(m+1)-ao(m-1))
            enddo

            df=0.
         
            do k=1,ie      ! *****
             if (ein(k).ge.edmin.and.ein(k).le.edmax) then
               ! calculate DD, Daoao*Gjac
               do m=0,ik+1
                 ao_d=ao(m)*180./pi
                 if (ao_d.lt.cPA(1)) ao_d=cPA(1)
                 if (ao_d.gt.cPA(ipa)) ao_d=cPA(ipa)
                 cDaoao=0.
                 if (ichor.ge.1.and.density(i,j).le.densityP) then
                    if (ein(k).ge.ckeV(1).and.ein(k).le.ckeV(iwc)) &
                    call lintp3IM(cOmpe,ckeV_log,cPA,cDaa,ipc,iwc, &
                                 ipa,ompe1,ein_log(k),ao_d,cDaoao)
                 endif
                 uDaoao=0.
                 if (iUBC.eq.1.and.density(i,j).le.densityP) then
                    if (ein(k).ge.ukeV(1).and.ein(k).le.ukeV(iwu)) &
                     call lintp3IM(uOmpe,ukeV_log,cPA,uDaa,ipu,iwu, &
                                 ipa,ompe1,ein_log(k),ao_d,uDaoao)
                 endif
                 hDaoao=0.
                 if(ihiss.ge.1.and.ompe1.ge.2..and.density(i,j).gt.densityP)then
                     if (ein(k).ge.hkeV(1).and.ein(k).le.hkeV(iwh))&
                     call lintp3IM(hOmpe,hkeV_log,cPA,hDaa,iph,iwh, &
                                 ipa,ompe1,ein_log(k),ao_d,hDaoao)
                 endif
                 Daoao=cDaoao*Cfactor+uDaoao*Ufactor+hDaoao*Hfactor
                 ! Next 5 are subject to remove: BEGIN TESTING:
                 !if ( ro(i,j).ge.3.0 .and. ro(i,j).le.5.0 ) then 
                 !  Daoao=0.01   ! artificially force diffusion time to be 100sec
                 !else
                 ! Daoao=0.
                 !endif
                 ! END TESTING
 
                 DD(m)=Daoao*Gjac(m)
               enddo

               ! calculate up and um
               u_mx=0.
               do m=1,ik
                  factor_1=dao(m)*(ao(m)-ao(m-1))
                  factor1=dao(m)*(ao(m+1)-ao(m))
                  DDm=0.5*(DD(m)+DD(m-1))
                  DDp=0.5*(DD(m)+DD(m+1))
                  um(m)=dt*DDm/factor_1/Gjac(m)
                  up(m)=dt*DDp/factor1/Gjac(m)
                  if (factor_1.eq.0 .or. factor1.eq.0 .or. Gjac(m).eq.0) then
                    write(*,*) 'WARNING: CIMI: Diffuseaa: Null in denominator!'
                  ! if factor1 or factor_1 eq 0 check fied tracing and Bo and sinA grid
                  endif
                  ump_mx=max(abs(up(m)),abs(um(m)))
                  if (ump_mx.gt.u_mx) u_mx=ump_mx
               enddo

               ! reduce time step if u_mx > u_max
               irun=ifix(u_mx/u_max)+1
               do m=1,ik
                  um(m)=um(m)/irun
                  up(m)=up(m)/irun
               enddo
               u_mx=u_mx/irun     ! new u_mx < u_max

               ! determine the implicitness, xlam
               if (u_mx.le.1.) then
                  xlam=0.5              ! Crank-Nicolson
               else
                  u_mx_log=log10(u_mx)
                  xlam=0.5*u_mx_log/u_max_log+0.5
               endif
               alam=1.-xlam

               ! calculate a1d, b1d, c1d
               do m=1,ik
                  a1d(m)=-xlam*um(m)
                  b1d(m)=1.+xlam*(um(m)+up(m))
                  c1d(m)=-xlam*up(m)
               enddo

               ! start diffusion in ao
               f0(1:ik)=f2d(k,1:ik)
               do nn=1,irun
                  f0(0)=f0(1)
                  f0(ik+1)=f0(ik)
                  fr(1)=alam*um(1)*f0(0)+(1.-alam*(up(1)+um(1)))*f0(1)+ &
                        alam*up(1)*f0(2)+xlam*um(1)*f0(0)
                  do m=2,ik-1    ! calculate the RHS of matrix equation
                     fr(m)=alam*um(m)*f0(m-1)+(1.-alam*(up(m)+um(m)))*f0(m)+ &
                           alam*up(m)*f0(m+1)
                  enddo
                  fr(ik)=alam*um(ik)*f0(ik-1)+(1.-alam*(up(ik)+um(ik)))*f0(ik) &
                         +alam*up(ik)*f0(ik+1)+xlam*up(ik)*f0(ik+1)
                  call tridagIM(a1d,b1d,c1d,fr,f0(1:ik),ik,ier)
                  if (ier.eq.1) call CON_STOP('ERROR in CIMI: Diffusionaa: tridag failed,ier=1')
                  if (ier.eq.2) call CON_STOP('ERROR in CIMI: Diffusionaa: tridag failed,ier=2')
               enddo

               do m=1,ik
                  df(k,m)=f0(m)-f2d(k,m)    ! df is differential psd
               enddo

             endif
            enddo         ! end of do k=1,ie      ! *****

            ! map psd back to M grid
            do m=1,ik
               do k=1,ie
                  df1(k)=df(k,m)
               enddo
               do k=1,iw2(n,m)
                  call lintpIM(ein_log,df1,ie,ekevlog(k,m),dpsd)
                  f2(n,i,j,k,m)=f2(n,i,j,k,m)+xjac(n,i,k)*dpsd
                  if (f2(n,i,j,k,m).lt.0.) f2(n,i,j,k,m)=0.
               enddo
            enddo
        endif    ! end of if (ompe1.ge.cOmpe(1).and.ompe1.le.cOmpe(ipc).and...

     enddo
   enddo         ! end of j loop

end subroutine diffuse_aa

!****************************************************************************
!                            diffuse_EE
!  Routine calculates the change of electron distributions due to
!  diffusion in E.
!****************************************************************************
  subroutine diffuse_EE(f2,dt,xmm,xjac,iw2,iba)

  use DensityTemp
  use ModMPI
  use ModCrcmGrid,   ONLY: MinLonPar,MaxLonPar

!  use constants
!  use cimigrid_dim CHECK je and ig
  implicit none
  integer i,j,k,m,k1,k2,iww,irun,nrun,ier,n
  real xjac(nspec,ir,iw),xmm(nspec,0:iw+1)  ! xjac and xmm have different dimension from CIMI !
  integer iba(ip),iw2(nspec,ik)
  real dt,Eo,f2(nspec,ir,ip,iw,ik),fl(iw), &
       xlam,alam,um(iw),up(iw),ompe1,ro1,PA1, &
       fr(iw),a1d(iw),b1d(iw),c1d(iw),&
       factor_1, &
       Hfactor,Upower0, &
       factor1,gjac_1,gjac1,gjac(iw),E1(0:iw),Ep,DEE(0:iw+1), &
       DDm,DDp,DEEc,DEEu,DEEh,u_mx,ump_mx,Enor(iw),dEn(iw), & ! CHECKED
       WM1,Wo,Cfactor,Ufactor,edmin,edmax, & ! CHECKED
       ckeV_log(iwc),ukeV_log(iwu),hkeV_log(iwh),Ep_log,densityP ! CHECKED

  Eo=511.               ! electron rest energy in keV
  xlam=0.5              ! implicitness in solving diffusion equation
  alam=1.-xlam
  Upower0=10000.        ! coeff based on UB chorus power of (100pT)^2
  densityP=2.e7         ! plasmaspheric density (m^-3) to define plasmapause

  ckeV_log(:)=log10(ckeV(:))
  ukeV_log(:)=log10(ukeV(:))
  hkeV_log(:)=log10(hkeV(:))

! determine edmax and edmin
  if (iUBC.eq.1) edmax=ukeV(iwu)
  if (ihiss.ge.1) edmax=hkeV(iwh)
  if (ichor.ge.1) edmax=ckeV(iwc)
  if (ihiss.ge.1) edmin=hkeV(1)
  if (ichor.ge.1) edmin=ckeV(1)
  if (iUBC.eq.1) edmin=ukeV(1)

 n=nspec  ! use last index to access electrons; different from CIMI

!   do j=1,ip
    do j=MinLonPar,MaxLonPar
     do i=1,iba(j)
        ro1=ro(i,j)
        Cfactor=(BLc0/bo(i,j))*CHpower(i,j)/Cpower0
        Ufactor=(BLu0/bo(i,j))*UBCpower(i,j)/Upower0
        Hfactor=(BLh0/bo(i,j))*HIpower(i,j)/Hpower0
        ompe1=ompe(i,j)                         ! fpe/fce

        if (testDiff_EE) ompe1=cOmpe(1)

        if (ompe1.ge.cOmpe(1).and.ompe1.le.cOmpe(ipc).and.ro1.ge.r_wave) then
           do m=1,ik
              PA1=asin(y(i,j,m))*180./pi        ! pitch angle in degree
              if (PA1.lt.cPA(1)) PA1=cPA(1)
              if (PA1.gt.cPA(ipa)) PA1=cPA(ipa)
              do k=1,iw
                 Enor(k)=ekev(n,i,j,k,m)/Eo
                 gjac(k)=(Enor(k)+1.)*sqrt(Enor(k)*(Enor(k)+2.))
              enddo

              ! determine k1 and k2, corresponding to edmin and edmax
              k1=iw2(n,m)
              findk1: do k=2,iw2(n,m)
                 if (edmin.le.ekev(n,i,j,k,m)) then
                    k1=k
                    if (k1.eq.1) k1=2
                    exit findk1
                 endif
              enddo findk1
              k2=1
              findk2: do k=iw2(n,m)-1,1,-1
                 if (edmax.ge.ekev(n,i,j,k,m)) then
                    k2=k
                    exit findk2
                 endif
              enddo findk2
              iww=k2-k1+1

!    if (i.eq.30 .and. j.eq.5) write(*,*) 'm',m,'energy range in D_ee',ekev(n,i,j,k1,m),ekev(n,i,j,k2,m),k1,k2

              ! find DEE/E^2 and dEn

             ! find DEE/E^2 and dEn
              Wo=Eo*1.6e-16/bm(i,j,m)   ! normalization factor of mag moment
              do k=k1-1,k2
                 WM1=sqrt(xmm(n,k)*xmm(n,k+1))
                 E1(k)=sqrt(2.*WM1/Wo+1.)-1.    ! normalized kinetic energy
                 Ep=E1(k)*Eo                    ! Ep in keV
                 Ep_log=log10(Ep)
                 if (k.ge.k1) dEn(k)=E1(k)-E1(k-1)
                 DEEc=0.
                 if (ichor.ge.1.and.density(i,j).le.densityP) then
                    if (Ep.ge.ckeV(1).and.Ep.le.ckeV(iwc)) call lintp3IM(cOmpe, &
                        ckeV_log,cPA,cDEE,ipc,iwc,ipa,ompe1,Ep_log,PA1,DEEc)
                 endif
                 DEEu=0.
                 if (iUBC.eq.1.and.density(i,j).le.densityP) then
                    if (Ep.ge.ukeV(1).and.Ep.le.ukeV(iwu)) call lintp3IM(uOmpe, &
                        ukeV_log,cPA,uDEE,ipu,iwu,ipa,ompe1,Ep_log,PA1,DEEu)
                 endif
                 DEEh=0.
                 if(ihiss.ge.1.and.ompe1.ge.2..and.density(i,j).gt.densityP)then
                    if (Ep.ge.hkeV(1).and.Ep.le.hkeV(iwh)) call lintp3IM(hOmpe, &
                        hkeV_log,cPA,hDEE,iph,iwh,ipa,ompe1,Ep_log,PA1,DEEh)
                 endif
                 DEE(k)=DEEc*Cfactor+DEEu*Ufactor+DEEh*Hfactor ! this is DEE/E^2
                   
                 if (testDiff_EE) then 
                  if ( ro(i,j).ge.3.0 .and. ro(i,j).le.5.0 ) then
                   DEE(k)=0.01   ! artificially force diffusion time to be 1/DEE sec;
                  ! over 1/DEE the psd should evolve at energy range E0~0.5MeV;
                  ! this is how the current setup is normalized. 
                  else
                   DEE(k)=0.
                  endif
                 endif

              enddo

              ! find the Courant number
              u_mx=0.
              do k=k1,k2
                 factor_1=dEn(k)*(Enor(k)-Enor(k-1))   ! normalized factor
                 factor1=dEn(k)*(Enor(k+1)-Enor(k))    !
                 gjac_1=(E1(k-1)+1.)*sqrt(E1(k-1)*(E1(k-1)+2.))
                 gjac1=(E1(k)+1.)*sqrt(E1(k)*(E1(k)+2.))
                 DDm=DEE(k-1)*E1(k-1)*E1(k-1)*gjac_1
                 DDp=DEE(k)*E1(k)*E1(k)*gjac1
                 if ( (factor_1.eq.0) .or. (factor1.eq.0) .or. (gjac_1.eq.0) &
                 .or. (gjac1.eq.0))  write(*,*) 'WARNING in CIMI: diffuse_EE:um,up:Null in denominator!!!'
                 um(k)=dt*DDm/factor_1/gjac(k)
                 up(k)=dt*DDp/factor1/gjac(k)
                 ump_mx=max(abs(up(k)),abs(um(k)))
                 if (ump_mx.gt.u_mx) u_mx=ump_mx
              enddo

              ! reduce time step size if up or um is too large
              irun=ifix(u_mx)+1
              do k=k1,k2
                 um(k)=um(k)/irun
                 up(k)=up(k)/irun
                 a1d(k)=-xlam*um(k)
                 b1d(k)=1.+xlam*(um(k)+up(k))
                 c1d(k)=-xlam*up(k)
              enddo
              ! Start diffusion in E
              do k=k1-1,k2+1
                 fl(k)=f2(n,i,j,k,m)/xjac(n,i,k)   ! fl is psd
               if (xjac(n,i,k).eq.0) write(*,*) 'WARNING in CIMI: diffuse_EE: xjac:Null in denominator!!!'  
              enddo
              ! force flat psd across the boundaries; this is done to force zero
              ! flux across the boundaries and close boundaries
                fl(k1-1)=fl(k1+1)
                fl(k1)=fl(k1+1)
                fl(k2+1)=fl(k2-1)
                fl(k2)=fl(k2-1)       
       
              do nrun=1,irun
                 fr(k1)=alam*um(k1)*fl(k1-1)+(1.-alam*(up(k1)+um(k1)))*&
                        fl(k1)+alam*up(k1)*fl(k1+1)+xlam*um(k1)*fl(k1-1)
                 do k=k1+1,k2-1    ! calculate the RHS of matrix equation
                    fr(k)=alam*um(k)*fl(k-1)+(1.-alam*(up(k)+um(k)))*&
                          fl(k)+alam*up(k)*fl(k+1)
                 enddo
                 fr(k2)=alam*um(k2)*fl(k2-1)+(1.-alam*(up(k2)+um(k2)))*&
                        fl(k2)+alam*up(k2)*fl(k2+1)+xlam*up(k2)*fl(k2+1)
               !  if (i.eq.30 .and. j.eq.5) write(*,*) 'nrun',nrun,'fr1',fr(k1),fr(k1+1),'fr2',fr(k2-1),fr(k2)           
               !  if (i.eq.30 .and. j.eq.5) write(*,*) 'nrun',nrun,'fl1',fl(k1),fl(k1+1),'fl2',fl(k2-1),fl(k2)
               !  if (i.eq.30 .and. j.eq.5) write(*,*) 'nrun',nrun,'f2_right',f2(n,i,j,k1,m),&
               !  f2(n,i,j,k1+1,m),'f2_left',f2(n,i,j,k2-1,m),f2(n,i,j,k2,m)

                 call tridagIM(a1d(k1),b1d(k1),c1d(k1),fr(k1),fl(k1),iww,ier)
              if (ier.eq.1) call CON_STOP('ERROR in CIMI: diffuse_EE:tridag failed,ier=1')
              if (ier.eq.2) call CON_STOP('ERROR in CIMI: diffuse_EE:tridag failed,ier=2')
              enddo

              ! get back f2
              ! if (i.eq.30 .and. j.eq.5 .and. m.eq.5) write(*,*) 'psd, ', fl(k1:k2)
              do k=k1+1,k2-1 ! change made by NB to avoid two 'ghost cells'. Originally was do k=k1,k2
              ! do k=k1,k2
                 if (fl(k).lt.0.) fl(k)=0.
                 f2(n,i,j,k,m)=fl(k)*xjac(n,i,k)
              enddo

           enddo       ! end of m=1,ik
        endif          ! end of if (ompe1.ge.cOmpe(1).and.ompe1.le.cOmpe(ipc)..)

     enddo       ! end of i loop
   enddo         ! end of j loop

end subroutine diffuse_EE

 EndModule
