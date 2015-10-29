Module DensityTemp
use ModCrcmGrid,ONLY: ir=>np, ip=>nt
implicit none
real,parameter ::  density(ir,ip)=10.
EndModule


!!!!!!!!!!!!!!!!  Main wave module !!!!!!!!!!!!!1
Module ModWaveDiff

use ModCrcmGrid,ONLY: ir=>np, ip=>nt, iw=>nm , ik=>nk, ie=>neng
use ModCrcmPlanet,ONLY: nspec
use ModFieldTrace, ONLY: ekev,tya,y=>sinA,Bo,ro,xmlto,lintpIM
use ModNumConst,       ONLY: pi => cPi


implicit none

logical            :: UseWaveDiffusion = .false.
logical            :: UseHiss  = .false.
logical            :: UseChorus  = .false.
logical            :: UseChorusUB = .false.
real               :: DiffStartT = 0.
character(len=100)  :: HissWavesD = 'Hiss_UCLA_v1'
character(len=12)  :: ChorusWavesD= 'LBchorus_mer'
character(len=12)  :: ChorusUpperBandD = 'UBchorus_mer'


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

 integer j,k,m,ichor,iUBC,ihiss
 real Eo,c_cgs,pa,daa,dap,dpp,daE,dEE,Dnorm,cMeV,Lc0,Lh0
 real E_cgs,cE2,EEo,E2Eo,daan,dapn,dppn,daa0,dap0,dpp0
 real,parameter :: Lu0=6.
 character header*80
 real rc

! list of available wave diffusion coefficients models:
! name in PARAM file       file name 
! LBchorus__QZ             D_LBchorus_QZ.dat
! LBchorus_mer             D_LBchorus_merge.dat
! Hiss_USLA_v1             D_hiss_UCLA.dat
! Hiss__Albert             D_hiss_Albert.dat

 ihiss=0
 if (UseHiss) then
   if ( HissWavesD .eq.'Hiss_UCLA_v1') ihiss=2
 else 
   ihiss=1
 endif

 ichor=0
 if (UseChorus) then 
   if (ChorusWavesD.eq.'LBchorus_mer') ichor=2
 else 
  ichor=1
 endif

 iUBC=0
 if (UseChorusUB) iUBC=1

 write(*,*) '*** Wave model is on ***'
 write(*,*) 'IfUseHiss:', UseHiss
 if (UseHiss) write(*,*) 'Dcoef for Hiss model:', HissWavesD

 write(*,*) 'IfUseChorus Lower Band:', UseChorus
 if (UseChorus)  write(*,*) 'Dcoef for Chorus model:', ChorusWavesD

 write(*,*) 'IfUseChorus Upper Band:', UseChorusUB
 if (UseChorusUB)  write(*,*) 'Dcoef for Upper Band Chorus:', ChorusUpperBandD

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

! Read LB chorus ckeV(0.1keV-10MeV), Daa, DEE, DaE at L = 6.5
  write(*,*) 'Reading chorus files:',ichor
  if (ichor.eq.1) open(unit=UnitTmp_,file='IM/D_LBchorus_QZ.dat',status='old')
  if (ichor.eq.2) open(unit=UnitTmp_,file='IM/D_LBchorus_merge.dat',status='old')
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

! Read UB chorus Daa, DEE, DaE at L = 6.0
  BLu0=xme/(Lu0*re_m)**3       ! dipole B at Lu0
  uDaa(:,:,:)=0.
  uDaE(:,:,:)=0.
  uDEE(:,:,:)=0.
  if (iUBC.eq.1) then
     write(*,*) 'reading UBchorus'
      open(unit=UnitTmp_,file='IM/D_UBchorus.dat',status='old')
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

! Read hiss Daa, DEE, DaE at L = 5.5

  write(*,*) 'reading hiss waves',ihiss
  if (ihiss.eq.2) open(unit=UnitTmp_,file='IM/D_hiss_UCLA.dat',status='old')
  if (ihiss.eq.1) open(unit=UnitTmp_,file='IM/D_hiss_Albert.dat',status='old')
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
           stop
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


!!!!!!!!!!!!!!!!!!!   main diffusion aa subroutine  !!!!!!!!!!!!

   subroutine diffuse_aa(f2,dt,xjac,iba,iw2)
      use DensityTemp
      use ModMPI
      use ModCrcmGrid,   ONLY: MinLonPar,MaxLonPar
      implicit none
      integer,parameter :: ie=40
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
        ompe1=ompe(i,j)                           ! fpe/fce
        Cfactor=(BLc0/bo(i,j))*CHpower(i,j)/Cpower0
        Ufactor=(BLu0/bo(i,j))*UBCpower(i,j)/Upower0
        Hfactor=(BLh0/bo(i,j))*HIpower(i,j)/Hpower0

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



   subroutine WaveDipoleTya
! calculates Tya=1/R_0 ( \int_0^sm ds/cos(a_0) ) in the dipole tilt; only for
! wave test

   end subroutine WaveDipoleTya

   subroutine WaveTest_fpe_fce
! calculates ratio fpe/fce with dipole tilt and exponential plasmasphere

   end subroutine WaveTest_fpe_fce 


Endmodule

