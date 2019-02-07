!-------------------------------------------------------------------------------
!                            plasmasphere_cimi.f90
! 
! A code to calculate plaamasphere flux tube content. The code is coupled with
! the cimi code, which provides the grid and the convection potential.
!
! Created on 1 February 2019 by Mei-Ching Fok. Code 673, NASA GSFC.
!-------------------------------------------------------------------------------

Module ModPlasmaphere
  implicit none
  private

  !constants for the plasmaphere module
  real,parameter :: pi=3.14159265358979
  real,parameter :: re_m=6.375e6        ! earth's radius (m)
  real,parameter :: DipM=8.e15          ! Dipole moment in T m^3
  real,parameter :: rc=1.0157           ! radial distance of ionosphere in RE
 

  !grid and saturation density values
  integer,parameter :: nl=209,np=192
  integer ibp(np)
  real xlatp(nl),pphi(np),rop(nl,np),phip(nl,np),dphip(nl,np),volumep(nl,np),&
       Nsat(nl,np),Nion(nl,np),potentp(nl,np),denSat(nl,np)

  !plasmasphere density on plasmasphere grid
  real :: densityp(nl,np)

  !plasmasphere density on CIMI grid
  real,allocatable,public :: PlasDensity_C(:,:)

  !copy of CIMI grid for use in interpolation
  integer :: nLatCimi,nLonCimi
  real, allocatable :: xlatCimi(:),phiCimi(:),roCimi(:,:),volumeCimi(:,:),&
       potentCimi(:,:),ibCIMI(:)

  
  ! Setup the trough by trough(6.6)*dipolevolumep(6.6) as in pbo_2.f
  ! number of particle per unit magnetic flux at trough 
  real, parameter :: NTro=9.40e20   
  
  real Bi(nl),dilatp(nl),dlatp(nl)

  character*1 Ndon(nl,np),Sdon(nl,np)
contains

  !============================================================================
  subroutine unit_test_plasmasphere 
    integer i,j,n,Nstep,Nstep4
    real xlat1,xlat2,dlat,xlatp1,dphi,volume0,volume1,coslat,Lshell(nl),L7, &
         tmin,tmax,tsec,Kpi,Kpf,xKp,dKpdt,PCP,&
         potent1,delt
    
    ! Setup grid (or provided in cimi)
    xlat1=20.
    xlat2=71.
    xlat1=xlat1*pi/180.
    xlat2=xlat2*pi/180.
    dlat=(1./cos(xlat2)-1./cos(xlat1))/(nl-1)
    do i=1,nl
       xlatp1=1./cos(xlat1)+(i-1)*dlat
       xlatp(i)=acos(1./xlatp1) 
    enddo
    dphi=2.*pi/np
    do j=1,np
       pphi(j)=(j-1)*dphi
    enddo
    
    ! Setup ibp, rop, phip and volume assuming dipole (or provided in cimi)
    ibp(:)=nl
    volume0=9.14e-4*re_m**4/DipM
    do i=1,nl
       coslat=cos(xlatp(i))
       Lshell(i)=rc/coslat/coslat
       L7=Lshell(i)**7  
       volume1=volume0*L7
       do j=1,np
          rop(i,j)=Lshell(i)
          phip(i,j)=pphi(j)
          volumep(i,j)=volume1
       enddo
    enddo
    
    
    ! More setup and output the initial density
    delt=30.      ! call plasmasphere every delt seconds
    tmin=0.
    tmax=86400.   
    Nstep=ifix((tmax-tmin)/delt)
    Nstep4=Nstep/4
    Kpi=1.
    Kpf=8.
    dKpdt=2.*(Kpf-Kpi)/(tmax-tmin)       ! assume Kp changing linearly
    xKp=Kpi
    ! xKp=3.   
    tsec=tmin
    open(unit=1,file='plasmasphere_cimi.txt')
    write(1,*) nl,np,'     ; nl,np'
    write(1,*) rop
    write(1,*) phip
    write(1,*) tsec,xKp,'     ; tsec,Kp'
    write(1,*) denSat
    write(*,*) 'tsec,Kp ',tsec,xKp
    
    ! Run the model with Vollend-Stern field (or potential in cimi)
    do n=1,Nstep
       tsec=tmin+n*delt
       if (n.le.Nstep4) xKp=xKp+dKpdt*delt
       if (n.gt.Nstep4) xKp=xKp-dKpdt*delt
       if (xKp.lt.Kpi) xKp=Kpi
       !  xKp=3.
       PCP=xkp*1.4e4               ! assume PCP increasing linearly with Kp
       do i=1,nl
          potent1=PCP*(Lshell(i)/Lshell(nl))**2     ! VS with shielding factor 2
          do j=1,np
             potentp(i,j)=potent1*sin(pphi(j))
          enddo
       enddo
       call plasmasphere(delt)
       ! output densityp
       if (mod(tsec,1800.).eq.0) then
          write(*,*) 'tsec,Kp ',tsec,xKp
          write(1,*) tsec,xKp,'     ; tsec,Kp'
          do i=1,nl
             do j=1,np
                densityp(i,j)=Nion(i,j)/volumep(i,j)
             enddo
          enddo
          write(1,*) densityp
       endif
    enddo
    close(1)
    
  end subroutine unit_test_plasmasphere
  
!-------------------------------------------------------------------------------
  subroutine plasmasphere(delt) 
!-------------------------------------------------------------------------------
! A plasmasphere model calculates number of plasmasphere ions per unit magnetic
! flux, Nion. This code is similar to Dan Ober's model (pbo_2.f).
!
! INPUT
! potentp(nl,np): potential in Volt at (xlatp,pphi)
! delt: simulation time in second of this call
!
! INPUT/OUTPUT
! Nion: number of ions per unit magnetic flux.    
!-------------------------------------------------------------------------------
    
    integer n,nstep,nrun

    real dtmax,dt,vl(nl,np),vp(nl,np)
    real dt1,cl(nl,np),cp(nl,np)    
    
    ! determine time step 
    dtmax=30.                    ! maximum time step in second
    nstep=ceiling(delt/dtmax)
    dt=delt/nstep                ! time step in second
    
    ! calculate convection velocity and Courant numbers
    call Vconvect(dt, dt1,nrun,vl,vp,cl,cp)
    
    ! time loop to update Nion    
    do n=1,nstep    
       call trough(nl,np)
       call drift_pl(nrun,vl,cl,cp,dt1)
       call refill_loss(dt)
    enddo
    
  end subroutine plasmasphere
  !==========================================================================
  subroutine init_plasmasphere_grid(nLatIn,nLonIn,xlatIn,phiIn,ibIn)
    integer, intent(in) :: nLatIn,nLonIn
    real,    intent(in) :: xlatIn(nLatIn),phiIn(nLonIn)
    integer, intent(in) :: ibIn(nLonIn)

    real xlatl,xlatu 
    real dphi,delt
    real DipMr3,coslat2,coslat,phi1,xsm
    integer i,j
    !--------------------------------------------------------------------------
    !save cimi grid dimension
    nLatCimi=LatIn
    nLonCimi=nLonIn

    !allocate cimi grid
    if (allocated(xlatCimi)) deallocate(xlatCimi)
    allocate(xlatCimi(nLatCimi))
    if (allocated(phiCimi)) deallocate(phiCimi)
    allocate(phiCimi(nLonCimi))
    if (allocated(ibCimi)) deallocate(ibCimi)
    allocate(ibCimi(nLonCimi))

    !allocate cimi arrays to hold input and output variables
    if (allocated(roCimi)) deallocate(roCimi)
    allocate(roCimi(nLatCimi,nLonCimi))
    
    if (allocated(volumeCimi)) deallocate(volumeCimi)
    allocate(volumeCimi(nLatCimi,nLonCimi))

    if (allocated(roCimi)) deallocate(roCimi)
    allocate(roCimi(nLatCimi,nLonCimi))

    if (allocated(potentCimi)) deallocate(potentCimi)
    allocate(potentCimi(nLatCimi,nLonCimi))

    if (allocated(PlasDensity_C)) deallocate(PlasDensity_C)
    allocate(PlasDensity_C(nLatCimi,nLonCimi))

    
    !save cimi grid
    xlatCimi = xlatIn
    phiCimi  = phiIn
    ibCimi   = ibIn
    
    !set up plasmasphere grid
    xlat1=20.
    xlat2=71.
    xlat1=xlat1*pi/180.
    xlat2=xlat2*pi/180.
    dlat=(1./cos(xlat2)-1./cos(xlat1))/(nl-1)
    do i=1,nl
       xlatp1=1./cos(xlat1)+(i-1)*dlat
       xlatp(i)=acos(1./xlatp1) 
    enddo
    dphi=2.*pi/np
    do j=1,np
       pphi(j)=(j-1)*dphi
    enddo
    xlat1=20.
    xlat2=71.
    xlat1=xlat1*pi/180.
    xlat2=xlat2*pi/180.
    dlat=(1./cos(xlat2)-1./cos(xlat1))/(nl-1)
    do i=1,nl
       xlatp1=1./cos(xlat1)+(i-1)*dlat
       xlatp(i)=acos(1./xlatp1) 
    enddo
    dphip=2.*pi/np
    do j=1,np
       phip(j)=(j-1)*dphip
    enddo

    ! Setup grid intervals
    xlatl=1.5*xlatp(1)-0.5*xlatp(2)       ! lower boundary of xlatp(1)
    xlatu=1.5*xlatp(nl)-0.5*xlatp(nl-1)   ! upper boundary of xlatp(1)
    do i=1,nl-1
       dilatp(i)=xlatp(i+1)-xlatp(i)
       if (i.gt.1) dlatp(i)=0.5*(xlatp(i+1)-xlatp(i-1))
    enddo
    dilatp(nl)=2.*(xlatu-xlatp(nl))      
    dlatp(1)=xlatp(1)-xlatl+0.5*dilatp(1) 
    dlatp(nl)=xlatu-xlatp(nl)+0.5*dilatp(nl-1)
    
    ! find Bi, Ndon, Sdon
    DipMr3=DipM/(rc*re_m)**3
    do j=1,np
       do i=1,ibp(j)
          ! find Ndon
          coslat=cos(xlatp(i))
          coslat2=coslat*coslat
          Bi(i)=DipMr3*sqrt(4.-3.*coslat2)         ! B in T at northern ionosphere
          phi1=pphi(j)+pi
          xsm=rc*coslat*cos(phi1)
          if (xsm.gt.0.) Ndon(i,j)='d'
          if (xsm.le.0.) Ndon(i,j)='n'
          Sdon(i,j)=Ndon(i,j)
       enddo
    enddo


    !set initial plasmasphere condition
    if (IsRestart) then
       ! Setup initial flux tube content assuming saturated condition
       Nsat(:,:)=0.
       do j=1,np
          do i=1,ibp(j)
             denSat(i,j)=1.e6*10.0**(-0.3145*rop(i,j)+3.9043)   ! saturated density
             Nsat(i,j)=denSat(i,j)*volumep(i,j)
          enddo
       enddo
       Nion(:,:)=Nsat(:,:)
    else
       call load_restart
    endif
    
  end subroutine init_plasmasphere_grid

  !=========================================================================
  subroutine load_restart
    integer :: iLat,iLon
    do iLon=1,np
       do iLat=1,nl
  end subroutine load_restart
  !==========================================================================
  subroutine cimi_put_to_plasmasphere(nLatIn,nLonIn,roIn,volumeIn,potentIn)
    integer, intent(in) :: nLatIn, nLonIn
    real,    intent(in) :: roIn(nLatIn,nLonIn),volumeIn(nLatIn,nLonIn),&
         potentIn(nLatIn,nLonIn)
    
    !-------------------------------------------------------------------------
    roCimi     = roIn
    volumeCimi = volumeIn
    potentCimi = potentIn

  end subroutine cimi_put_to_plasmasphere
  
  !==========================================================================
  subroutine interpolate_cimi_to_plasmasphere
    use ModInterpolate, ONLY: bilinear
    use ModNumConst,    ONLY: cDegToRad,cPi,cTwoPi
    integer :: iLat, iLon
    real :: LatLon_D(2)
    !--------------------------------------------------------------------------
    do iLon = 1, np
       do iLat = 1, nl
          LatLon_D(1) = xlatp(iLat)
          LatLon_D(2) = phip(iLon)!modulo(phip(iLon)+cPi,cTwoPi)

          rop(iLat,iLon) = &
               bilinear(roCimi,1,nLatCimi,1,nLonCimi,LatLon_D, &
               xlatCimi,phiCimi,DoExtrapolate=.true.)
          volumep(iLat,iLon) = &
               bilinear(volumeCimi,1,nLatCimi,1,nLonCimi,LatLon_D, &
               xlatCimi,phiCimi,DoExtrapolate=.true.)
          potentp(iLat,iLon) = &
               bilinear(potentCimi,1,nLatCimi,1,nLonCimi,LatLon_D, &
               xlatCimi,phiCimi,DoExtrapolate=.true.)
       enddo
    enddo

    !interpolate ibCimi to ibp
    !first find find latitude boundary array for CIMI
    if (.not.allocated(LatTmp))allocate(LatTmp(nLonCimi))
    do iLon = 1, nLonCimi
       LatTmp(iLon)=xlatCimi(ibCimi(iLon))
       
    enddo

    !for each longitude, interpolate to find boundary lat
    do iLon=1,np
       phi(iLon)
       
    end do
  end subroutine interpolate_cimi_to_plasmasphere
     
  
  !==========================================================================
  subroutine interpolate_plasmasphere_cimi
    use ModInterpolate, ONLY: bilinear
    use ModNumConst,    ONLY: cDegToRad,cPi,cTwoPi

    integer :: iLat, iLon
    real :: LatLon_D(2)
    !--------------------------------------------------------------------------
    do iLon = 1, nLonCimi
       do iLat = 1, nLatCimi
          LatLon_D(1) = xlatCimi(iLat)
          LatLon_D(2) = phiCimi(iLon)!modulo(phip(iLon)+cPi,cTwoPi)

          PlasDensity_C(iLat,iLon) = &
               bilinear(densityp,1,nl,1,np,LatLon_D, &
               xlatp,phip,DoExtrapolate=.true.)
       enddo
    enddo
              
  end subroutine interpolate_plasmasphere_cimi


!-------------------------------------------------------------------------------
  subroutine Vconvect(dt,dt1,nrun,vl,vp,cl,cp)
!-------------------------------------------------------------------------------
! Routine calculates the convection velocities vl and vp
! Input: nl,np,ibp,xlatp,dilatp,dphi,DipM,rc,re,dt
! Output: dt1,nrun,vl,vp,cl,cp

  implicit none

  integer i,j,j0,j2,i0,i2,nrun
  real vl(nl,np),vp(nl,np)
  real cl(nl,np),cp(nl,np),dlatp(nl)
  real pi,dphi2,kfactor,cor,ksai,xlatp1,ksai1,sf0,sf2,dlat2,dt,dt1,cmax,cmx

  pi=acos(-1.)
  dphi2=dphip*2.
  kfactor=DipM/rc/re_m
  cor=2.*pi/86400.                  ! corotation speed in rad/s

! Find vl, vp
  do i=1,nl
     ksai=kfactor*sin(2.*xlatp(i))
     xlatp1=xlatp(i)+0.5*dilatp(i)    
     ksai1=kfactor*sin(2.*xlatp1)         ! ksai at i+0.5
     do j=1,np
        j0=j-1
        if (j0.lt.1) j0=j0+np
        j2=j+1
        if (j2.gt.np) j2=j2-np

        ! calculate vl
        if (ibp(j0).gt.i.and.ibp(j2).gt.i) then
           sf0=0.5*(potentp(i,j0)+potentp(i+1,j0))
           sf2=0.5*(potentp(i,j2)+potentp(i+1,j2))
           vl(i,j)=-(sf2-sf0)/dphi2/ksai1        ! vl at (i+0.5,j)
        else
           vl(i,j)=vl(i-1,j)
        endif

        ! calculate vp
        if (ibp(j2).gt.i.and.ibp(j).gt.i) then
           i0=i-1
           if (i.eq.1) i0=1
           i2=i+1
           if (i.eq.nl) i2=nl
           dlat2=xlatp(i2)-xlatp(i0)
           sf0=0.5*(potentp(i0,j2)+potentp(i0,j))
           sf2=0.5*(potentp(i2,j2)+potentp(i2,j))
           vp(i,j)=cor+(sf2-sf0)/dlat2/ksai       ! vp@(i,j+0.5)
        else
           vp(i,j)=vp(i-1,j)
        endif
     enddo          ! end of j loop
  enddo             ! end of i loop
  
! Find nrun, new dt (dt1), cl and cp
  cmax=0.
  do i=1,nl
     do j=1,np
        cl(i,j)=dt/dlatp(i)*vl(i,j)
        cp(i,j)=dt/dphip*vp(i,j)
        cmx=max(abs(cl(i,j)),abs(cp(i,j))) 
        cmax=max(cmx,cmax) 
     enddo
  enddo
  nrun=ifix(cmax/0.25)+1     ! nrun to limit the Courant number
  dt1=dt/nrun
  if (nrun.gt.1) then       ! reduce cl and cp if nrun > 1
     cl(1:nl,1:np)=cl(1:nl,1:np)/nrun
     cp(1:nl,1:np)=cp(1:nl,1:np)/nrun
  endif

  end subroutine Vconvect


!-------------------------------------------------------------------------------
  subroutine drift_pl(nrun,vl,cl,cp,dt1)
!-------------------------------------------------------------------------------
! Routine updates Nion due to drift (convection+corotation)
    ! Input: nl,np,ibp,nrun,vl,cl,cp,dt1,dlatp,Nsat
    ! Input/Output: Nion

  implicit none

  integer i,j,j_1,n,nrun
  real vl(nl,np),dt1,f2(nl,np),fb1
  real cl(nl,np),cp(nl,np),fal(nl,np),fap(nl,np),sin2l(nl)

! calculate f2
  Nion(1,1:np)=Nsat(1,1:np)    ! saturation density for first cell
  do i=1,nl
     sin2l(i)=sin(2.*xlatp(i))
     f2(i,1:np)=sin2l(i)*Nion(i,1:np)
  enddo

! flux at outer boundary
  fb1=NTro*sin2l(nl)

! Update f2 by drift
  do n=1,nrun
     call inter_flux(nl,np,ibp,cl,cp,fb1,f2,fal,fap)     
     do j=1,np
        j_1=j-1
        if (j_1.lt.1) j_1=j_1+np
        do i=2,ibp(j)              
           f2(i,j)=f2(i,j)+dt1/dlatp(i)*(vl(i-1,j)*fal(i-1,j)-vl(i,j)*fal(i,j))&
                   +cp(i,j_1)*fap(i,j_1)-cp(i,j)*fap(i,j)  
           if (f2(i,j).lt.0.) then
              write(*,*) ' Error: f2(i,j).lt.0.'
              stop
           endif
        enddo
     enddo
  enddo 

! get tbe new Nion
  do j=1,np
     Nion(1:nl,j)=f2(1:nl,j)/sin2l(1:nl)
  enddo

  end subroutine drift_pl


!-------------------------------------------------------------------------------
  subroutine refill_loss(dt)
!-------------------------------------------------------------------------------
! Routine updates Nion due to refilling and loss.
! On the nightside, dN/dt = -N/(Bi*tau) and N = No*exp(-dt/tau)
! on the dayside, dN/dt = Fmax*(Nsat-N)/Nsat and 
!                 N = Nsat - (Nsat - No)*exp(-dt*Fmax/Nsat/Bi)
!
! Input: nl,np,ibp,dt,velume,Nsat,Bi,Ndon,Sdon
! Input/Output: Nion

  implicit none

  integer i,j
  real dt
  real tau,Fmax,factorN,factorD,tFNB

  Fmax=2.e12          ! limiting refilling flux in particles/m**2/sec
  tau=10.*86400.      ! nightside downward diffusion lifetime in second
  factorN=exp(-dt/tau)

  do j=1,np
     do i=1,ibp(j)
        ! at Northern ionosphere
        if (Ndon(i,j).eq.'d') then          ! dayside refilling
           tFNB=-dt*Fmax/Nsat(i,j)/Bi(i)
           factorD=exp(tFNB)
           Nion(i,j)=Nsat(i,j)-(Nsat(i,j)-Nion(i,j))*factorD
        else
           Nion(i,j)=Nion(i,j)*factorN      ! nightside diffusion
        endif

        ! at Southern ionosphere
        if (Sdon(i,j).eq.'d') then          ! dayside refilling
           tFNB=-dt*Fmax/Nsat(i,j)/Bi(i)
           factorD=exp(tFNB)
           Nion(i,j)=Nsat(i,j)-(Nsat(i,j)-Nion(i,j))*factorD
        else
           Nion(i,j)=Nion(i,j)*factorN      ! nightside diffusion
        endif

        ! limit Nion to Nsat
        if (Nion(i,j).gt.Nsat(i,j)) Nion(i,j)=Nsat(i,j)
     enddo
  enddo

  end subroutine refill_loss


!*******************************************************************************
  subroutine trough(nl,np)
!*******************************************************************************
! Routine makes sure the plasmasphere density is not lower than the trough
! density.
!
! Input: nl,np
! Input/Output: Nion

  implicit none

  integer nl,np,i,j


  do j=1,np
     do i=1,nl
        if (Nion(i,j).lt.NTro) Nion(i,j)=NTro
     enddo
  enddo

  end subroutine trough


!*******************************************************************************
      subroutine inter_flux(nl,np,ibp,cl,cp,fb1,f2,fal,fap)
!*******************************************************************************
!  Routine calculates the inter-flux, fal(i+0.5,j) and fap(i,j+0.5), using
!  2nd order flux limited scheme with super-bee flux limiter method.
!
!  Input: nl,np,ibp,cl,cp,fb1,f2
!  Output: fal,fap

   implicit none

   integer nl,np,i,j,j_1,j1,j2
   integer ibp(np),ibm
   real cl(nl,np),cp(nl,np),f2(nl,np),fal(nl,np),fap(nl,np),fb1, &
        fwbc(0:nl+2,np),xsign,fup,flw,x,r,xlimiter,corr

      fwbc(1:nl,1:np)=f2(1:nl,1:np)     ! fwbc is f2 with boundary condition

! Set up boundary condition
       fwbc(0,1:np)=f2(1,1:np)
       fwbc(nl+1:nl+2,1:np)=fb1

! find fa*
      do j=1,np
         j_1=j-1
         j1=j+1
         j2=j+2
         if (j_1.lt.1) j_1=j_1+np
         if (j1.gt.np) j1=j1-np
         if (j2.gt.np) j2=j2-np
         ibm=max(ibp(j),ibp(j1))
         do i=1,ibm    
            ! find fal
            xsign=sign(1.,cl(i,j))
            fup=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i+1,j)   ! upwind
            flw=0.5*(1.+cl(i,j))*fwbc(i,j)+0.5*(1.-cl(i,j))*fwbc(i+1,j)   ! LW
            x=fwbc(i+1,j)-fwbc(i,j)
            if (abs(x).le.1.e-27) fal(i,j)=fup
            if (abs(x).gt.1.e-27) then
               if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i-1,j))/x
               if (xsign.eq.-1.) r=(fwbc(i+2,j)-fwbc(i+1,j))/x
               if (r.le.0.) fal(i,j)=fup
               if (r.gt.0.) then
                  xlimiter=max(min(2.*r,1.),min(r,2.))
                  corr=flw-fup
                  fal(i,j)=fup+xlimiter*corr
                  if (fal(i,j).lt.0.) fal(i,j)=fup
               endif
            endif
            ! find fap
            xsign=sign(1.,cp(i,j))
            fup=0.5*(1.+xsign)*fwbc(i,j)+0.5*(1.-xsign)*fwbc(i,j1)   ! upwind
            flw=0.5*(1.+cp(i,j))*fwbc(i,j)+0.5*(1.-cp(i,j))*fwbc(i,j1)   ! LW
            x=fwbc(i,j1)-fwbc(i,j)
            if (abs(x).le.1.e-27) fap(i,j)=fup
            if (abs(x).gt.1.e-27) then
               if (xsign.eq.1.) r=(fwbc(i,j)-fwbc(i,j_1))/x
               if (xsign.eq.-1.) r=(fwbc(i,j2)-fwbc(i,j1))/x
               if (r.le.0.) fap(i,j)=fup
               if (r.gt.0.) then
                  xlimiter=max(min(2.*r,1.),min(r,2.))
                  corr=flw-fup
                  fap(i,j)=fup+xlimiter*corr
                  if (fap(i,j).lt.0.) fap(i,j)=fup
               endif
            endif
         enddo              ! end of do i=1,nl
      enddo                 ! end of do j=1,np

      end subroutine inter_flux
