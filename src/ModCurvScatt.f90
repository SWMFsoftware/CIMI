!*******************************************************************************
!                           ModCurvScatt.f90
! Module field line curvature (FLC) scattering.
!*******************************************************************************
  module ModCurvScatt
    use cimigrid_dim, only: ir,ip

    logical :: DoFlcLoss

    real,dimension(ir,ip) :: pc_flc,zeta1,zeta2
    !!real Daa_flc(ns,ir,ip,ie,0:ik+1),ekev_flc(ns,ir,ip,ie)
    !real Daa_flc(ir,ip,ie,0:ik+1),ekev_flc(ir,ip,ie)
    real tau_flc(ns,ir,ip)

  contains
!***********************************************************************
!                        calc_FLC_para      
!  Routine calculates parameters required to calculate pitch-angle
!   diffusion coefficient (Daa) due to field line curvature 
!   (curvature scattering).
!***********************************************************************
    ! for coupled model, this subroutine should take lat and lon index
    ! and the min curvature (Rc), associated d2Rc/ds2, d2B/ds2, and B
    ! in region around min B location (9 points)
    ! calculates zeta1 and zeta2 which are needed for scattering calcuation
    subroutine calc_FLC_para(i,j,npf,ibmin,req,dss,xa1,ya1,za1,bba)
 
  use constants, only: EM_speed,re_m 
  !!use CurvScatt, only: pc_flc,zeta1,zeta2
  implicit none
  
  integer,intent(in) :: i,j,npf,ibmin
  real,intent(in) :: req
  real,dimension(npf),intent(in) :: dss,xa1,ya1,za1,bba
  real,dimension(9) :: xa,ya,za,ds,Rc
  real Rc0,Rc1,Rc2,ds0,ds1,ds2,drC1,drC2,B0,dB1,dB2
  integer i10,i01,is,is1,is2

  ! (1) calculate curvature at bo
    i10=ibmin+4 ! high point around min b
    i01=ibmin-4 ! low point around min b
    xa(1:9)=xa1(i01:i10)
    ya(1:9)=ya1(i01:i10)
    za(1:9)=za1(i01:i10)
    ds(1:9)=dss(i01:i10)
    call calc_FLC(9,xa,ya,za,ds,Rc)

    !locate minimum of field line curvature
    is=minloc(Rc(2:8),dim=1)
    is=is+1
    Rc0=Rc(is)

    ! (2) calculate pc_flc,zeta1,zeta2 at bo
    ! pc_flc is the momentum which has a gyroradius that matches
    ! the minimum field line curvature
    ! Zeta1 and Zeta2 from Young et al 2002; 2008 are parameters
    ! which govern the calculation of the scattering
    ! note that the two curvature values for points around the
    ! min curvature are also saved for
    ! future calculation of second derivative of curvature
    pc_flc(i,j)=bba(ibmin)*EM_speed*Rc0*re_m*1.e-3  ! pc corresponding to FLC in keV
    is1=is-1
    is2=is+1
    Rc1=Rc(is1)
    Rc2=Rc(is2)
    ds1=(is1)
    ds2=(is2)
    ds0=0.5*(ds1+ds2)

    !When r less than 4Re (or whatever value) just use the dipole curvature scattering values
    !for zeta. Otherwise calculate 
    if (req.gt.4.) then   ! assume a dipole < 3 RE
       drC1=Rc0-Rc1
       drC2=Rc2-Rc0
     zeta1(i,j)=Rc0*(drC2/ds2-drC1/ds1)/ds0       ! Zeta1 as in Young et al. 2002; 2008
       B0=bba(is)
       dB1=B0-bba(is1)
       dB2=bba(is2)-B0           
     zeta2(i,j)=Rc0*Rc0/B0*(dB2/ds2-dB1/ds1)/ds0  ! Zeta2 as in Young et al. 2002; 2008
    else
       zeta1(i,j)=0.6666667
       zeta2(i,j)=1.
    endif

  end subroutine calc_FLC_para

  
!***********************************************************************
!                           calc_FLC      
!  Routine calculates field line curvature 
!***********************************************************************
    subroutine calc_FLC(n,xa,ya,za,ds,Rc)

    implicit none
    
    integer,intent(in) :: n
    real,dimension(n),intent(in) :: xa,ya,za,ds
    real,intent(out) :: Rc(n)
    integer i,i1
    real ds0,ds1,ds2,dx1,dx2,dy1,dy2,dz1,dz2,ddx,ddy,ddz

    ds1=ds(1)
    dx1=xa(2)-xa(1)
    dy1=ya(2)-ya(1)
    dz1=za(2)-za(1)
    do i=2,n-1
       i1=i+1
       ds2=ds(i)
       ds0=0.5*(ds2+ds1)
       dx2=xa(i1)-xa(i)
       dy2=ya(i1)-ya(i)
       dz2=za(i1)-za(i)
       ddx=abs(dx2/ds2-dx1/ds1)/ds0
       ddy=abs(dy2/ds2-dy1/ds1)/ds0
       ddz=abs(dz2/ds2-dz1/ds1)/ds0
       Rc(i)=1./(ddx+ddy+ddz)
       ds1=ds2
       dx1=dx2
       dy1=dy2
       dz1=dz2
    enddo
    Rc(1)=Rc(2)
    Rc(n)=Rc(n-1)

  end subroutine calc_FLC
  


!***********************************************************************
!                        FLC_loss     
!  Routine calculates characteristic lifetime due to 
!   field line curvature scattering and applies loss to f2
!***********************************************************************
  subroutine FLC_loss

  use cimigrid_dim, only: ns,ir,ip,iw,ik
  use cread1, only: rc ! loss cone altitude (could be input)
  use cread2, only: ijs,js,dt ! ijs= nSpecies, js=iSpecies could be inputs
  use cinitial, only: f2
  use cfield, only: iba,ekev,Tbounce,y,rmir !y=sinA this would ModFieldTrace
                                            ! Tbounce needs to be saved in coupled model
  !!use ModCurvScatt, only: ie,pc_flc,zeta1,zeta2,Daa_flc,ekev_flc
  implicit none
  
  real ekev0,Tbounce0(iw,ik),sina,cosa,a0
  real E02(4),pc,&
       e,e2,e3,c,a1,a2,De,w,b,Amax,tb,D0,Ni0,Dmax,&
       Daa_flc1,tau_flc1
  integer n,i,j,k,m,m_lc,k1

  !this should be set in ModPlanet.f90 and used here
  E02(1)=938272.0813*2      ! p+ rest mass*2
  E02(2)=16.*E02(1)         ! O+ rest mass*2
  E02(3)=4.*E02(1)          ! He+ rest mass*2
  E02(4)=1021.997692        ! e- rest mass*2 in keV


  do j=1,ip
  do i=1,iba(j)

! find the index m adjacent to loss cone 
  do m=ik,1,-1
     if (rmir(i,j,m).ge.rc) exit
  enddo
  m_lc=m+1
  sina=y(i,j,m_lc)
  a0=asin(sina)
  cosa=cos(a0)

  tau_flc(:,i,j)=864000.  !! for print only

  !for each species calculate the coefficients as in Young et al
  ! to determine the scattering liftime and apply the loss
  do n=1,ijs
     if (js(n).ge.4.and.js(n).le.5) then
        do k=1,iw
           if (ekev(n,i,j,k,m_lc).ge.1000.) exit
        enddo
     else
        do k=1,iw
           if (ekev(n,i,j,k,m_lc).ge.10.) exit
        enddo
     endif
     k1=k
     do k=1,iw
        ekev0=ekev(n,i,j,k,m_lc)
        pc=sqrt(ekev0*(ekev0+E02(js(n))))
        tb=Tbounce(n,i,j,k,m_lc)
        ! coefficients as in Young et al. 2002; 2008
        ! e=RatioGyroToCurvature or RatioGyroCurv
        e=pc/pc_flc(i,j)                         ! epsilon
        ! Ratio of Gyroradius to curvature is assumed to be limited
        ! as in Young et al as they assume slow scattering limit
        ! for applicability of diffusion approach. We assume
        ! same limit as in that paper.
        if (e.gt.0.584) e=0.584
        if (e.gt.0.05) then
           e2=e*e
           e3=e2*e
           !these variable names are identical to Young et al paper
           ! but Amax would be A0
           w=1.051354+0.13513581*e-0.50787555*e2        ! omega
           c=1.0663037-1.0944973/e+0.016679378/e2-0.00049938987/e3
           a1=-0.35533865+0.12800347/e+0.0017113113/e2    
           a2= 0.23156321+0.15561211/e-0.001860433/e2    
           b=-0.51057275+0.93651781/e-0.0031690066/e2    
           De=-0.49667826-0.00819418/e+0.0013621659/e2  ! D(epsilon)
           Amax=exp(c)*(zeta1(i,j)**a1*zeta2(i,j)**a2+De)

           ! now use coef. to calculate lifetime and apply loss
           Ni0=sin(w*a0)*cosa**b
           D0=0.5*Amax*Amax/tb
           Dmax=4.9348/tb    ! da0=0.5pi during 0.5*tb
           Daa_flc1=D0/Ni0/Ni0*(sin(w*a0)/sina)**2*cosa**(2.*b-2.)
           if (Daa_flc1.gt.Dmax) Daa_flc1=Dmax
           tau_flc1=0.5*abs(log(sina))/Daa_flc1    ! life time
           f2(n,i,j,k,:)=f2(n,i,j,k,:)*exp(-dt/tau_flc1)
        else 
           w=0.
           b=0.
           Amax=0. 
           tau_flc1=864000.
        endif
        
        if (e>0.05) then
        
        else

        endif
        if (k.eq.k1) tau_flc(n,i,j)=tau_flc1 
     enddo    ! end of k

  enddo       ! end of n
  enddo       ! end of i
  enddo       ! end of j

  end subroutine FLC_loss

!***********************************************************************
!                           write_FLC     
! 
!  Routine prints ratio of 10 keV proton gyroradius to
!   field line curvature radius
!***********************************************************************
  subroutine write_FLC(t,tstep)
  use cread1, only:outname
  use cread2, only: tstart,ijs,js
  use cfield, only: ro,xmlto
  !!use ModCurvScatt, only: pc_flc
  implicit none

  real,intent(in) :: t,tstep
  real E02(4),E1,pc
  integer n
  character outname1*40

  E1=10.     ! in keV

  E02(1)=938272.0813*2      ! p+ rest mass*2
  E02(2)=16.*E02(1)         ! O+ rest mass*2
  E02(3)=4.*E02(1)          ! He+ rest mass*2
  E02(4)=1021.997692        ! e- rest mass*2 in keV

  pc=sqrt(E1*(E1+E02(1)))   ! pc for E1 proton
  
  if (t.eq.tstart) then 
     call system('mkdir -p FLC')
     open(unit=61,file='FLC/'//trim(outname)//'.parm',&
          form='unformatted',status='replace')
     write(61) tstart,tstep,ir,ip,ijs
     write(61) js(1:ijs)
     close(61)
  endif
  write(outname1,'(a,i8.8,a)') 'FLC/'//trim(outname)//&
        '_',int(t),'.flc'
  open(unit=60,file=trim(outname1),&
       form='unformatted',status='replace')
  write(60) t
  write(60) ro(:,:)
  write(60) xmlto(:,:)
  write(60) pc/pc_flc(:,:)
  do n=1,ijs
     write(60) tau_flc(n,:,:)
  enddo
  close(60)
  open(unit=62,file='FLC/'//trim(outname)//'.end',&
       form='unformatted',status='replace')
  write(62) t
  close(62)
  
  
  end subroutine write_FLC

  end module           
  
!  !***********************************************************************
!  !                        calc_Daa_FLC      
!  !  Routine calculates parameters required to calculate pitch-angle
!  !   diffusion coefficient (Daa) due to field line curvature 
!  !   (curvature scattering).
!  !***********************************************************************
!  ! This subroutine not used in coupled model, but may save for the future
!  ! Daa calculation
!  subroutine calc_Daa_FLC(i,j,ekev0,Tbounce0,sina)
!
!  use cimigrid_dim, only: iw,ik
!  !!use ModCurvScatt, only: ie,pc_flc,zeta1,zeta2,Daa_flc,ekev_flc
!  implicit none
!  
!  integer,intent(in) :: i,j
!  real,intent(in) :: ekev0(iw,ik),Tbounce0(iw,ik),sina(0:ik+1)
!  real E02,emin,emax,emin_log,emax_log,dlogE,ein_log,ein(ie),pc,&
!       e,e2,e3,e_1,e_2,c,a1,a2,De,b0,w0,tb,D0,Ni0,Amx0,&
!       cosa(0:ik+1),a0(0:ik),cosa_bar,Dmax
!  real,dimension(ie,ik) :: w,b,Amax
!  integer k,m,m_bar
!
!  E02=1021.997692      ! e- rest mass*2 in keV
!
!  emin=minval(ekev0(1 ,1:ik))
!  emax=maxval(ekev0(iw,1:ik))
!  emin_log=log10(emin)
!  emax_log=log10(emax)
!  dlogE=(emax_log-emin_log)/(ie-1)
!  do k=1,ie
!     ein_log=emin_log+float(k-1)*dlogE
!     ein(k)=10.**ein_log
!     ekev_flc(i,j,k)=ein(k)
!  enddo
!
!  do m=0,ik+1
!     a0(m)=asin(sina(m))
!     cosa=cos(a0(m))
!  enddo
!
!  do m=1,ik
!     do k=1,ie
!        pc=sqrt(ein(k)*(ein(k)+E02))
!        ! coefficients as in Young et al. 2002; 2008
!        e=pc/pc_flc(i,j)                         ! epsilon
!        if (e.gt.0.584) e=0.584
!        if (e.gt.0.05) then
!           e2=e*e
!           e3=e2*e
!           w(k,m)=1.051354+0.13513581*e-0.50787555*e2        ! omega
!           c=1.0663037-1.0944973/e+0.016679378/e2-0.00049938987/e3
!           a1=-0.35533865+0.12800347/e+0.0017113113/e2    
!           a2= 0.23156321+0.15561211/e-0.001860433/e2    
!           b(k,m)=-0.51057275+0.93651781/e-0.0031690066/e2    
!           De=-0.49667826-0.00819418/e+0.0013621659/e2  ! D(epsilon)
!           Amax(k,m)=exp(c)*(zeta1(i,j)**a1*zeta2(i,j)**a2+De)
!        else 
!           w(k,m)=0.
!           b(k,m)=0.
!           Amax(k,m)=0. 
!        endif
!     enddo
!  enddo 
!
!  do k=1,ie
!     m_bar=maxloc(Amax(k,:),dim=1)
!     cosa_bar=cosa(m_bar)
!     do m=1,ik
!        Amx0=Amax(k,m)
!        if (Amx0.ne.0.) then
!           w0=w(k,m)
!           b0=b(k,m)
!           Ni0=sin(w0*a0(m_bar))*cosa_bar**b0
!           call lintp(ekev0(:,m),Tbounce0(:,m),iw,ein(k),tb)    ! tau_b
!           D0=0.5*Amx0*Amx0/tb
!           Dmax=4.9348/tb    ! da0=0.5pi during 0.5*tb
!           Daa_flc(i,j,k,m)=D0/Ni0/Ni0*(sin(w0*a0(m))/sina(m))**2*cosa(m)**(2.*b0-2.)
!           if (Daa_flc(i,j,k,m).gt.Dmax) Daa_flc(i,j,k,m)=Dmax
!        else
!           Daa_flc(i,j,k,m)=0.
!        endif
!     enddo
!     Daa_flc(i,j,k,0)=Daa_flc(i,j,k,1)
!     Daa_flc(i,j,k,ik+1)=Daa_flc(i,j,k,ik)
!  enddo
!
!  end subroutine calc_Daa_FLC
