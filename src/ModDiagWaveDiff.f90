!*******************************************************************************
!
!                               ModDiagDiff.f90  
!
! This module contains subroutines that calculates momentum space diffusion
!  in (Q1, Q2) coordinates.
!  where (Q1, Q2) is found such that cross-diffusion coef., DQ1Q2 and DQ2Q1 are
!  zero [Albert and Young, 2005].
!
!  Q1 = K
!  const Q2 curves: dE/da0 = DaE/Daa
!
!  D_Q1Q1 = DKK = Daa(dK/da)^2
!  D_Q2Q2 = (DEE - DaE^2/Daa)(dQ2/dE)^2
! 
! Created by Suk-Bin Kang, NASA/GSFC Geospace Physics Lab, in Jan - Mar 2018.
!
!*******************************************************************************

  module diagDiffCoef 
        use cimigrid_dim, only: ir,ip,ik
        integer,parameter :: iq=80
        real VarQ(0:iq+1)
        real,allocatable,dimension(:,:,:) :: E_Q2c,E_Q2h,&
             cDqq1,cDqq2,hDqq1,hDqq2,cSign,hSign
        real,dimension(ir,ip,iq,0:ik+1) :: Dqq1K,Dqq2K,dQ2dEK,E_Q2K
        real PSD_Q(ir,ip,iq,ik)
        integer ipc0(ir,ip),iph0(ir,ip),iLpp(ip)

  contains
! ************************************************************************
!                            diffuse_Q2
!  Routine calculates the change of electron distributions due to
!  diffusion in Q2.
! ************************************************************************
  subroutine diffuse_Q2
  use constants, only: e_mass,EM_speed,echarge
  use cimigrid_dim, only: ir,ip,ik
  use cgrid, only: xjac
  use cinitial, only: f2
  use cread2, only: dt,js,ijs,ichor,ihiss
  use cfield, only: iba,bo,ro
  use cWpower, only: CHpower,HIpower,&
                     Cpower0,Hpower0,BLc0,BLh0
  !!use diagDiffCoef, only: iq,VarQ,Dqq2K,Dqq2K,dQ2dEK,E_Q2K,PSD_Q,&
  !!                        ipc0,iph0,iLpp
  implicit none

  real,dimension(iq,ik) :: Dqq2,dQ2dE,E1,PSD
  real,dimension(iq) :: G_Q2,Dqq2m,Dqq2g
  real,dimension(iq) :: b1,b2,b3,b4,fnew
  real dU_Q2,jac0,E0,E02,ompe1,Po
  integer n,i,j,k,m,nel
  logical DoDiff

!(1) Find the "n" corresponds to electrons, nel 
  nel=0
  do n=1,ijs
     if (js(n).eq.4) nel=n
  enddo

!(2) Determine dU_Q2
  dU_Q2=log(VarQ(2)/VarQ(1))
  jac0=dt/dU_Q2/dU_Q2
  E0=e_mass*EM_speed**2/echarge*1.e-3 ! in keV
  E02=2.*E0

  do j=1,ip             ! start of j
     do i=1,iba(j)      ! start of i
        DoDiff=.False.
        if (i.gt.iLpp(j).and.ichor.ge.1) then
!(3) outside plasmasphere
           Po=(BLc0/bo(i,j))*CHpower(i,j)/Cpower0
           DoDiff=.True.
        else if (ihiss.ge.1) then
!(4) inside plasmasphere
           Po=(BLh0/bo(i,j))*HIpower(i,j)/Hpower0
           DoDiff=.True.
        endif
        if (DoDiff) then
           Dqq2(1:iq,1:ik)=Dqq2K(i,j,1:iq,1:ik)*Po
           dQ2dE(1:iq,1:ik)=dQ2dEK(i,j,1:iq,1:ik)
           E1(1:iq,1:ik)=E_Q2K(i,j,1:iq,1:ik)
           PSD(1:iq,1:ik)=PSD_Q(i,j,1:iq,1:ik)
           do m=1,ik      ! start of m
!(5) calulate Jacobian G_Q2
              G_Q2(:)=(E1(:,m)+E0)*sqrt(E1(:,m)*(E1(:,m)+E02))/dQ2dE(:,m)*VarQ(1:iq)
!(6) calculate Dq2q2 at each mid-grid point
              Dqq2g(:)=Dqq2(:,m)*G_Q2(:)
              do k=1,iq-1   ! start of k
                 Dqq2m(k)=0.5*(Dqq2g(k)+Dqq2g(k+1))    ! arithmetic mean
              enddo       ! end of k
!(7) Set BTCS (Backward Time and Cetral Space) variables
              do k=2,iq-1   ! start of k
                 b1(k)=-Dqq2m(k-1)/G_Q2(k)*jac0
                 b3(k)=-Dqq2m(k  )/G_Q2(k)*jac0
                 b2(k)=1.-b1(k)-b3(k)
                 b4(k)=PSD(k,m)
              enddo       ! end of k
!(8) Apply boundary condition
              b4(2)=b4(2)-b1(2)*PSD(1,m)            ! const. at lowest V
              b4(iq-1)=b4(iq-1)-b3(iq-1)*PSD(iq,m)  ! const. at highest V
!(9) Compute diffusion using Thomas Algorithm
              call tridiagonal(b1(2:iq-1),b2(2:iq-1),b3(2:iq-1),b4(2:iq-1),iq-2,&
                               fnew(2:iq-1))
!(10) Make sure f2 > 0
              do k=2,iq-1   ! start of k
                 if (fnew(k).gt.1.e90) then
                    write(*,*) ' PSD > 1.e90 !'
                    write(*,'(a,4i3)') 'i,j,k,m =',i,j,k,m
                    write(*,'(1p3E15.7)') PSD(k-1:k+1,m)
                    write(*,'(1pE15.7)') fnew(k)
                    write(*,'(1p2E15.7)') Dqq2m(k-1:k)/G_Q2(k)
                    stop
                 endif
                 if (fnew(k).ne.fnew(k).or.fnew(k).lt.-1.e-50) then
                     write(*,'(a,i3,a,i3,a,i3,a,i3,a)') &
                     'At i =',i,' j=',j,', k=',k,', m =',m,' in diffusion_Q2'
                     write(*,*) fnew(k)
                     write(*,'(1p10E11.3)') fnew(2:iq-1)
                     write(*,*) ' f2(t-1)'
                     write(*,'(1p10E11.3)') PSD(:,m)
                     write(*,*) ' D_Q2Q2(t)'
                     write(*,'(1p10E11.3)') Dqq2m(1:iq-1)/G_Q2(1:iq-1)
                     stop
                 endif
              enddo       ! end of k
!(11) Update PSD_Q
              PSD_Q(i,j,2:iq-1,m)=fnew(2:iq-1)
           enddo          ! end of m
        endif
     enddo             ! end of i
  enddo                ! end of j
  
  end subroutine diffuse_Q2


! ************************************************************************
!                            diffuse_Q1
!  Routine calculates the change of electron distributions due to
!  diffusion in Q1(=K).
! ************************************************************************
  subroutine diffuse_Q1
  use constants, only: e_mass,EM_speed,echarge,re_m
  use cimigrid_dim, only: ir,ip,ik
  use cgrid, only: xjac,xk
  use cinitial, only: f2
  use cread2, only: dt,js,ijs,ichor,ihiss
  use cfield, only: iba,bo,bm,ro,tya,y
  use cWpower, only: CHpower,HIpower,&
                     Cpower0,Hpower0,BLc0,BLh0
  !!use diagDiffCoef, only: iq,Dqq1K,PSD_Q,&
  !!                        ipc0,iph0,iLpp
  implicit none

  real,dimension(iq,ik) :: PSD
  real,dimension(iq,0:ik+1) :: Daa1
  real,dimension(0:ik+1) :: VarK,G_a,bm1,G_Q1,dKda_K,Dqq1m,Dqq1g
  real,dimension(ik) :: b1,b2,b3,b4,fnew
  real dU_Q1,jac0,Po,cosa,sina
  integer n,i,j,k,m,nel
  logical DoDiff

!(1) Find the "n" corresponds to electrons, nel 
  nel=0
  do n=1,ijs
     if (js(n).eq.4) nel=n
  enddo

!(2) Determine dU_Q2
  VarK(0:ik+1)=xk(0:ik+1)/re_m  ! K in T^0.5 RE
  dU_Q1=log(VarK(2)/VarK(1))
  jac0=dt/dU_Q1/dU_Q1

  do j=1,ip             ! start of j
     do i=1,iba(j)      ! start of i
        DoDiff=.False.
        if (i.gt.iLpp(j).and.ichor.ge.1) then
!(3) outside plasmasphere
           Po=(BLc0/bo(i,j))*CHpower(i,j)/Cpower0
           DoDiff=.True.
        else if (ihiss.eq.1) then
!(4) inside plasmasphere
           Po=(BLh0/bo(i,j))*HIpower(i,j)/Hpower0
           DoDiff=.True.
        endif
        if (DoDiff) then
           Daa1(1:iq,0:ik+1)=Dqq1K(i,j,1:iq,0:ik+1)*Po
           PSD(1:iq,1:ik)=PSD_Q(i,j,1:iq,1:ik)
!(5) calulate Jacobian G_Q1
           bm1(1:ik)=bm(i,j,1:ik)
           bm1(0)=bo(i,j)/y(i,j,0)**2
           bm1(ik)=bo(i,j)/y(i,j,ik+1)**2
           do m=0,ik+1
              sina=y(i,j,m)
              cosa=sqrt(1.-sina*sina)
              G_a(m)=tya(i,j,m)*sina*cosa                            ! T(a0)sin2a0/2
              dKda_K(m)=sqrt(bm1(m))*cosa/sina*tya(i,j,m)/VarK(m) ! -dK/da/K
              G_Q1(m)=G_a(m)/dKda_K(m)
           enddo
           do k=1,iq        ! start of k
!(6) calculate Dq1q1 at each mid-grid point
              Dqq1g(0:ik+1)=Daa1(k,0:ik+1)*G_a(0:ik+1)*dKda_K(0:ik+1) ! DKK*G_Q1*K2
              do m=0,ik     ! start of m
                 Dqq1m(m)=0.5*(Dqq1g(m)+Dqq1g(m+1))              ! arithmetic mean
              enddo         ! end of m
!(7) Set BTCS (Backward Time and Cetral Space) variables
              do m=1,ik     ! start of k
                 b1(m)=-Dqq1m(m-1)/G_Q1(m)*jac0
                 b3(m)=-Dqq1m(m)/G_Q1(m)*jac0
                 b2(m)=1.-b1(m)-b3(m)
                 b4(m)=PSD(k,m)
              enddo         ! end of m
!(8) Apply boundary condition
              b2(1)=b1(1)+b2(1)                     ! zero gradient at lowest K
              !!b4(ik)=b4(ik)-b3(ik)*PSD(k,ik)        ! const. at highest K
              b4(ik-1)=b4(ik-1)-b3(ik-1)*PSD(k,ik)  ! const. at highest K
!(9) Compute diffusion using Thomas Algorithm
              !!call tridiagonal(b1,b2,b3,b4,ik,fnew)
              call tridiagonal(b1(1:ik-1),b2(1:ik-1),b3(1:ik-1),&
                               b4(1:ik-1),ik-1,fnew(1:ik-1))
              fnew(ik)=PSD(k,ik)
!(10) Make sure f2 > 0
              do m=1,ik     ! start of m
                 if (fnew(m).ne.fnew(m).or.fnew(m).lt.-1.e-50) then
                     write(*,'(a,i3,a,i3,a,i3,a,i3,a)') &
                     'At i =',i,' j=',j,', k=',k,', m =',m,' in diffusion_Q1'
                     write(*,*) fnew(m)
                     write(*,'(1p10E11.3)') fnew(1:ik)
                     write(*,*) ' f2(t-1) in #s3/m6/kg3'
                     write(*,'(1p10E11.3)') PSD(k,:)
                     write(*,*) ' DKK(t) in 1/s'
                     write(*,'(1p10E11.3)') Dqq1m(0:ik)/G_Q1(0:ik)
                     write(*,*) ' a0(t) in deg'
                     write(*,'(10f11.4)') asin(y(i,j,0:ik+1))*57.2957795
                     write(*,*) ' Daa(t) in 1/s'
                     write(*,'(1p10E11.3)') Daa1(k,0:ik+1)
                     write(*,*) ' T(y) in RE'
                     write(*,'(1p10E11.3)') tya(i,j,0:ik+1)
                     stop
                 endif
              enddo       ! end of k
!(11) Update PSD_Q
              PSD_Q(i,j,k,1:ik)=fnew(1:ik)
           enddo          ! end of k
        endif
     enddo             ! end of i
  enddo                ! end of j
  
  end subroutine diffuse_Q1


! ************************************************************************
!                          tridiagonal      
!
!  Subroutine solves inversion of tridiagonal matrix 
!   through Thomas's tridiagonal matrix algorithm.
!
!  b1*fi-1 + b2*fi + b3*fi+1 = b4,  i= 1,2,3,..,n 
!
!  Inputs : b1,b2,b3,b4,n
!  Outputs : f
! ************************************************************************
      subroutine tridiagonal(b1,b2,b3,b4,n,f)
      implicit none

      integer,intent(in) :: n
      real,intent(in) :: b1(n),b2(n),b3(n),b4(n)
      real,intent(out) :: f(n)
      real c0,c1(n),c2(n)
      integer i
 
      c1(1)=b3(1)/b2(1)
      c2(1)=b4(1)/b2(1)
      do i=2,n
         c0=b2(i)-b1(i)*c1(i-1)
         c1(i)=b3(i)/c0
         c2(i)=(b4(i)-b1(i)*c2(i-1))/c0
      enddo      

      f(n)=c2(n)
      do i=n-1,1,-1
         f(i)=c2(i)-c1(i)*f(i+1)
      enddo

      end subroutine tridiagonal


! ************************************************************************
!                        find_ompe_index       
!
! Subroutine finds indices of closest ompe 
! ************************************************************************
  subroutine find_ompe_index
  use cimigrid_dim, only: ir,ip,ipc,iph
  use cread2, only: ichor,ihiss
  use cfield, only: iba
  use waveDiffCoef, only: cOmpe,hOmpe
  use cWpower, only: ompe
  !!use diagDiffCoeff, only:ipc0,iph0
  implicit none

  real ompe1
  integer i,j,ipc1,iph1

  do j=1,ip             ! start of j
     do i=1,iba(j)      ! start of i
        ompe1=ompe(i,j)                         ! fpe/fce
        if (ichor.ge.1) then
           call locate1(cOmpe,ipc,ompe1,ipc1)
           if (ipc1.eq.0) ipc1=1
           if (ipc1.gt.ipc) ipc1=ipc
           if ((cOmpe(ipc1+1)-ompe1).lt.(ompe1-cOmpe(ipc1))) ipc1=ipc1+1
           ipc0(i,j)=ipc1
        endif
        if (ihiss.ge.1) then
           call locate1(hOmpe,iph,ompe1,iph1)
           if (iph1.eq.0) iph1=1
           if (iph1.gt.iph) iph1=iph
           if ((hOmpe(iph1+1)-ompe1).lt.(ompe1-hOmpe(iph1))) iph1=iph1+1
           iph0(i,j)=iph1
        endif
     enddo              ! end of i
  enddo                 ! end of j
 
  end subroutine find_ompe_index
  
  
!*****************************************************************************
!                           calc_Dqq    
!  Routine find grids (Q1,Q2) and calculates diffusion coefficients in (Q1,Q2).
!  DQQ = (Daa - DaE^2/Daa)(dQ/dE)**2 = (DKK - 2DKE^2/DKK)*(dQ/dE)**2
!
!  Q1=a0 (or K)
!  Q2 is chosen such that Dq1q2=0 [Albert and Young 2005; Camperol et al. 2013]
!
!  Note
!  All Daa,DaE,DEE due to all plasma waves at the same location should be 
!   taken into account, to determine (Q1,Q2)
!*****************************************************************************
  subroutine calc_DQQ
  use constants, only: pi
  use cimigrid_dim, only: ir,ip,ik,ipc,iwc,iph,iwh,ipa
  use cread2, only: ichor,ihiss
  use waveDiffCoef, only: ckeV,hkeV,cDEE,cDaa,cDaE,hDEE,hDaa,hDaE,cPA,cOmpe
  !!use diagDiffCoef, only: iq,VarQ,E_Q2c,E_Q2h,cDqq1,cDqq2,hDqq1,hDqq2,&
  !!                        cSign,hSign,Dqq1K,Dqq2K,dQ2dEK,E_Q2K
  implicit none

  real Qmin,Qmax,rQ,Eq(0:iq+1),logEq(0:iq+1)
  real,dimension(ipc,0:iq+1,ipa) :: cDaa0,cDaE0,cDEE0
  real,dimension(iph,0:iq+1,ipa) :: hDaa0,hDaE0,hDEE0
  integer j,k,m,k1,k2

!(1) Set variable Q2
   Qmin=1.      ! at k = 1  
   Qmax=1.e4    ! at k = iq+1
   rQ=10.**(log10(Qmax/Qmin)/iq)
   VarQ(0)=Qmin/rQ
   do k=1,iq+1
      VarQ(k)=VarQ(k-1)*rQ
   enddo
   Eq(0:iq+1)=VarQ(0:iq+1)
   logEq(:)=log(Eq(:))


  if (ichor.ge.1) then
     allocate (cSign(ipc,iwc,ipa))
     cSign(:,:,:)=sign(1.,cDaE(:,:,:))
     where (cDaE(:,:,:).eq.0.) cSign(:,:,:)=1.

!(2) interpolate diffusion coef. to same grids of Q
     call interpol_D_coef(cDaa,cDaE,cDEE,cSign,ckeV,logEq,ipc,iwc,iq,cDaa0,cDaE0,cDEE0)
     allocate (E_Q2c(ipc,0:iq+1,ipa))
     allocate (cDqq1(ipc,iq,ipa),cDqq2(ipc,iq,ipa))

!(3) calculate energy crossponding to constant Q
     call calc_Q_curve(cDaa0,cDaE0,Eq,cPA*pi/180,ipc,0,iq+1,ipa,E_Q2c)

!(4) calculate Daa and Dqq2 at the fixed grids of (a0,Q2)
     call DaEtoDQQ(cDaa0,cDEE0,cDaE0,Eq,E_Q2c,VarQ,ipc,iq,iq,cDqq1,cDqq2)
  endif
  if (ihiss.eq.1) then
     allocate (hSign(iph,iwh,ipa))
     hSign(:,:,:)=sign(1.,hDaE(:,:,:))
     where (hDaE(:,:,:).eq.0.) hSign(:,:,:)=1.
     call interpol_D_coef(hDaa,hDaE,hDEE,hSign,hkeV,logEq,iph,iwh,iq,hDaa0,hDaE0,hDEE0)
     call calc_Q_curve(hDaa0,hDaE0,Eq,cPA*pi/180,iph,0,iq+1,ipa,E_Q2h)
     call DaEtoDQQ(hDaa0,hDEE0,hDaE0,Eq,E_Q2h,VarQ,ipc,iq,iq,hDqq1,hDqq2)
  endif
  
  call write_Q2info(cDEE,cDaE,cDaa0,cDaE0,cDEE0,cDqq1,cDqq2,ckeV,Eq,E_Q2c,ipc,iwc,iq)
  call interpol_D_coefK
  
  end subroutine calc_DQQ


!*****************************************************************************
!                           calc_Q_curve    
!  Routine calculates energy E(a,Q) corresponding to contant Q curve.
!*****************************************************************************
  subroutine calc_Q_curve(Daa,DaE,keV,a0,ip,iw0,iw,ipa,E)
  use constants, only: pi
  implicit none

  integer,intent(in) :: ip,iw0,iw,ipa
  real,dimension(ip,iw0:iw,ipa),intent(in) :: Daa,DaE
  real,intent(in) :: keV(iw0:iw),a0(ipa)
  real,intent(out) :: E(ip,iw0:iw,ipa)
  real,allocatable :: da(:),keV1(:),dEda(:,:)
  real r,dEdam,dEda1
  integer j,k,m,m0,k1,k2,kk,iww

  iww=iw-iw0+1
  allocate (da(ipa-1))
  allocate (keV1(iww),dEda(iww,ipa))
 
  r=0.99  ! minimum ratio difference between adjacent grids
  do m=1,ipa-1
     da(m)=a0(m+1)-a0(m)
  enddo
  keV1(1:iww)=keV(iw0:iw)

  do j=1,ip
     dEda(:,:)=0.
     where(Daa(j,:,:).ne.0.) &
        dEda(1:iww,1:ipa)=DaE(j,iw0:iw,1:ipa)/Daa(j,iw0:iw,1:ipa)
     do m=1,ipa
        dEda(:,m)=dEda(:,m)*keV1(:)
     enddo
     do k=iw0,iw
        call runge_kutta(keV(k),da,keV1,a0,dEda(:,:),iww,ipa,E(j,k,:))
     enddo
! make sure E is monotonically increasing
     do m=2,ipa
        m0=m-1
        do k=iw-1,iw0,-1
           k2=k+1
           if (E(j,k,m).ge.E(j,k2,m)) then
              write(*,*) ' E(a0,Q) calculation'
              write(*,'(a,3i3,a)') ' At j,k,m =',j,k,m,',  E(k-1)>E(k)'
              write(*,'(a,1p2E11.3)') ' E(j,k,1),E(j,k+1,1) =',E(j,k:k2,1)
              write(*,'(a,1p2E11.3)') ' E(j,k,m),E(j,k+1,m) =',E(j,k:k2,m)
              write(*,'(a,1p2E11.3)') ' DaE/Daa(k,m-1:m) (in keV) =',dEda(k,m0:m)*keV(k)
              write(*,'(a,1p2E11.3)') ' DaE/Daa(k+1,m-1:m) (in keV) =',dEda(k+1,m0:m)*keV(k+1)
              if (E(j,k2,m)/E(j,k,m).ge.0.97) then
                 !!E(j,k,m)=0.97*E(j,k2,m)
              else
                 stop
                 write(*,'(f6.2,a)') E(j,k2,m)/E(j,k,m)*100.,' %'
              endif
              stop
              !!E1=E(j,k2,m)*r
              !!E(j,k,m)=E1
           endif
        enddo
     enddo
  enddo
  
  deallocate (da,keV1,dEda)

  end subroutine calc_Q_curve


!*****************************************************************************
!                           DaEtoDQQ  
!  Routine calculates DQQ(a0,Q) at fixed grids of (a0,Q)
!   from Daa,DEE,DaE at fixed grids of (a0,E)
!*****************************************************************************
     subroutine DaEtoDQQ(Daa,DEE,DaE,keV,E,VarQ,ip,iw,iq,Dqq1,Dqq2)
     use cimigrid_dim, only: ipa
     implicit none

     integer,intent(in) :: ip,iw,iq
     real,dimension(ip,0:iw+1,ipa),intent(in) :: Daa,DEE,DaE
     real,dimension(ip,0:iq+1,ipa),intent(in) :: E
     real,intent(in) :: keV(0:iw+1),VarQ(0:iq+1)
     real,dimension(ip,iq,ipa),intent(out) :: Dqq1,Dqq2
     real,dimension(iw+2) ::  keVlog,logDaa,logDEE,logDaE
     real,dimension(0:iq+1) :: logE,VarQ2,E2
     real dQ(iq),Daa1,DaE1,DEE1,dQdE,DaE2,DaE2max,Dqq0,&
          logDaa1,logDaE1,logDEE1
     integer k1,iw1,j,k,m

     Dqq1(:,:,:)=0.
     Dqq2(:,:,:)=0.
     keVlog(1:iw+2)=log(keV(0:iw+1))
     do k=1,iq
        dQ(k)=VarQ(k+1)-VarQ(k-1)
     enddo
     VarQ2(0:iq+1)=VarQ(0:iq+1)*VarQ(0:iq+1) 

     do j=1,ip
        do m=1,ipa
           E2(0:iq+1)=E(j,0:iq+1,m)*E(j,0:iq+1,m)
           logE(0:iq+1)=log(E(j,0:iq+1,m))
           logDaa(:)=-75.
           logDaE(:)=-75.
           logDEE(:)=-75.
           where(Daa(j,0:iw+1,m).ne.0.) &
                logDaa(1:iw+2)=log(Daa(j,0:iw+1,m))
           where(DaE(j,0:iw+1,m).ne.0.) &
                logDaE(1:iw+2)=log(abs(DaE(j,0:iw+1,m))) !  in s^-1
           where(DEE(j,0:iw+1,m).ne.0.) &
                logDEE(1:iw+2)=log(DEE(j,0:iw+1,m))      !  in s^-1
           do k=1,iq
              if (logE(k).ge.keVlog(1).and.logE(k).le.keVlog(iw+2)) then
                 call lintp(keVlog(:),logDaa(:),iw+2,logE(k),logDaa1)
                 call lintp(keVlog(:),logDaE(:),iw+2,logE(k),logDaE1)
                 call lintp(keVlog(:),logDEE(:),iw+2,logE(k),logDEE1)
                 if (logDaa1.gt.-70.) then
                    Daa1=exp(logDaa1)
                 else
                    Daa1=0.
                 endif
                 if (logDEE1.gt.-70.) then
                    DEE1=exp(logDEE1)
                 else
                    DEE1=0.
                 endif
                 if (Daa1.gt.0..and.DEE1.gt.0.) then
                    if (DaE1.gt.-70.) then
                       DaE1=exp(logDaE1)
                    else
                       DaE1=0.
                    endif
                 else
                    DaE1=0.
                 endif
                 Dqq1(j,k,m)=Daa1
                 if (Daa1.eq.0.) then
                    Dqq2(j,k,m)=DEE1*E2(k)/VarQ2(k)
                 else
                    DaE2=DaE1*DaE1
                    DaE2max=Daa1*DEE1
                    if (DaE2.gt.DaE2max) DaE2=DaE2max
                    DQQ0=(DaE2max-DaE2)/Daa1
                    dQdE=dQ(k)/(E(j,k+1,m)-E(j,k-1,m))
                    Dqq2(j,k,m)=DQQ0*dQdE*dQdE*E2(k)/VarQ2(k) ! Dqq2/Q^2 in 1/s
                 endif
                 if (Dqq2(j,k,m).lt.0.) then
                    write(*,'(a,3i3,a,1p5E11.3)')&
                     ' At j,k,m =',j,k,m,', DQQ,Daa,DEE,DaE =',&
                     Dqq2(j,k,m),Daa1,DEE1,DaE1,DEE1-DaE1*DaE1/Daa1
                    stop
                 endif
              endif
           enddo
        enddo
     enddo

     end subroutine DaEtoDQQ


!*****************************************************************************
!                         interpol_D_coef
!  Routine interpolates Daa,DaE,DEE to Q2(m=1) energy grids
!   where Q2 = E at a0 = a0(m=1)
!*****************************************************************************
     subroutine interpol_D_coef(Daa,DaE,DEE,xSign,keV,logE,ip,iw,iq,Daa0,DaE0,DEE0)
     use waveDiffCoef, only: cPA,ipa
     implicit none

     integer,intent(in) :: ip,iw,iq
     real,dimension(ip,iw,ipa),intent(in) :: Daa,DaE,DEE,xSign
     real,intent(in) :: keV(iw),logE(0:iq+1)
     real,dimension(ip,0:iq+1,ipa),intent(out) :: Daa0,DaE0,DEE0
     real,dimension(iw) :: keVlog,logDoo
     !!real,dimension(iw) :: b,c,d
     real minDoo,sign1,DaEmax,dE,dE2
     integer j,k,m,kk

     keVlog(:)=log(keV(:))
     Daa0(:,:,:)=0.   
     DaE0(:,:,:)=0.   
     DEE0(:,:,:)=0.   
     do j=1,ip
        do m=1,ipa
!(1) interpolate Daa
           minDoo=minval(Daa(j,:,m),Daa(j,:,m).gt.0.)
           minDoo=log(minDoo)
           logDoo(:)=-115.
           where(Daa(j,:,m).gt.0.) logDoo(:)=log(Daa(j,:,m))
           do k=0,iq+1
           ! linear interpolation
              call lintp(keVlog,logDoo,iw,logE(k),Daa0(j,k,m))
              if (Daa0(j,k,m).ge.minDoo) then
                 Daa0(j,k,m)=exp(Daa0(j,k,m))
              else
                 Daa0(j,k,m)=0.
              endif
           enddo
!(2) interpolate DEE
           minDoo=minval(DEE(j,:,m),DEE(j,:,m).gt.0.)
           minDoo=log(minDoo)
           logDoo(:)=-115.
           where(DEE(j,:,m).gt.0.) logDoo(:)=log(DEE(j,:,m))
           do k=0,iq+1
           ! linear interpolation
              call lintp(keVlog,logDoo,iw,logE(k),DEE0(j,k,m))
              if (DEE0(j,k,m).ge.minDoo) then
                 DEE0(j,k,m)=exp(DEE0(j,k,m))
              else
                 DEE0(j,k,m)=0.
              endif
           enddo
!(3) interpolate DaE
           minDoo=minval(abs(DaE(j,:,m)),abs(DaE(j,:,m)).gt.0.)
           minDoo=log(minDoo)
           logDoo(:)=-115.
           where(abs(DaE(j,:,m)).gt.0.) logDoo(:)=log(abs(DaE(j,:,m)))
           do k=0,iq+1
           ! linear interpolation
              call lintp(keVlog,logDoo,iw,logE(k),DaE0(j,k,m))
              call lintp(keVlog,xSign(j,:,m),iw,logE(k),sign1)
              if (DaE0(j,k,m).ge.minDoo) then
                 DaE0(j,k,m)=sign(exp(DaE0(j,k,m)),sign1)
              else
                 DaE0(j,k,m)=0.
              endif
           ! DaE < sqrt(Daa * DEE)
              if (Daa0(j,k,m).eq.0..or.DEE0(j,k,m).eq.0.) then
                 DaE0(j,k,m)=0.
              else
                 DaEmax=sqrt(Daa0(j,k,m)*DEE0(j,k,m))
                 if (abs(DaE0(j,k,m)).gt.DaEmax) &
                    DaE0(j,k,m)=sign(DaEmax,sign1)
              endif
           enddo
        enddo
     enddo
     end subroutine interpol_D_coef


!*****************************************************************************
!                         interpol_D_coefK
!  Routine interpolates Dqq1,Dqq2 from a0=1,...,89 to K grids
!      and calculates jacobian dQ2/dE
!*****************************************************************************
  subroutine interpol_D_coefK
  use constants, only: pi
  use cimigrid_dim, only: ir,ip,ik,ipc,iph,ipa
  use cread2, only: ichor,ihiss
  use cfield, only: y,iba,ro
  use waveDiffCoef, only: cPA,cOmpe
  use cWpower, only: rppa
  !!use diagDiffCoef, only: iq,cDqq1,cDqq2,hDqq1,hDqq2,Dqq1K,&
  !!                        Dqq2K,E_Q2c,E_Q2h,dQ2dEK,E_Q2K,VarQ,&
  !!                        ipc0,iph0,iLpp
  implicit none

  real a0(0:ik+1),dQ(iq)
  integer i,j,k,m,l,m1,m2,ipc1,iph1
 
  call find_ompe_index
  do k=1,iq
     dQ(k)=VarQ(k+1)-VarQ(k-1)
  enddo

  do j=1,ip
     call locate1(ro(1:iba(j),j),iba(j),rppa(j),iLpp(j))
!(1) interpolate begins when L>Lpp for chorus
     do i=iLpp(j)+1,iba(j)
        if (ichor.eq.1) then
           ipc1=ipc0(i,j)
           a0(0:ik+1)=asin(y(i,j,0:ik+1))*180./pi
           m1=ik+1
           m2=0
           do m=0,ik+1
              if (a0(m).ge.cPA(1).and.a0(m).le.cPA(ipa)) then
                 if (m.le.m1) m1=m
                 if (m.ge.m2) m2=m
!(2) interpolate Dqq1, Dqq2 to y grids
                 do k=1,iq
                    call lintp(cPA,cDqq1(ipc1,k,:),ipa,a0(m),Dqq1K(i,j,k,m))
                    call lintp(cPA,cDqq2(ipc1,k,:),ipa,a0(m),Dqq2K(i,j,k,m))
                    call lintp(cPA,E_Q2c(ipc1,k,:),ipa,a0(m),E_Q2K(i,j,k,m))
                 enddo   ! end of k
              endif
           enddo         ! end of m
!(3) assign Dqq1, Dqq2, where <1 deg and >89 deg
           do m=0,m1
              Dqq1K(i,j,:,m)=Dqq1K(i,j,:,m1)
              Dqq2K(i,j,:,m)=Dqq2K(i,j,:,m1)
              E_Q2K(i,j,:,m)=VarQ(1:iq)
           enddo         ! end of m
           do m=m2,ik+1
              Dqq1K(i,j,:,m)=Dqq1K(i,j,:,m2)
              Dqq2K(i,j,:,m)=Dqq2K(i,j,:,m2)
              E_Q2K(i,j,:,m)=E_Q2K(i,j,:,m2)
           enddo         ! end of m
!(4) Calculate jacobian dQ2/dE
           do m=0,ik+1
              do k=2,iq-1
                 dQ2dEK(i,j,k,m)=dQ(k)/(E_Q2K(i,j,k+1,m)-E_Q2K(i,j,k-1,m))
              enddo      ! end of k
              dQ2dEK(i,j,1,0:ik+1)=dQ2dEK(i,j,2,0:ik+1)
              dQ2dEK(i,j,iq,0:ik+1)=dQ2dEK(i,j,iq-1,0:ik+1)
           enddo         ! end of m
        endif
     enddo               ! end of i
  enddo                  ! end of j
  do j=1,ip
!(5) interpolate begins when L>Lpp for hiss
     do i=i,iLpp(j)
        if (ihiss.eq.1) then
           iph1=iph0(i,j)
           a0(:)=asin(y(i,j,:))*180./pi
           m1=ik+1
           m2=0
           do m=0,ik+1
              if (a0(m).ge.cPA(1).and.a0(m).le.cPA(ipa)) then
                 if (m.le.m1) m1=m
                 if (m.ge.m2) m2=m
!(6) interpolate Dqq1, Dqq2 to y grids
                 do k=1,iq
                    call lintp(cPA,hDqq1(iph1,k,:),ipa,a0(m),Dqq1K(i,j,k,m))
                    call lintp(cPA,hDqq2(iph1,k,:),ipa,a0(m),Dqq2K(i,j,k,m))
                    call lintp(cPA,E_Q2h(ipc1,k,:),ipa,a0(m),E_Q2K(i,j,k,m))
                 enddo   ! end of k
              endif
           enddo         ! end of m
!(7) assign Dqq1, Dqq2, where <1 deg and >89 deg
           do m=0,m1
              Dqq1K(i,j,:,m)=Dqq1K(i,j,:,m1)
              Dqq2K(i,j,:,m)=Dqq2K(i,j,:,m1)
              E_Q2K(i,j,:,m)=VarQ(1:iq)
           enddo         ! end of m
           do m=m2,ik+1
              Dqq1K(i,j,:,m)=Dqq1K(i,j,:,m2)
              Dqq2K(i,j,:,m)=Dqq2K(i,j,:,m2)
              E_Q2K(i,j,:,m)=E_Q2K(i,j,:,m2)
           enddo         ! end of m
!(8) Calculate jacobian dQ2/dE
           do m=0,ik+1
              do k=2,iq-1
                 dQ2dEK(i,j,k,m)=dQ(k)/(E_Q2K(i,j,k+1,m)-E_Q2K(i,j,k-1,m))
              enddo   ! end of k
              dQ2dEK(i,j,1,0:ik+1)=dQ2dEK(i,j,2,0:ik+1)
              dQ2dEK(i,j,iq,0:ik+1)=dQ2dEK(i,j,iq-1,0:ik+1)
           enddo         ! end of m
        endif
     enddo               ! end of i
  enddo                  ! end of j
  
  end subroutine interpol_D_coefK


!*****************************************************************************
!                           mapPSDtoQ 
!  Routine maps PSD to fixed grids of (K,Q) from (K,V)
!*****************************************************************************
  subroutine mapPSDtoQ
  use cimigrid_dim, only: ir,ip,iw,ik,ns
  use cread2, only: ijs,js,ichor,ihiss
  use cgrid, only: xjac
  use cfield, only: iba,y,ekev
  use cinitial, only: f2
  !!use diagDiffCoef, only: iq,E_Q2K,PSD_Q,iLpp
  implicit none

  real logEQ,logPSD1,EQ,dE,dE2
  real,dimension(iw) :: ekev0,psd0,logPSD,logEkeV,b,c,d
  integer n,i,j,k,m,nel,k1,k2,kk
  logical DoMapping

!(1) Find the "n" corresponds to electrons, nel 
  nel=0
  do n=1,ijs
     if (js(n).eq.4) nel=n
  enddo

  do j=1,ip
     do i=1,iba(j)
        DoMapping=.False.
        if (ichor.eq.1.and.i.gt.iLpp(j)) then
           DoMapping=.True.
        else if (ihiss.eq.1.and.i.le.iLpp(j)) then
           DoMapping=.True.
        endif
        if (DoMapping) then
           do m=1,ik
              ekev0(:)=ekev(nel,i,j,1:iw,m)
              logEkeV(:)=log(ekev0)
              psd0(:)=f2(nel,i,j,1:iw,m)/xjac(nel,i,1:iw,m)
              logPSD(:)=log(psd0)
              where(psd0(:).lt.1.e-50) logPSD(:)=-50.
!(1) find cubic spline interpolation coefficients.
              call spline (logEkeV, logPSD, b, c, d, iw)
              call locate1(E_Q2K(i,j,:,m),iq,ekev0(1 ),k1)
              call locate1(E_Q2K(i,j,:,m),iq,ekev0(iw),k2)
              k1=k1+1
              do k=k1,k2
!(2) map PSD to fixed Q2 grids
                 EQ=E_Q2K(i,j,k,m)
                 logEQ=log(EQ)
                 call locate1(logEkeV,iw,logEQ,kk)
                 dE=logEQ-logEkeV(kk)
                 dE2=dE*dE
                 logPSD1=logPSD(kk)+b(kk)*dE+c(kk)*dE2+d(kk)*dE2*dE
                 if (logPSD1.lt.min(logPSD(kk),logPSD(kk+1)).or.&
                     logPSD1.gt.max(logPSD(kk),logPSD(kk+1))) &
                    call lintp(logEkeV,logPSD,iw,logEQ,logPSD1)
                 PSD_Q(i,j,k,m)=exp(logPSD1)
              enddo      ! end of k
!(3) fill PSD beyond interpolation ranges
              PSD_Q(i,j,1:k1,m)=PSD_Q(i,j,k1,m)
              PSD_Q(i,j,k2:iq,m)=PSD_Q(i,j,k2,m)
           enddo         ! end of m
        endif            ! end of DoMapping
     enddo               ! end of i
  enddo                  ! end of j

  end subroutine mapPSDtoQ


!*****************************************************************************
!                           mapPSDtoE 
!  Routine maps PSD to fixed grids of (K,V) from (K,Q)
!*****************************************************************************
  subroutine mapPSDtoE
  use cimigrid_dim, only: ir,ip,iw,ik,ns
  use cread2, only: ijs,js,ichor,ihiss
  use cgrid, only: xjac
  use cfield, only: iba,y,ekev
  use cinitial, only: f2
  !!use diagDiffCoef, only: iq,E_Q2K,PSD_Q,iLpp
  implicit none

  real logEkeV,logPSD1,ekev0,dE,dE2
  real,dimension(iq) :: logPSD,EQ,logEQ,b,c,d
  integer n,i,j,k,m,nel,kk
  logical DoMapping

!(1) Find the "n" corresponds to electrons, nel 
  nel=0
  do n=1,ijs
     if (js(n).eq.4) nel=n
  enddo

  do j=1,ip
     do i=1,iba(j)
        DoMapping=.False.
        if (ichor.eq.1.and.i.gt.iLpp(j)) then
           DoMapping=.True.
        else if (ihiss.eq.1.and.i.le.iLpp(j)) then
           DoMapping=.True.
        endif
        if (DoMapping) then
           do m=1,ik
!(1) map PSD to fixed E grids for chorus L>Lpp
              EQ(:)=E_Q2K(i,j,1:iq,m)
              logEQ(:)=log(EQ)
              logPSD(:)=log(PSD_Q(i,j,1:iq,m))
              where(logPSD(:).lt.-50.) logPSD(:)=-50.
              call spline (logEQ, logPSD, b, c, d, iq)
              do k=2,iw    ! keep PSD(k=1) unchanged during mapping
                 ekev0=ekev(nel,i,j,k,m)
                 if (ekev0.ge.EQ(1).and.ekeV0.le.EQ(iq)) then
                    logEkeV=log(ekev0)
                    call locate1(EQ,iq,ekev0,kk)
                    dE=logEkeV-logEQ(kk)
                    dE2=dE*dE
                    logPSD1=logPSD(kk)+b(kk)*dE+c(kk)*dE2+d(kk)*dE2*dE
                    if (logPSD1.lt.min(logPSD(kk),logPSD(kk+1)).or.&
                        logPSD1.gt.max(logPSD(kk),logPSD(kk+1))) &
                       call lintp(logEq,logPSD,iq,logEkeV,logPSD1)
                    f2(nel,i,j,k,m)=exp(logPSD1)*xjac(nel,i,k,m)
                 endif
              enddo      ! end of k
           enddo         ! end of m
        endif            ! end of DoMapping
     enddo            ! end of i
  enddo               ! end of j

  end subroutine mapPSDtoE


!*****************************************************************************
!                           runge_kuntta    
!  Routine calculates y(t) integrating dy/dt and dt using Runge-Kuntta method
!  
!  y(t+dt) = y(t) + (k1 + 2*k2 + 2*k3 + k4) * dt/6
!  k1=dy/dt(t, y(t))
!  k2=dy/dt(t + 0.5*dt, y(t) + 0.5*dt*k1)
!  k3=dy/dt(t + 0.5*dt, y(t) + 0.5*dt*k2)
!  k4=dy/dt(t + dt, y(t) + dt*k3)
!*****************************************************************************
   subroutine runge_kutta(y0,dt,y,t,dydt,ny,nt,y_out)
   implicit none
   integer,intent(in) :: ny,nt
   real,intent(in) :: y0,dydt(ny,nt),dt(nt),y(ny),t(nt)
   real,intent(out) :: y_out(nt)
   integer i,i_1
   real y1,dt1,dt2,k1,k2,k3,k4

   y_out(1)=y0
   do i=2,nt
      i_1=i-1
      y1=y_out(i_1)
      dt1=dt(i_1)
      dt2=0.5*dt1
      call lintp2(y,t,dydt,ny,nt,y1,t(i_1),k1)
      call lintp2(y,t,dydt,ny,nt,y1+dt2*k1,t(i_1)+dt2,k2)
      call lintp2(y,t,dydt,ny,nt,y1+dt2*k2,t(i_1)+dt2,k3)
      call lintp2(y,t,dydt,ny,nt,y1+dt1*k3,t(i_1)+dt1,k4)
      y_out(i)=y1+dt1*(k1 + 2.*k2 + 2*k3 + k4)/6.
   enddo

   end subroutine runge_kutta


!*****************************************************************************
!                        write_Q2info
!  Routine prints Q2 grid information.      
!   (1) constant Q1 in (a0,E), setting Q2=E
!   (2) constant Q2 in (a0,E), setting Q1=a0
!   (3) Daa, Dq2,q2
!*****************************************************************************
  subroutine write_Q2info(DEE,DaE,Daa0,DaE0,DEE0,Dqq1,Dqq2,keV,Eq,E,ip,iw,iq)
  use waveDiffCoef, only: ipa,cPA,cOmpe,hOmpe
  !!use diagDiffCoef, only: VarQ
  implicit none

  integer,intent(in) :: ip,iw,iq
  real,intent(in),dimension(ip,iw,ipa) :: DEE,DaE
  real,intent(in),dimension(ip,0:iq+1,ipa) :: Daa0,DaE0,DEE0,E
  real,intent(in),dimension(ip,iq,ipa) :: Dqq1,Dqq2
  real,intent(in) :: keV(iw),Eq(0:iq+1)
  real a(iw),DaEEE(iw,ipa),dE
  integer j,k,m
  
!!!! write Dqq
  open(unit=48,file='D_LBchorus_QQ.dat')
  write(48,'(f7.0,f7.2)') 10000.,6.5 
  write(48,*) ip,iq,ipa
  write(48,'(10f7.2)') cOmpe
  write(48,'(8f12.5)') VarQ(1:iq)  
  do j=1,ip
     do k=1,iq
        write(48,'(f7.2,f12.5,a)') cOmpe(j),VarQ(k),'     ! fpe/fce   E(keV)'
        write(48,*) 'alpha0    Daa/p^2     DQQ/Q^2  (sec^-1)'
        do m=1,ipa
           write(48,'(f5.1,1p2E18.5)') cPA(m),Dqq1(j,k,m),Dqq2(j,k,m)
        enddo
     enddo
  enddo
  close(48)
  write(*,*) 'write Dq1q1, Dq2q2'
     
  open(unit=60,file='constQ.dat')
  write(60,*) ip,iw,ipa
  write(60,'(1p10E11.3)') keV(:)
  do j=1,ip
     DaEEE(:,:)=DaE(j,:,:)/DEE(j,:,:)
     where (DEE(j,:,:).eq.0.) DaEEE(:,:)=0.
     do m=1,ipa
        a(1)=cPA(m)
        do k=2,iw
           dE=(keV(k)-keV(k-1))
           a(k)=a(k-1)+0.5*(DaEEE(k,m)/keV(k)+DaEEE(k-1,m)/keV(k-1))*dE
        enddo
        write(60,'(1p10E12.4)') a(:)
     enddo
  enddo
  close(60)
  open(unit=60,file='constQ2.dat')
  write(60,*) ip,iq,ipa
  write(60,'(1p10E11.3)') Eq(1:iq)
  do j=1,ip
     do k=1,iq
        write(60,'(1p10E12.4)') E(j,k,:)
     enddo
  enddo
  close(60)

  end subroutine write_Q2info


! ************************************************************************
!                        calc_num_ptl          
!
! Subroutine calculates a value proportional to the number of particles
! ************************************************************************
  subroutine calc_num_ptl(f2,x,n)
  implicit none

  integer,intent(in) :: n
  real,intent(in) :: f2(n),x(n)
  real ptl
  integer i

  ptl=0.
  do i=2,n
     ptl=ptl+0.5*(f2(n)+f2(n-1))*(x(n)-x(n-1))
  enddo
  ptl=ptl+f2(1)*x(1)
  write(*,'(a,1pE11.3)') '# particle =',ptl

  end subroutine calc_num_ptl
  
  end module
