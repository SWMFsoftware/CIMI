!***********************************************************************
!                              calc_Lstar1  
!
! Routine modified from calc_Lstar. This routine calculates L * at given 
! (xlati,xmlt) in a given magnetic field configuration for the CIMI model 
! of Mei-Ching H. Fok at (Code 673) NASA/GSFC.
! Calculation of magnetic flux and L* is made at ionophere
!
! 
! In put: xlati,xmlt,bo,iba (field parameters)
! Output: L*, max L*
!
! Creaded by Suk-Bin Kang (Code 673) at NASA/GSFC on November 12, 2015
!
! Subroutines
! - calc_lstar1: calculates L* for K = 0 (or 90deg) particles.
! - calc_lstar2: calculates L* with K dependence.
! - locateB: is modified from locate1 and calculates latitudinal index
!            of bo or bm excluding magnetic island.
! - dBdMLAT: calculates magnetic flux at the ionosphere
!
! Parameters
! - BE: magnetic field strength at the Earth surface at magnetic equator.
! - Bflux: magnetic flux enclosed by drift path at each grid.
! - dBdlat: analytic integration of B * dMLAT
! - Lstar: out L*
! - Lstar_max: L* of last closed drift shell.
! - dmlt: MLT grid size at ionosphere.
! - mlt1,mlt2,dp,dmlt1,j1,j2: variables to be used in dmlt.
! - xlat_i, dBdlat_i: interpolated MLAT and dBdlat, respectively.
! - b1,b2,mlat2: variables to be used in xlat_i and dBdlat_i.
! - logBo_i: bo from i=1 to iba(j).
! - jdawn,jdusk,jnoon: MLT indices for dawn,dusk, and noon, respectively.
! - B_island: contains information if bo(i,j) is magnetic island
!           if 0, then normal, if 1, then magnetic island.
! - dBdMLAT: function that calculates integral B * dMLAT
!
! Modification history
!  * November 16, 2015 - correct interpolation of xlat_i between ii and ii+1
!                        to ii and ii1, skipping magnetic island.
!  * November 19, 2015 - To avoid magnetic island
!                        locateB: ju=j(ii)+1 -> j(ii)-1
!                                 jl=j(1)    -> j(1)+1
!  * November 20. 215  - direct linear interpolation -> log linear interpol.
!  * February 11, 2016 - allocable logBo_i(:) -> real logBo_i(ir)
!                      - remove allocate and deallocate logBo_i
!  * March 21, 2016    - Add calc_lstar2.f90
!  * May 22, 2018      - Added verbose 
!**************************************************************************
Module ModLstar

  use ModCimiGrid,	ONLY: np, nt, nk, nLonPar, MinLonPar, MaxLonPar, &
       MinLonGlobal, MaxLonGlobal, iProc, nProc, iComm
  use ModMpi
  
  logical	:: &
       DoVerboseLstar = .false.
  
  real, allocatable :: Lstar_C(:,:), Lstar_max, Lstarm(:,:,:), Lstarm_max(:)

  ! Full-longitude L*, only filled on iProc 0 by gather_Lstar
  real, allocatable :: LstarGlobal_C(:,:), LstarmGlobal(:,:,:)

  ! Fields gathered to the full longitude range 1:nt on every processor
  real,    allocatable :: boAll(:,:), bmAll(:,:,:)
  integer, allocatable :: ibaAll(:)

  ! Count and 0-based global start of every processor's longitude block,
  ! only used as MPI_allgatherv counts/displacements
  integer, allocatable :: nLonProc(:), iLonOffProc(:)
  
  public :: calc_Lstar1
  public :: calc_Lstar2
  public :: gather_Lstar
  
contains
  
!**************************************************************************
  subroutine gather_Lstar_field
!
! Gathers the local bo, bm and iba to boAll(np,nt), bmAll(np,nt,nk) and
! ibaAll(nt) on every processor. bo and bm are packed into one buffer so
! that a single MPI_allgatherv moves both.
!**************************************************************************
  use ModCimiTrace, ONLY: bo, bm, iba
  implicit none

  integer :: iLonStart, m, iError
  integer :: iCount(0:nProc-1), iDisp(0:nProc-1)
  ! Slot 0 holds bo, slots 1:nk hold bm, so every processor's longitude
  ! block is contiguous in memory
  real    :: BufferSend(np,0:nk,nLonPar), BufferRecv(np,0:nk,nt)

  ! Size and 0-based global start of every processor's longitude block
  if (.not.allocated(nLonProc)) then
     allocate(nLonProc(0:nProc-1), iLonOffProc(0:nProc-1))
     iLonStart = MinLonGlobal - 1
     call MPI_allgather(nLonPar, 1, MPI_INTEGER, nLonProc, 1, MPI_INTEGER, &
          iComm, iError)
     call MPI_allgather(iLonStart, 1, MPI_INTEGER, iLonOffProc, 1, MPI_INTEGER, &
          iComm, iError)
     if (sum(nLonProc).ne.nt) then
        if (iProc.eq.0) write(*,*) 'ModLstar: local longitude counts sum to', &
             sum(nLonProc), ' but nt =', nt
        call MPI_abort(iComm, 1, iError)
     endif
  endif

  if (.not.allocated(boAll))  allocate(boAll(np,nt))
  if (.not.allocated(bmAll))  allocate(bmAll(np,nt,nk))
  if (.not.allocated(ibaAll)) allocate(ibaAll(nt))

  BufferSend(:,0,:) = bo(1:np,MinLonPar:MaxLonPar)
  do m = 1, nk
     BufferSend(:,m,:) = bm(1:np,MinLonPar:MaxLonPar,m)
  enddo

  iCount = np*(nk+1)*nLonProc
  iDisp  = np*(nk+1)*iLonOffProc
  call MPI_allgatherv(BufferSend, np*(nk+1)*nLonPar, MPI_REAL, &
       BufferRecv, iCount, iDisp, MPI_REAL, iComm, iError)

  boAll = BufferRecv(:,0,:)
  do m = 1, nk
     bmAll(:,:,m) = BufferRecv(:,m,:)
  enddo

  call MPI_allgatherv(iba(MinLonPar:MaxLonPar), nLonPar, MPI_INTEGER, &
       ibaAll, nLonProc, iLonOffProc, MPI_INTEGER, iComm, iError)

  end subroutine gather_Lstar_field

!**************************************************************************
  subroutine gather_Lstar
!
! Gathers the local Lstar_C and Lstarm to LstarGlobal_C and LstarmGlobal
! on iProc 0 with a single MPI_gatherv (same packing as gather_Lstar_field)
!**************************************************************************
  implicit none

  integer :: m, iError
  integer :: iCount(0:nProc-1), iDisp(0:nProc-1)
  real    :: BufferSend(np,0:nk,nLonPar), BufferRecv(np,0:nk,nt)

  BufferSend(:,0,:) = Lstar_C(1:np,MinLonPar:MaxLonPar)
  do m = 1, nk
     BufferSend(:,m,:) = Lstarm(1:np,MinLonPar:MaxLonPar,m)
  enddo

  iCount = np*(nk+1)*nLonProc
  iDisp  = np*(nk+1)*iLonOffProc
  call MPI_gatherv(BufferSend, np*(nk+1)*nLonPar, MPI_REAL, &
       BufferRecv, iCount, iDisp, MPI_REAL, 0, iComm, iError)

  if (iProc.eq.0) then
     LstarGlobal_C = BufferRecv(:,0,:)
     do m = 1, nk
        LstarmGlobal(:,:,m) = BufferRecv(:,m,:)
     enddo
  endif

  end subroutine gather_Lstar

!**************************************************************************
  subroutine calc_Lstar1
    
  use ModNumConst, 	ONLY: pi => cPi
  use ModCimiPlanet,	ONLY: re_m, xme => dipmom, rc
  use ModCimiTrace,	ONLY: ro, rb
  use ModCimiGrid,	ONLY: xlati => xlatr, ir => np, ip => nt, xmlt
  implicit none

  ! Bflux is only needed for the longitudes owned by this processor
  real BE, Bflux(ir,MinLonPar:MaxLonPar), dBdlat(ir), dmlt(ip), mlt1, mlt2, &
       dp, dmlt1, xlat_i, dBdlat_i, b1, b2, mlat2, logBo, logBo_i(ir), &
       LstarMaxLoc
  ! j0 = local longitude index, jg = its global index, jl = local index of
  ! a given global longitude
  integer i,j,j0,jg,jl,ib0,ii,ii1,jdawn,jdusk,jnoon,j1,j2,B_island(ir,ip),iError
  !real,external :: dBdMLAT

  ! Gather the drift-path data (all longitudes) onto every processor

  if(.not. allocated(Lstar_C)) allocate(Lstar_C(np,MinLonPar:MaxLonPar), Lstar_max, &
          Lstarm(np,MinLonPar:MaxLonPar,nk), Lstarm_max(nk), LstarGlobal_C(np,nt), LstarmGlobal(np,nt,nk))

  call gather_Lstar_field

  jdawn=nint(float(ip)*0.75)+1     ! MLT index at dawn (global)
  jdusk=nint(float(ip)*0.25)+1     ! MLT index at dusk (global)
  jnoon=1      ! MLT index at noon (global)
  dp=24./ip                        ! MLT grid size
  ib0=11          ! ro(ib0,*)=2.1, inner boundary of L* calculation

  BE=xme/re_m**3                    ! magnetic field at L* = 1 in


! Calulates dmlt at MLT grid 
  dmlt(1:ip)=24./ip

! Find magnetic island
! Needed for ALL longitudes on every processor, since the flux integral
! below loops over all j and looks at B_island(ii1,j). The cost is only
! O(ir*ip), negligible compared to the flux integral, so it is computed
! redundantly instead of communicated.
  B_island(1:ir,1:ip)=0    ! 0=normal, 1=island
  do j=1,ip
     ii=1000
     do i=ibaAll(j)-1,1,-1
        if (ii.lt.i) cycle
        if (boAll(i,j).lt.boAll(i+1,j)) then
           B_island(i,j)=1          
           ii=i
           do while ( boAll(ii-1,j).lt.boAll(i+1,j) .and. ii.gt.1)
             ii=ii-1
             B_island(ii,j)=1
           enddo
           ! print once, from the processor that owns longitude j
           if ( DoVerboseLstar .and. j.ge.MinLonGlobal .and. j.le.MaxLonGlobal ) &
                write(*,'("IM0: bo(i,j) has a island at, j = ",i3,"  i = ",i3," -",i3 )') j,ii,i
        endif
     enddo
  enddo
   
! Caculates L* at each grid, only for the longitudes owned by this processor
  Bflux=0.  ! Initialize magnetic fluxes
  Lstar_C(1:ib0-1,MinLonPar:MaxLonPar)=ro(1:ib0-1,MinLonPar:MaxLonPar)  ! set L*=ro at i<ib0
  do j0=MinLonPar,MaxLonPar      ! local MLT grid
     jg = MinLonGlobal + j0 - 1            ! global MLT index of j0
     do i=ib0,ibaAll(jg)
        Bflux(i,j0)=0.             ! initialization of magnetic flux  
        if (B_island(i,jg).eq.1) then 
           lstar_C(i,j0)=-1.*rb      
           cycle                   ! skipping magnetic island
        endif
        logBo=log10(boAll(i,jg))
        do j=1,ip                  ! drift path spans ALL (global) longitudes
           logBo_i(1:ibaAll(j))=log10(boAll(1:ibaAll(j),j))
           ! Find the MLAT index for the same bo
           IF ( ( j .NE. jg ) .AND. &
                ( logBo .LT. logBo_i( 1 ) ) .AND. &
                ( logBo .GE. logBo_i( ibaAll( j ) ) ) ) &
                CALL locateB( logBo_i( 1 : ibaAll( j ) ), ibaAll( j ), logBo, ii )
           if (logBo.ge.logBo_i(1)) ii=1
           if (logBo.lt.logBo_i(ibaAll(j))) ii=ibaAll(j)
           if (j.eq.jg) ii=i
           if (ii.lt.ibaAll(j)) then
              ii1=ii+1
              do while (ii1.lt.ibaAll(j) .and. B_island(ii1,j).eq.1 ) 
                 ii1=ii1+1
              enddo
              b2=logBo_i(ii1)-logBo
              b1=logBo-logBo_i(ii)
              xlat_i=(xlati(ii)*b2+xlati(ii1)*b1)/(b2+b1) ! interpolates xlati
           else
              xlat_i=xlati(ibaAll(j))                   
           endif
           dBdlat_i=dBdMLAT(BE,rc,xlat_i)
           Bflux(i,j0)=Bflux(i,j0)+dBdlat_i*dmlt(j)*pi/12. ! analytic intergration of magnetic flux
        enddo          ! end of j
        lstar_C(i,j0)=2.*pi*BE/Bflux(i,j0)
     enddo             ! end of i
     if (ibaAll(jg).lt.ir) lstar_C(ibaAll(jg)+1:ir,j0)=rb     ! arbitrary big L* 
  enddo                ! end of j0

  ! Determines L* max
  ! The noon value and the dawn-dusk values live on different processors, so
  ! each processor takes the min over the points it owns and the result is
  ! combined with an allreduce. Same bounds as the serial code: noon value,
  ! then j = jdawn,jdusk.
  LstarMaxLoc = huge(LstarMaxLoc)
  if (jnoon.ge.MinLonGlobal .and. jnoon.le.MaxLonGlobal) then
     jl = jnoon - MinLonGlobal + 1                    ! local index of global jnoon
     LstarMaxLoc = Lstar_C(ibaAll(jnoon),jl)          ! L* at noon magnetopause
  endif
  do j=max(jdawn,MinLonGlobal),min(jdusk,MaxLonGlobal)   ! find L* max from dawn to dusk
     jl = j - MinLonGlobal + 1                        ! local index of global j
     if (LstarMaxLoc.gt.Lstar_C(ibaAll(j),jl)) LstarMaxLoc = Lstar_C(ibaAll(j),jl)
  enddo
  call MPI_allreduce(LstarMaxLoc, Lstar_max, 1, MPI_REAL, MPI_MIN, iComm, iError)

  if ( DoVerboseLstar ) then
     
     if (iProc.eq.0) write(*,'("IM0: L*max =",f8.2)') Lstar_max
     ! print once, from the processor that owns global longitude 1
     if (1.ge.MinLonGlobal .and. 1.le.MaxLonGlobal) then
        jl = 1 - MinLonGlobal + 1
        write(*,'("IM0: L* =",10f8.2)') (Lstar_C(i,jl),i=1,ir)
     endif

  endif

  end subroutine calc_Lstar1

!**************************************************************************
  subroutine calc_Lstar2
!
! Routine calcuates L* with K (second adiabatic invariant) dependence.
!**************************************************************************

  use ModNumConst,	ONLY: pi=>cPi
  use ModCimiPlanet,	ONLY: re_m, xme => dipmom, rc
  use ModCimiTrace,	ONLY: ro, rb
  use ModCimiGrid, 	ONLY: xlati=>xlatr, ir=>np, ip=>nt, xmlt, ik=>nk
  implicit none

  real BE, Bflux(ir,MinLonPar:MaxLonPar), dBdlat(ir), dmlt(ip), mlt1, mlt2, &
       dp, dmlt1, xlat_i, dBdlat_i, b1, b2, mlat2, logBm, logbm_i(ir), &
       LstarMaxLoc(ik)
  integer i,j,m,j0,jg,jl,ib0,ii,ii1,jdawn,jdusk,jnoon,j1,j2,B_island(ir,ip),iError
 
  ! Gather the drift-path data (all longitudes) onto every processor
  call gather_Lstar_field

  jdawn=nint(float(ip)*0.75)+1     ! MLT index at dawn (global)
  jdusk=nint(float(ip)*0.25)+1     ! MLT index at dusk (global)
  jnoon=1      ! MLT index at noon (global)
  dp=24./ip                        ! MLT grid size
  ib0=11          ! ro(ib0,*)=2.1, inner boundary of L* calculation

  BE=xme/re_m**3                    ! mangetic field at L* = 1 in

! Calulates dmlt at MLT grid 
  dmlt(1:ip)=24./ip

  do m=1,ik
! Find magnetic island (all longitudes, computed redundantly on each
! processor, see calc_Lstar1)
     B_island(1:ir,1:ip)=0    ! 0=normal, 1=island
     do j=1,ip
        ii=1000
        do i=ibaAll(j)-1,1,-1
           if (ii.lt.i) cycle
           if (bmAll(i,j,m).lt.bmAll(i+1,j,m)) then
              B_island(i,j)=1          
              ii=i
              do while (bmAll(ii-1,j,m).lt.bmAll(i+1,j,m) .and. ii.gt.1)
                ii=ii-1
                B_island(ii,j)=1
              enddo
              if ( DoVerboseLstar .and. j.ge.MinLonGlobal .and. j.le.MaxLonGlobal ) &
                   write(*,'("IM0: bm(i,j,m) has a island at ,m,j = ",2i3,"  i = ",i3," -",i3 )') m,j,ii,i
           endif
        enddo
     enddo
   
! Caculates L* at each grid, only for the longitudes owned by this processor
     Bflux=0.  ! Initialize magnetic fluxes
     Lstarm(1:ib0-1,MinLonPar:MaxLonPar,m)=ro(1:ib0-1,MinLonPar:MaxLonPar)  ! set L*=ro at i<ib0
     do j0=MinLonPar,MaxLonPar      ! local MLT grid
        jg = MinLonGlobal + j0 - 1            ! global MLT index of j0
        do i=ib0,ibaAll(jg)
           Bflux(i,j0)=0.             ! initialization of magnetic flux  
           if (B_island(i,jg).eq.1) then 
              lstarm(i,j0,m)=-1.*rb      
              cycle                   ! skipping magnetic island
           endif
           logBm=log10(bmAll(i,jg,m))
           do j=1,ip                  ! drift path spans ALL (global) longitudes
              logbm_i(1:ibaAll(j))=log10(bmAll(1:ibaAll(j),j,m))            ! Drift path is estimated with bm 
              if (j.ne.jg.and.logBm.lt.logbm_i(1).and.logBm.ge.logbm_i(ibaAll(j))) &
              call locateB(logbm_i(1:ibaAll(j)),ibaAll(j),logBm,ii) ! Find the MLAT index for the same bo
              if (logBm.ge.logbm_i(1)) ii=1
              if (logBm.lt.logbm_i(ibaAll(j))) ii=ibaAll(j)
              if (j.eq.jg) ii=i
              if (ii.lt.ibaAll(j)) then
                 ii1=ii+1
                 do while (B_island(ii1,j).eq.1)
                    ii1=ii1+1
                 enddo
                 b2=logbm_i(ii1)-logBm
                 b1=logBm-logbm_i(ii)
                 xlat_i=(xlati(ii)*b2+xlati(ii1)*b1)/(b2+b1) ! interpolates xlati
              else
                 xlat_i=xlati(ibaAll(j))                   
              endif
              dBdlat_i=dBdMLAT(BE,rc,xlat_i)
              Bflux(i,j0)=Bflux(i,j0)+dBdlat_i*dmlt(j)*pi/12. ! analytic intergration of magnetic flux
           enddo ! end of j
           lstarm(i,j0,m)=2.*pi*BE/Bflux(i,j0)
        enddo    ! end of i
        if (ibaAll(jg).lt.ir) lstarm(ibaAll(jg)+1:ir,j0,m)=rb     ! arbitrary big L* 
     enddo       ! end of j0

  ! Determines L* max (local part; combined across processors after the m loop)
     LstarMaxLoc(m) = huge(LstarMaxLoc)
     if (jnoon.ge.MinLonGlobal .and. jnoon.le.MaxLonGlobal) then
        jl = jnoon - MinLonGlobal + 1                    ! local index of global jnoon
        LstarMaxLoc(m) = Lstarm(ibaAll(jnoon),jl,m)      ! L* at noon magnetopause
     endif
     do j=max(jdusk,MinLonGlobal),min(jdawn,MaxLonGlobal)   ! find L* max from dawn to dusk
        jl = j - MinLonGlobal + 1                        ! local index of global j
        if (LstarMaxLoc(m).gt.Lstarm(ibaAll(j),jl,m)) LstarMaxLoc(m)=Lstarm(ibaAll(j),jl,m)
     enddo
  enddo          ! end of k

  ! One reduction for all K values instead of one per m
  call MPI_allreduce(LstarMaxLoc, Lstarm_max, ik, MPI_REAL, MPI_MIN, iComm, iError)

  if ( DoVerboseLstar ) then
     
     if (iProc.eq.0) write(*,'("IM0: L*max =",8f8.2)') (Lstarm_max(m),m=1,ik)
     ! print once, from the processor that owns global longitude 1
     if (1.ge.MinLonGlobal .and. 1.le.MaxLonGlobal) then
        jl = 1 - MinLonGlobal + 1
        write(*,'("IM0: L* =",10f8.2)') (Lstarm(i,jl,1),i=1,ir)
     endif

  endif

end subroutine calc_Lstar2

!**************************************************************************
      subroutine locateB(xx,n,x,j)
!  
!  Routine modified from "locate1" in the CIMI model
!
!  Routine is return a value of j such that x is between xx(j) and xx(j+1),
!   excluding small humps and dips.
!  Routine is still valid if xx is a  monotonically increasing 
!   or decreasing function including small humps or dips by passing those.
!
!**************************************************************************

      implicit none
 
      integer n,i,ii,iii,j,jl,ju,jm,j1,j2
      integer js(n)
      real xx(n),x

!  Make sure xx is increasing or decreasing monotonically
      js(1:n)=1000     ! arbitrary big value 
      ii=0        ! dummy index to count magnetic island
      i=n
      do
         i=i-1
         if (xx(i).le.xx(i+1)) then
            ii=ii+1
            js(ii)=i
            do while (xx(i-1).lt.xx(i+1))
               ii=ii+1
               i=i-1
               js(ii)=i
            enddo
         endif
         if (i.le.1) exit
      enddo

      if (js(1).gt.999) then
         jl=0
         ju=n+1
         do while (ju-jl.gt.1)
            jm=(ju+jl)/2
            if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
               jl=jm
            else
               ju=jm
            endif
         enddo
         j=jl
      endif
      if (js(1).le.n) then
         jl=0
         ju=js(ii)
         do while (ju-jl.gt.1)
            jm=(ju+jl)/2
            if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
               jl=jm
            else
               ju=jm
            endif
         enddo
         j1=jl
         jl=js(1)
         ju=n+1
         do while (ju-jl.gt.1)
            jm=(ju+jl)/2
            if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
               jl=jm
            else
               ju=jm
            endif
         enddo
         j2=jl
         ! avoid a hump or dip
         j=j1
         if (x.le.xx(js(ii)-1).and.x.gt.xx(js(1)+1)) j=js(ii)-1
         if (x.le.xx(js(1)+1)) j=j2
      endif

      end subroutine locateB


!**************************************************************************
!                             dBdMLAT
!  Function calculates interation of dipole Br * dMLAT * r**2 at ionosphere
!  form mlat1 to MLAT = 90 degree
!**************************************************************************
     function dBdMLAT(BE,rc,mlat)

     implicit none
     real BE,rc,mlat,dBdMLAT

     dBdMLAT=BE/rc*(1.-sin(mlat)**2)

     end function dBdMLAT
  end Module ModLstar