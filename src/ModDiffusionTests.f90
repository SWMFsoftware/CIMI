module ModDiffusionTests

    implicit none

    private ! except

    public :: diffusion_test_run

contains
!==============================================================================
subroutine diffusion_test_run(DtAdvance)
  use ModCimiPlanet, ONLY: &
        nspec, amu_I
  use ModConst, ONLY: &
        cLightSpeed, cElectronCharge
  use ModImTime, ONLY: &
        StartTime, CurrentTime, iCurrentTime_I
  use ModTimeConvert, ONLY: &
        time_real_to_int
  use ModCimi, ONLY: &
        f2, dt, Time, energy
  use ModCimiGrid, ONLY: &
        MinLon, MaxLon, MinLonPar, MaxLonPar, &
        MinLonGlobal, MaxLonGlobal, &
        nt, np, neng, npit, nm, nk, phi, xlat, xmlt
  use ModCimiInitialize, ONLY: &
        xmm, xk, dmm, xjac
  use ModCimiTrace, ONLY: &
        fieldpara, brad => ro, bo, irm, ekev, iba, vel, &
        alscone, iw2, xmlto, ionization, sinA, pp, bm
  use ModDiagDiff, ONLY: &
        calc_DQQ, calc_Q2_only, interpol_D_coefK, &
        mapPSDtoQ, mapPSDtoE, diffuse_Q1, diffuse_Q2, &
        init_diag_diff, Eq
  use CIMI_waves, ONLY:	&
        UseDiagDiffusion, UseKpIndex, ReadDiffCoef, WavePower,&
        TypeChorus, TypeHiss, HIpower, CHpower
  use whistler_model, ONLY: init_whistlers, WPowTotal,&
        update_wave_power
  use ModImIndices, ONLY: &
        interpolate_ae, UseAeKyoto, UseKpApF107IndicesFile, &
        get_im_indices_Kp
  use ModIndicesInterfaces, ONLY: &
        get_kp
  use ModUtilities, ONLY: &
        CON_STOP
  use ModPlasmasphere, ONLY: &
        PlasDensity_C
  use ModCimiMethods, ONLY: &
        initial_f2, lossconeIM
  use ModCimiPlot, ONLY: &
        DoSavePlot, DoSaveFlux,DoSavePSD,&
        Cimi_plot_fls, Cimi_plot_psd
  use ModWaveDiff, ONLY:   &
        diffuse_aa, diffuse_EE

  implicit none

  ! regular variables
  real,intent(in):: DtAdvance
  integer, save :: ib0(nt)
  real:: flux(nspec,np,nt,neng,npit), psd(nspec,np,nt,nm,nk)
  real::  AE_temp = 0., Kp_temp = 0.
  integer:: iLat, iLon, iSpecies, iSat, iOperator
  logical:: IsFirstCall =.true.
  integer:: iError
  integer:: n,i,j,k,m

  dt = DtAdvance

  ! Update CurrentTime and iCurrentTime_I
  CurrentTime = StartTime+Time
  call time_real_to_int(CurrentTime,iCurrentTime_I)

  ! do field line integration and determine vel, ekev, momentum (pp), etc.
  call timing_start('cimi_fieldpara')
  !if ( IsFirstCall ) then
    call fieldpara(Time,dt,cLightSpeed,cElectronCharge,xlat,xmlt,phi,xk,&
                   .false.)
    ! TODO implement a better solution here
    irm(1) = 1 ! overwrite irm since it doesnt work right with single point
    iba(1) = 1 ! overwrite iba since it doesnt work right with single point
  ! overwrite iw2, it depends on irm so needs to be recalculated
  if ( IsFirstCall ) then
    do n=1,nspec
        do k=1,nk
            iw2(m,n)=nw
            find_iw2: do m=1,nw
                if (ekev(irm(1),1,m,k,n).gt.energy(neng,n)) then
                    iw2(m,n)=m
                    exit find_iw2
                endif
            enddo find_iw2
        enddo
    enddo ! nspec
  endif
  call timing_stop('cimi_fieldpara')

  if ((.not. IsFirstCall) .and. UseDiagDiffusion) then
    ! interpolate a0 of const. Q2 curve correspoding K
    call timing_start('cimi_interpol_D_coefK')
    call interpol_D_coefK
    call timing_stop('cimi_interpol_D_coefK')
  endif

  ! set the boundary and temperature and density (also sets the interior
  ! density and temperature for I.C. but that is only if initial_f2 is called)
  !call cimi_set_boundary_empirical
  if (UseKpIndex) then
    if (UseKpApF107IndicesFile) then
        call get_im_indices_Kp(CurrentTime, Kp_temp)
    else
        call get_kp(CurrentTime, Kp_temp, iError)
    endif
  elseif(UseAeKyoto) then
      call interpolate_ae(CurrentTime, AE_temp)
  else
      call CON_STOP('IM error: Kp not used.'//&
                    'Kp is needed for simple plasmasphere.')
  endif
  Kp_temp = 4.0

  ! Set up simple plasmasphere
  call timing_start('simple_plasmasphere')
  call simple_plasmasphere(Kp_temp)
  call timing_stop('simple_plasmasphere')


  !  read wave diffcoef 
  if (IsFirstCall)then
      call ReadDiffCoef(np,MinLon,MaxLon)
  endif

  call WavePower(np,Time,PlasDensity_C(:,MinLonGlobal:MaxLonGlobal),&
                 AE_temp,Kp_Temp,bo,brad,xmlto,iba,&
                 MinLon,MaxLon,MinLonPar,MaxLonPar)


  if (IsFirstCall .and. UseDiagDiffusion) then
    call init_diag_diff
    select case(TypeHiss)
    case('UCLA')
        call calc_DQQ
    case('BU')
        call calc_Q2_only
        call timing_start('cimi_interpol_D_coefK')
        call interpol_D_coefK
        call timing_stop('cimi_interpol_D_coefK')
    end select
    call init_whistlers(pp(1,1,:,1,3),bo(1,1))
  endif

  ! setup initial distribution
  if (IsFirstCall) then
     call initial_f2(nspec,np,nt,amu_I,xjac,ib0)
     IsFirstCall=.false.
  endif

  call cimi_output_flux(flux,psd)

  call update_wave_power(pp(1,1,:,1,3),sinA,&
                         psd(1,1,:,:,3)*1.61e19,&
                         PlasDesnsity_C(:,1,1),bo(1,1),ro(1,1)*6371e3,dt)

  write(*,*) 'HIpower, CHpower, Localpower: ',&
              HIpower(1,1), CHpower(1,1), WPowTotal

  HIpower(1,1) = WpowTotal
  CHpower(1,1) = WpowTotal

  call timing_start('cimi_WaveDiffusion')
  if ( UseDiagDiffusion ) then

    ! map PSD from M to Q2
    call timing_start('cimi_mapPSDtoQ')
    call mapPSDtoQ
    call timing_stop('cimi_mapPSDtoQ')

    ! calculate diffusion in Q2 at fixed Q1
    call timing_start('cimi_Diffuse_Q2')
    call diffuse_Q2
    call timing_stop('cimi_Diffuse_Q2')

    ! calculate diffusion in Q1 at fixed Q2
    call timing_start('cimi_Diffuse_Q1')
    call diffuse_Q1
    call timing_stop('cimi_Diffuse_Q1')

    ! map back PSD from Q2 to M
    call timing_start('cimi_mapPSDtoE')
    call mapPSDtoE
    call timing_stop('cimi_mapPSDtoE')
  else
    call timing_start('cimi_Diffuse_aa')
    call diffuse_aa(f2,dt,xjac,iba,iw2)
    call timing_stop('cimi_Diffuse_aa')

    call timing_start('cimi_Diffuse_EE')
    call diffuse_EE(f2,dt,xmm,xjac,iw2,iba)
    call timing_stop('cimi_Diffuse_EE')
  endif

  call timing_stop('cimi_WaveDiffusion')

  call timing_start('cimi_lossconeIM')
  call lossconeIM
  call timing_stop('cimi_lossconeIM')

  ! save the initial point (needs to be fixed and so not called)
  !if ( DoSavePlot .and. mod(Time,50.0).eq.0) then
  if ( DoSavePlot ) then
     call timing_start('cimi_output')
     call cimi_output_flux(flux,psd)
     call timing_stop('cimi_output')

     if ( mod(Time,50.0).eq.0) then
        do iSpecies = 1, nspec
            if ( DoSaveFlux( iSpecies ) ) &
                call Cimi_plot_fls( flux( iSpecies, :, :, :, : ), &
                                    iSpecies, time )
            if ( DoSavePSD( iSpecies ) ) &
                call Cimi_plot_psd( psd( iSpecies, :, :, :, : ), &
                                    iSpecies, time, xmm, xk )
        enddo
     endif 
  endif

  ! Update time
  Time = Time+dt

end subroutine diffusion_test_run
!==============================================================================
subroutine cimi_output_flux(flux,psd)

  use ModCimiPlanet, ONLY: nspec, amu_I
  use ModConst,	ONLY: cLightSpeed, cElectronCharge
  use ModCimi, ONLY: f2,energy
  use ModCimiGrid, ONLY: nt, np, neng, npit, nm, nk, sinAo,&
                         MinLonPar,MaxLonPar
  use ModCimiInitialize, ONLY: xmm, dmm, xjac
  use ModCimiTrace, ONLY: ekev, sinA, pp, iba
  use ModCimiUtil, ONLY: lintp2aIM

  implicit none

  real, intent(inout):: &
       flux(:,:,:,:,:), psd(:,:,:,:,:)
  ! Locals
  real:: ekev2D(nm,nk), sinA1D(nk), aloge(neng,nspec), flux2D(nm,nk)

  real:: amu1, flx,flx_lo,psd1
  integer:: i, j, k, m, n

  !----------------------------------------------------------------------------
  ! Set floor values
  flux=1e-50
  psd =1e-50
  ! Get log of energy
  do n=1,nspec
    aloge(1:neng,n) = log10(energy(n,1:neng))
  enddo
  longloop:do j=MinLonPar,MaxLonPar
    latloop:do i=1,iba(j)
        speciesloop:do n=1,nspec
            amu1=amu_I(n)**1.5
            do k=1,nk
              sinA1D(k) = sinA(i,j,k)
            enddo
            ! Loop over M-K space for PSD output
            do k=1,nk
                do m=1,nm
                    psd1=f2(i,j,m,k,n)/1.e20/1.e19/xjac(i,m,n)! mug^-3cm^-6s^3
                    flx=psd1*(1.6e19*pp(i,j,m,k,n))*pp(i,j,m,k,n)
                    flux2D(m,k)=-50.
                    if (flx > 1.e-50) flux2D(m,k)=log10(flx)
                    ekev2D(m,k)=log10(ekev(i,j,m,k,n))
                    psd(i,j,m,k,n)=psd1
                enddo
            enddo

            ! Loop over E-alpha space for Flux output (interpolated)
            do k=1,npit
                do m=1,neng
                    ekev2D(m,k)=log10(ekev(i,j,m,k,n))
                    call lintp2aIM(ekev2D,sinA1D,flux2D,nm,nk,aloge(m,n),&
                                   sinAo(k),flx_lo)
                    flux(i,j,m,k,n)=10.**flx_lo
                enddo
            enddo
        enddo speciesloop
    enddo latloop
  enddo longloop
end subroutine cimi_output_flux
!==============================================================================
end module ModDiffusionTests

