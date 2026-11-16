module CIMI_LocalWave
    implicit none
    save
    private ! except
    
    ! Global Constants
    real, parameter, public :: cC = 2.9979e8       ! speed of light (m/s)
    real, parameter, public :: cMe = 9.1094e-31    ! electron mass (kg)
    real, parameter, public :: cPi = 3.14159265359
    real, parameter, public :: ckb = 1.381e-23
    real, parameter, public :: cEta0 = 8.854e-12   ! permittivity of free space
    real, parameter, public :: cQe = 1.602e-19     ! electron charge (coulombs)
    real, parameter, public :: ckeV2Joule = cQe*1e3! keV to Joules
    real, parameter, public :: ccm2m = 1e-4        ! cm^2 to m^2
    
    ! Grid & Simulation parameters (Manually set in this module)
    ! NOTE W_I indicates array with shape of of resonant frequency grid
    !      P_I indicates array with shape of perpendicular momentum grid
    integer,public :: nFreq=181, nPperp=100, nPpar=100
    real, allocatable,public :: wResW_I(:), pperpP_I(:), pparP_I(:),&
                                rIonoW_I(:)
    real, allocatable, save  :: lnp(:,:), lnf(:,:), sina(:)
    real :: FreqGridFactor=0.95
    real :: dlnpperp, dlnppar

    ! Grid & Simulation parameters (To be allocated outside the module)
    !integer,public :: n_energy, n_alpha
    real, public :: tSim, dt
    
    ! Global State Variables
    real, allocatable, public :: WPowW_I(:)          ! Wave power
    real, allocatable, public :: WPowGrownW_I(:)     ! Power generated
    real, allocatable, public :: WPowLostW_I(:)      ! Power lost
    real, public :: WPowTotal                        ! Wave power in pT^2
    real, parameter, public :: WPowFloor=1e-4

    ! Interal variables
    real, dimension(181) :: v_group, eta_r, a_r,  gain_rate, loss_rate
    
    public :: init_whistlers
    public :: update_wave_power
    public :: calc_derivativesPA
    public :: momentum, energy, wplasma, wgyro,&
              calc_kr, calc_gamma_lorentz, calc_pR, calc_delR,calc_norm1,&
              calc_growth_rate
    public :: write_resonances,write_verbose_output
    public :: unit_test_localwave

    ! Flag for berbose output files
    logical, save :: isFirstPSDWrite=.true.
    logical, save :: isFirstResonanceWrite=.true.
    
contains
    ! ======================================================================
    ! Initialize Whistlers
    ! ======================================================================
    subroutine init_whistlers(pp,sina,N,bo)
        use ModCimiUtil, ONLY: logspace

        real,intent(in) :: pp(:,:),sina(:)
        real, intent(in) :: N,bo
        integer :: ipperp,ippar,iFreq
        real :: wp,wg,ppmax

        ! allocate wave grid and momentum grid
        allocate(wResW_I(nFreq))
        allocate(rIonoW_I(nFreq))
        allocate(pperpP_I(nPperp))
        allocate(pparP_I(nPpar))
        allocate(WPowW_I(nFreq))
        allocate(WPowGrownW_I(nFreq))
        allocate(WPowLostW_I(nFreq))

        ! Populate the resonance frequency grid
        wg = wgyro(bo)
        wp = wplasma(N)
        do iFreq=1,nFreq
            wResW_I(iFreq) = (0.5*wg) * FreqGridFactor**(iFreq-1)
        enddo
        write(*,*) 'wplasma: ', wp
        write(*,*) 'wgyro: ', wg
        write(*,*) 'wr range: ',minval(wResW_I),maxval(wResW_I)

        ! Populate regular ppar and pperp grid
        ppmax = maxval(pp)
        call logspace(minval(pp),ppmax,nPpar,pparP_I,dlnppar)
        call logspace(minval(pp),ppmax,nPperp,pperpP_I,dlnpperp)

        WPowW_I = 1.0e-2 ! Small seed wave power

        ! Initialize rIonoW_I array (e.g., 5% reflection -> 0.05)
        rIonoW_I = 0.05
        
        call freq_to_energy()
    end subroutine init_whistlers

    ! ======================================================================
    ! Update Wave Power
    ! ======================================================================
    subroutine update_wave_power(pp,sina,PSD,&
                                 N_dens, bo, b_length,dt,isVerbose)
        ! Inputs
        real, dimension(:),intent(in) :: sina
        real, dimension(:,:),intent(in) :: pp,PSD
        real, intent(in) :: N_dens, bo, b_length
        real, intent(in) :: dt
        logical, intent(inout), optional :: isVerbose

        ! Locals
        integer :: i
        real,dimension(nFreq) :: growth_rate
        real :: wp, wg, k_res, gamma_L, p_res, del_R, a_crit

        ! Default to isVerbose=False
        if (.not.present(isVerbose)) isVerbose=.false.

        ! 0. update module level variable log(PSD) for use
        lnf = log(PSD)
        
        ! 1. Plasma properties (Functions I & II)
        wp = wplasma(N_dens)
        wg = wgyro(bo)
        
        ! Loop over all frequencies
        if (isVerbose) then
            write(*,"(7a11)") "wr/wg","eta","A","growth",&
                                                       "vgrp","loss","Wpow"
        endif
        do i = 1, nFreq

            ! Get the portion of the distribution in cyclotron resonance
            call calc_eta_rel(pp,sina,PSD,N_dens,wResW_I(i),wp,wg,eta_r(i))

            ! Get the relativistic pitch angle anisotropy factor
            call calc_a_rel(pp,sina,PSD,wResW_I(i),wp,wg,a_r(i))

            ! Get the growth rate
            growth_rate(i) = calc_growth_rate(wResW_I(i),wp,wg,eta_r(i),a_r(i))

            ! Set hard limits on growth_rate values
            !if (growth_rate(i).lt.0.0) growth_rate(i)=0.0
            !if (growth_rate(i).gt.5.0) growth_rate(i)=5.0

            ! Get the group velocity and loss rate
            ! NOTE growth is Wave amp ratio, power ~ A^2
            v_group(i) = calc_vgroup(wResW_I(i),wp,wg)
            gain_rate(i) = 2*growth_rate(i)
            loss_rate(i) = 2*v_group(i)/b_length * log(rIonoW_I(i))

            ! Update Wave power state
            WPowGrownW_I(i) = WPowW_I(i) * (exp(gain_rate(i)*dt)) - WPowW_I(i)
            WPowLostW_I(i) =  WPowW_I(i) * (exp(loss_rate(i)*dt)) - WPowW_I(i)
            WPowW_I(i) = WPowW_I(i) * exp((gain_rate(i)+loss_rate(i)) * dt)
            
            ! Set floor if  wave power drops too low
            if (WPowW_I(i) .lt. WPowFloor) WPowW_I(i) = WPowFloor

            if (isVerbose) then
                write(*,"(7e11.3)") wResW_I(i)/wg,eta_r(i),a_r(i),&
                                gain_rate(i),v_group(i),loss_rate(i),WPowW_I(i)
            endif
        end do

        call freq_to_energy()
    end subroutine update_wave_power

    ! ======================================================================
    ! Frequency to Energy Integration
    ! ======================================================================
    subroutine freq_to_energy()
        real, dimension(nFreq) :: frW_I
        integer :: iFreq
        ! Integrate over frequency space and obtain total wave power in pT^2
        frW_I = wResW_I/(2*cPi)
        WPowTotal = 0.0
        do iFreq=2,nFreq-1
            WPowTotal = WPowTotal + (WPowW_I(iFreq+1)+WPowW_I(iFreq-1))/2 * &
                                            (frW_I(iFreq-1)-frW_I(iFreq+1))
        enddo
    end subroutine freq_to_energy

    ! ======================================================================
    ! Portion of distribution in cyclotron resonance (Eq A2 integration)
    ! ======================================================================
    subroutine calc_eta_rel(pp,sina,PSD,N,&
                            wResW_I,wp,wg,&
                            eta_r)
        ! Inputs
        real, dimension(:), intent(in) :: sina       ! On k grid
        real, dimension(:,:), intent(in) :: pp, PSD  ! On m-k grid
        real, intent(in) :: N, wResW_I, wg, wp       ! Dens, Wave frequencies
        ! Outputs
        real, intent(out) :: eta_r
        ! NOTE Global variables modified: None
        !                     unmodified: pperpP_I, pperpP_I
        ! Locals
        real :: kr, p_res, del_R, dpperp
        real :: integral,dfdpparR,dfdpperpR,dfdpperpRmax,dfdpperpRmin
        integer :: j
        real :: f
        
        integral = 0.0
        ! Numerical integration over perpendicular momentum (pperp)
        !write(*,"(a4,5a12)") "j","pperpP_I(j)","f",&
        !                  "dpperp","dfdpperpR","integral"
        dfdpperpRmax = 0.0
        dfdpperpRmin = 0.0
        do j = 1,  SIZE(pperpP_I)-1
            dpperp = pperpP_I(j+1)-pperpP_I(j)
            ! Find strictly parallel resonant momentum at this freq & pperp
            kr = calc_kr(wResW_I,wp,wg)
            p_res = calc_pR(wResW_I,wp,wg,pperpP_I(j))
            del_R = calc_delR(wResW_I,wp,wg,pperpP_I(j))
            ! Get derivative interpolated to ppar = p_res
            call calc_derivativesPA(pp,sina,PSD/pp**2,&
                                    abs(p_res),pperpP_I(j),dlnppar,dlnpperp,&
                                    f,dfdpparR,dfdpperpR)
            if (dfdpperpR.gt.dfdpperpRmax) dfdpperpRmax=dfdpperpR
            if (dfdpperpR.lt.dfdpperpRmin) dfdpperpRmin=dfdpperpR
            ! Construct the integrand
            integral = integral + (pperpP_I(j)**2/del_R) * dpperp *dfdpperpR&
                                / (ccm2m*ckeV2Joule)
            !write(*,"(i4,5e12.3)") j, pperpP_I(j), f,&
            !                       dpperp, dfdpperpR, integral
        end do
        
        eta_r = (cPi * cMe * (wResW_I - wg)) / (N * kr) * integral
    end subroutine calc_eta_rel

    ! ======================================================================
    ! Anisotropy factor (Eq A3 integration)
    ! ======================================================================
    subroutine calc_a_rel(pp,sina,PSD,&
                          wResW_I,wp,wg,&
                          a_r)
        implicit none
        ! Inputs
        real, dimension(:), intent(in) :: sina        ! on k grid
        real, dimension(:,:), intent(in) :: pp,PSD    ! on m-k grid
        real, intent(in) :: wResW_I, wg, wp           ! Dens, Wave frequencies
        ! Outputs
        real, intent(out) :: a_r
        ! Locals
        real :: kr, p_res, del_R, gamma_L, dpperp
        real :: integral1,integral2,dfdpparR,dfdpperpR
        integer :: j
        real :: f

        integral1 = 0.0
        integral2 = 0.0
        ! Numerical integration over perpendicular momentum (pperp)
        do j = 1,  nPperp-1
            dpperp = pperpP_I(j+1)-pperpP_I(j)
            ! Find strictly parallel resonant momentum at this freq & pperp
            kr = calc_kr(wResW_I,wp,wg)
            p_res = calc_pR(wResW_I,wp,wg,pperpP_I(j))
            del_R = calc_delR(wResW_I,wp,wg,pperpP_I(j))
            gamma_L = calc_gamma_lorentz(wResW_I,wp,wg,pperpP_I(j))
            ! Get derivative interpolated to ppar = p_res
            call calc_derivativesPA(pp,sina,PSD/pp**2,&
                                    abs(p_res),pperpP_I(j),dlnppar,dlnpperp,&
                                    f,dfdpparR,dfdpperpR)
            ! Construct the integrand(s)
            integral1 = integral1 + &
                            (pperpP_I(j)**2/del_R) * dpperp * dfdpperpR&
                            / (ccm2m*ckeV2Joule)
            integral2 = integral2 + &
                            (pperpP_I(j)**2/(del_R*gamma_L)) * dpperp* &
                              (pperpP_I(j)*dfdpparR - abs(p_res)*dfdpperpR)&
                            / (ccm2m*ckeV2Joule)
        end do
        
        a_r = kr/(cMe*(wResW_I-wg)) * integral2 / integral1
    end subroutine calc_a_rel

    ! ======================================================================
    ! Utility Subroutines and Functions
    ! ======================================================================
    subroutine calc_derivativesPA(pp,sina,f,&
                                  ppar_0,pperp_0,dlnppar,dlnpperp,&
                                  f_0,dfdppar,dfdpperp)
    use ModCimiUtil, only: lintp2aIM
    real, intent(in)  :: pp(:,:),sina(:),f(:,:)
    real, intent(in)  :: ppar_0,pperp_0,dlnppar,dlnpperp
    real, intent(out) :: f_0,dfdppar,dfdpperp
    integer :: nx, ny
    real :: hpar0, hperp0, g0, gl, gr, gu, gd

    ! Derivative in log-log scale: g=log(f), h=log(x);
    !                       dg/dh =   (dg/dx)   / (dh/dx) 
    !                             = (1/f df/dx) / (1/x)
    !
    !            rearrange: df/dx = f/x * dg/dh
    !                               -----------------
    !                             = |exp(g-h) * dg/dh|
    !                               -----------------

    nx = size(f,1)
    ny = size(f,2)

    hpar0 = log(ppar_0)
    hperp0 = log(pperp_0)

    ! Interpolate g=log(f) at five locations
    g0 = g_at(hpar0, hperp0)
    gl = g_at(hpar0-dlnppar, hperp0)
    gr = g_at(hpar0+dlnppar, hperp0)
    gd = g_at(hpar0, hperp0-dlnpperp)
    gu = g_at(hpar0, hperp0+dlnpperp)

    ! Calculate derivatives on log-log scale
    f_0 = exp(g0)
    dfdppar   = exp(g0-hpar0)  * (gr-gl)/(2*dlnppar)
    dfdpperp  = exp(g0-hperp0) * (gu-gd)/(2*dlnpperp)

    contains
        real function g_at(hpar, hperp)
            real, intent(in) :: hpar, hperp
            real :: sina_target, lnp_target, g
            sina_target = (exp(2*hperp)) / (exp(2*hperp)+exp(2*hpar))
            lnp_target = 0.5 * log(exp(2*hpar)+exp(2*hperp))
            call lintp2aIM(log(pp), sina, log(f), nx, ny,&
                            lnp_target, sina_target, g)
            g_at = g
        end function g_at
    end subroutine

    !-------------------------------------------------------------------------
    elemental real function momentum(energy)
        real,intent(in) :: energy
        momentum = sqrt(energy*cMe*((energy/(cMe*cC**2))+2))
    end function momentum

    !-------------------------------------------------------------------------
    elemental real function energy(momentum)
        real,intent(in) :: momentum
        energy = sqrt( (cMe*cC**2)**2 + (momentum*cC)**2 ) - cMe*cC**2
    end function energy

    !-------------------------------------------------------------------------
    real function wplasma(N_dens)
        real, intent(in) :: N_dens
        wplasma = sqrt(N_dens * cQe**2 / (cMe*cEta0))
    end function wplasma

    !-------------------------------------------------------------------------
    real function wgyro(bo)
        real, intent(in) :: bo
        wgyro = abs((cQe * bo) / (cMe))
    end function wgyro

    !-------------------------------------------------------------------------
    real function calc_kr(w_r, wp, wg)
        real, intent(in) :: w_r, wp, wg
        calc_kr = (w_r / cC) * sqrt(1.0 + (wp**2) / (w_r * (wg - w_r)))
    end function calc_kr

    !-------------------------------------------------------------------------
    real function calc_gamma_lorentz(wResW_I, wp, wg, pperp)
        real, intent(in) :: wResW_I, wp, wg, pperp
        real :: kr, term1, term2, term3
        ! Calculate gamma_R based on momentum components
        kr = calc_kr(wResW_I,wp,wg)
        term1 = (cC * kr / wResW_I)**2 - 1.0
        term2 = 1.0 + (pperp / (cMe * cC))**2
        term3 = (wResW_I / wg)**2
    calc_gamma_lorentz = &
                (-1.0+(cC * kr/wResW_I) * sqrt(term1 * term2 * term3 + 1))&
                                 /(term1 * (wResW_I / wg))
    end function calc_gamma_lorentz

    !-------------------------------------------------------------------------
    real function calc_pR(wResW_I, wp, wg, pperp)
        real,intent(in) :: wResW_I,wp,wg,pperp
        real :: kr,gamma_L
        kr = calc_kr(wResW_I,wp,wg)
        gamma_L = calc_gamma_lorentz(wResW_I,wp,wg,pperp)
        calc_pR = (cMe / kr) * (gamma_L * wResW_I - wg)
    end function calc_pR

    !-------------------------------------------------------------------------
    real function calc_delR(wResW_I,wp,wg,pperp)
        real, intent(in) :: wResW_I,wp,wg,pperp
        real :: kr,gamma_L,pR

        kr = calc_kr(wResW_I,wp,wg)
        gamma_L = calc_gamma_lorentz(wResW_I,wp,wg,pperp)
        pR = calc_pR(wResW_I,wp,wg,pperp)
        calc_delR = 1.0 - (wResW_I * pR)/(cMe * cC**2 * kr * gamma_L)
    end function calc_delR

    !-------------------------------------------------------------------------
    real function calc_Acritical(wResW_I, wg)
        real, intent(in) :: wResW_I, wg
        calc_Acritical = wResW_I / (wg - wResW_I)
    end function calc_Acritical

    !-------------------------------------------------------------------------
    real function calc_vgroup(wResW_I, wp, wg)
        real, intent(in) :: wResW_I, wp, wg
        real :: kr, num, den
        kr = calc_kr(wResW_I,wp,wg)
        num = 2.0 * cC**2 * kr * wg * (1.0-wResW_I/wg)**2
        den = wp**2 + 2.0 * wg * wResW_I * (1.0 - wResW_I/wg)**2
        calc_vgroup = num / den
    end function calc_vgroup

    !-------------------------------------------------------------------------
    real function calc_growth_rate(wResW_I,wp,wg,eta_r, a_r)
        real, intent(in) :: wResW_I, wp, wg, eta_r, a_r
        real :: a_crit,num,denom
        a_crit = calc_Acritical(wResW_I,wg)
        num = cPi * wp**2 * eta_r * (a_r - a_crit)
        denom = 2 * wResW_I + (wp**2 * wg) / (wResW_I - wg)**2
        calc_growth_rate = num / denom
    end function calc_growth_rate

    !-------------------------------------------------------------------------
    real function calc_norm1(wResW_I,wp,wg,N)
        real, intent(in) :: wResW_I,wp,wg,N
        real :: kr
        kr = calc_kr(wResW_I,wp,wg)
        calc_norm1 = cPi*cMe*(wResW_I-wg)/(N*kr)
    end function calc_norm1

    !-------------------------------------------------------------------------
    real function calc_norm2(wResW_I,wp,wg)
        real, intent(in) :: wResW_I,wp,wg
        real ::  kr
        kr = calc_kr(wResW_I,wp,wg)
        calc_norm2 =  kr/(cMe*(wResW_I-wg))
    end function calc_norm2

    !-------------------------------------------------------------------------

    subroutine write_resonances(nsteps, tSim, N_dens, bo)
        use ModIoUnit, ONLY: UnitTmp_

        integer, intent(in) :: nsteps
        real, intent(in) :: tSim
        real, intent(in) :: N_dens,bo

        integer :: ifreq,ipperp
        real :: wg, wp
        real,dimension(nPperp) :: delR, gammaL, pr

        wp = wplasma(N_dens)
        wg = wgyro(bo)

        ! Write out pR and delR as function of resonant frequency
        if (isFirstResonanceWrite) then
            open(unit=UnitTmp_,&
                file='IM/plots/localwave_resonances.dat',&
                status='unknown')
            ! Header info
            write(UnitTmp_,*)"! ntimes, nfreq, npperp"
            write(UnitTmp_,"(3i6)")nsteps,nFreq,nPperp
            write(UnitTmp_,*)"! resonant frequencies"
            write(UnitTmp_,"(7e14.3e3)") wResW_I
            write(UnitTmp_,*)"! perpendicular momentum grid"
            write(UnitTmp_,"(7e14.3e3)") pperpP_I
        else
            open(unit=UnitTmp_,&
                file='IM/plots/localwave_resonances.dat',&
                status='old',&
                position='append')
        endif
        ! Entry info
        write(UnitTmp_,*)"! tSim, wp, wg"
        write(UnitTmp_,"(3f10.1)") tSim, wp, wg
        do ifreq=1,nFreq
            write(UnitTmp_,"(5a14)") "! wResW_I","kres","WPowW_I","eta_r","a_r"
            write(UnitTmp_,"(5e14.3e3)") wResW_I(ifreq),&
                                         calc_kr(wResW_I(ifreq),wp,wg),&
                                         WPowW_I(ifreq),&
                                         eta_r(ifreq), a_r(ifreq)
            do ipperp=1,nPperp
                gammaL(ipperp) = calc_gamma_lorentz(wResW_I(ifreq),wp,wg,&
                                                              pperpP_I(ipperp))
                pR(ipperp) = calc_pR(wResW_I(ifreq),wp,wg,pperpP_I(ipperp))
                delR(ipperp) = calc_delR(wResW_I(ifreq),wp,wg,pperpP_I(ipperp))
            enddo
            write(UnitTmp_,"(1a14)") "! gammaL"
            write(UnitTmp_,"(7e14.3e3)") gammaL
            write(UnitTmp_,"(1a14)") "! pR"
            write(UnitTmp_,"(7e14.3e3)") pR
            write(UnitTmp_,"(1a14)") "! delR"
            write(UnitTmp_,"(7e14.3e3)") delR
        enddo
        if (isFirstResonanceWrite) then 
            write(*,*)"SAVED: IM/plots/localwave_resonances.dat"
            isFirstResonanceWrite=.false.
        endif
        close(UnitTmp_)
    end subroutine write_resonances

    !-------------------------------------------------------------------------

    subroutine write_verbose_output(nsteps,tSim,N_dens, bo, pp, sina, PSD)
        ! Inputs
        integer, intent(in) :: nsteps
        real, intent(in) :: tSim
        real, intent(in) :: N_dens,bo
        real, dimension(:), intent(in) :: sina        ! on k grid
        real, dimension(:,:), intent(in) :: pp,PSD    ! on m-k grid
        ! Locals
        real,dimension(nPpar,nPperp) :: PSD_out,dfdppar_out,dfdpperp_out
        character(len=100) :: PSDFile,dfdpparFile,dfdpperpFile

        PSDFile = "IM/plots/localwave_PSDinterp.psd"
        dfdpparFile = "IM/plots/localwave_dfdppar.psd"
        dfdpperpFile = "IM/plots/localwave_dfdpperp.psd"

        call interpolate_full_grid(PSD_out,dfdppar_out,dfdpperp_out)
        call write_interp_psd(PSDFile, PSD_out)
        call write_interp_psd(dfdpparFile, dfdppar_out)
        call write_interp_psd(dfdpperpFile, dfdpperp_out)

        isFirstPSDWrite = .false.

    contains
        subroutine interpolate_full_grid(PSD_out, dfdppar_out, dfdpperp_out)
            use ModCimiUtil, only: lintp2aIM
            ! Outputs
            real,dimension(nPpar,nPperp),intent(out) :: PSD_out,&
                                                        dfdppar_out,&
                                                        dfdpperp_out
            ! Locals
            integer :: ippar,ipperp
            real :: f,logpsd
            do ipperp=1,nPperp
                do ippar=1,nPpar
                    ! Interpolate PSD to output grid
                    !call lintp2aIM(log(pp), sina, log(PSD/pp**2),&
                    !               size(PSD,1),size(PSD,2),&
                    !         log(hypot(pparP_I(ippar),pperpP_I(ipperp))),&
                    !           atan2(pperpP_I(ipperp),pparP_I(ippar)),&
                    !                logpsd)
                    ! Get derivatives on output grid
                    call calc_derivativesPA(pp,sina,PSD/pp**2,&
                       pparP_I(ippar),pperpP_I(ipperp),dlnppar,dlnpperp,& ! IN
                                        f,dfdppar_out(ippar,ipperp),&
                                        dfdpperp_out(ippar,ipperp))     ! OUT
                    PSD_out(ippar,ipperp) = f*(pperpP_I(ipperp)**2+&
                                               pparP_I(ippar)**2)

                    if (abs(PSD_out(ippar,ipperp)).lt.1e-50) then
                        PSD_out(ippar,ipperp)=PSD_out(ippar,ipperp)
                    endif
                    if (abs(dfdppar_out(ippar,ipperp)).lt.1e-50) then
                        dfdppar_out(ippar,ipperp)=sign(1e-50,&
                                                    dfdppar_out(ippar,ipperp))
                    endif
                    if (abs(dfdpperp_out(ippar,ipperp)).lt.1e-50) then
                        dfdpperp_out(ippar,ipperp)=sign(1e-50,&
                                                   dfdpperp_out(ippar,ipperp))
                    endif
            enddo
        enddo
        end subroutine interpolate_full_grid

        subroutine write_interp_psd(NameFile, PSD_arr)
            use ModIoUnit, ONLY: UnitTmp_

            ! Inputs
            character(len=100), intent(in) :: NameFile
            real, dimension(nPpar,nPperp), intent(in) ::PSD_arr
            ! Locals
            integer :: ippar
            !write(*,*) 'FILE NAME IS: ',trim(NameFile)

            ! Write out PSD_plot to file
            if (isFirstPSDWrite) then
                open(unit=UnitTmp_,&
                    file=trim(NameFile),&
                    status='unknown')
                ! Header info
                write(UnitTmp_,*)"! ntimes, nppar, npperp"
                write(UnitTmp_,"(3i6)")nsteps,nPpar,nPperp
                write(UnitTmp_,*)"! parallel momentum grid"
                write(UnitTmp_,"(7e14.3e3)") pparP_I
                write(UnitTmp_,*)"! perpendicular momentum grid"
                write(UnitTmp_,"(7e14.3e3)") pperpP_I
            else
                open(unit=UnitTmp_,&
                    file=trim(NameFile),&
                    status='old',&
                    position='append')
            endif
            ! Entry info
            write(UnitTmp_,*)"! tSim"
            write(UnitTmp_,"(f6.1)") tSim
            do ippar=1,nPpar
                write(UnitTmp_,"(7e14.3e3)") PSD_arr(ippar,:)
            enddo
            if (isFirstPSDWrite) write(*,*)"SAVED: ",trim(NameFile)
            close(UnitTmp_)
        end subroutine write_interp_psd

    end subroutine write_verbose_output

    ! ======================================================================
    ! Unit tests:
    !   > sets up synthetic PSD using analytical expression
    !   > adds 'unstable' gaussian population
    !   > initializes arrays using a prescribed grid
    !   > calcuates single value targets spot checking basic parameters
    !   > calculates integrated targets building up to growth rate
    !   > loops through some time steps allowing growth rate to affect Wpow
    !   > saves the PSD used in detail and as a simplified .fin
    ! ======================================================================
    subroutine unit_test_localwave
        use ModCimiUtil, ONLY: linspace,logspace
        integer :: step, iE, iA
        integer, parameter :: max_steps = 10
    
        ! Dummy parameters for execution
        real :: N_dens = 16e6       ! m^-3
        real :: bo = 149e-9             ! T (149 nT)
        real :: b_length = 5.0 * 6371e3     ! Interaction length (m)
        real :: E_min=1.e-2, E_max=1.e4
        real :: C1,kT,gamma1,S,E0,gamma2 ! analytical PSD parameters
        real :: PSD0_inj,ppar0_inj,pperp0_inj ! injection PSD
        logical :: doVerbose
        logical,save :: IsFirstTestCall = .true.
        integer, parameter :: nM=48, nK=48
        integer :: iM,iK
        real :: Ekev_MK(nM,nK),sina_K(nK),PSD_MK(nM,nK), PSD_inj(nM,nK)
        real :: pp_MK(nM,nK), ppar_MK(nM,nK), pperp_MK(nM,nK)
    
        ! Timing variables used just for the test
        dt = 1.0
        tSim = 0.0

        ! Analytical PSD function (stable) parameters
        C1 = 2.38e6
        kT = 0.001
        gamma1 = 0.978364
        E0 = 1748.57
        gamma2 = 7.03623
        S = 0.3

        ! Injection PSD function (gaussian) parameters
        PSD0_inj = 1e4
        ppar0_inj = momentum(ckeV2Joule*1e3) ! 1 MeV
        pperp0_inj = momentum(ckeV2Joule*1e3) ! 1 MeV

        ! Set energy and pitch angle grid (will be from CIMI)
        call linspace(0.01,0.99,nK,sina_K) ! sin(alpha) in linear scale
        do iK=1,nK
            ! Set Energy 
            call logspace(E_min,E_max,nM,Ekev_MK(:,iK)) ! Energy in log scale
            pp_MK(:,iK)  = momentum(Ekev_MK(:,iK)*ckeV2Joule)
            ppar_MK(:,iK)  = pp_MK(:,iK) * sqrt(1-sina_K(iK)**2)
            pperp_MK(:,iK) = pp_MK(:,iK) * sina_K(iK)
        enddo

        ! Get a sample distribution
        call synthetic_PSD(N_dens,bo,C1,kT,gamma1,S,E0,gamma2,Ekev_MK,sina_K,&
                            PSD_MK)

        ! Add an injection distribution
        call injection_PSD(ppar0_inj,pperp0_inj,PSD0_inj,&
                           ppar_MK, pperp_MK, PSD_inj)
        PSD_MK = PSD_MK + PSD_inj

        ! Spot check values and frequency integrals
        call single_value_tests
        call integration_tests

        ! Initialize wave grid
        call init_whistlers(pp_MK, sina_K, N_dens, bo)

        ! Main Time Loop
        call write_verbose_output(max_steps+1,0.0,N_dens,bo,&
                                  pp_MK,sina_K,PSD_MK)
        call write_test_psd(max_steps+1)
        call write_resonances(max_steps+1,tSim,N_dens,bo)
        write(*,"(a,f12.2,a,e12.3,a,e12.3)") "Time: ",tSim,&
                                             "   Max W: ",maxval(WPowW_I),&
                                             "   W_pT2: ",WPowTotal
        doVerbose=.true.
        do step = 1, max_steps
            tSim = tSim + dt
        
            call update_wave_power(pp_MK,sina_K,PSD_MK,N_dens,bo,b_length,dt,&
                               doVerbose)
            call write_verbose_output(max_steps+1,0.0,N_dens,bo,&
                                      pp_MK,sina_K,PSD_MK)
            call write_test_psd(max_steps+1)
            call write_resonances(max_steps+1,tSim,N_dens,bo)
        
            write(*,"(a,f12.2,a,e12.3,a,e12.3)") "Time: ",tSim,&
                                             "   Max W: ",maxval(WPowW_I),&
                                             "   W_pT2: ",WPowTotal
            doVerbose=.false.
            IsFirstTestCall=.false.
        end do

        ! Save the synthetic PSD (with injection) as .psd and .fin
        call write_test_fin
    contains
        !-----------------------------------------------------------------
        subroutine synthetic_PSD(N_dens,bo,C1,kT,gamma1,S,E0,gamma2,Ekev,&
                             sina,PSD)
            real,intent(in) :: N_dens,bo,C1,kT,gamma1,S,E0,gamma2
            real,intent(in) :: Ekev(:,:), sina(:)
            real,intent(out) :: PSD(:,:)

            integer :: iM,iK

            do iK = 1,SIZE(sina)
                do iM = 1,SIZE(Ekev,1)
                    PSD(iM,iK) =  &
                     C1*Ekev(iM,iK)*(kT*(gamma1+1)+Ekev(iM,iK))**(-gamma1-1)/&
                       (1+(Ekev(iM,iK)/E0)**(gamma2)) * sina(iK)**(2*S)
                enddo
            enddo
        end subroutine synthetic_PSD
        !-----------------------------------------------------------------
        subroutine injection_PSD(ppar0,pperp0,PSD0,ppar,pperp,PSD)
            real, intent(in) :: ppar0,pperp0,PSD0
            real, intent(in),dimension(:,:) :: ppar,pperp
            real, intent(out),dimension(:,:) :: PSD
            real :: pp0
            integer :: iM,iK

            pp0 = sqrt(ppar0**2+pperp0**2)
            do iM = 1,SIZE(ppar,1)
                do iK = 1,SIZE(ppar,2)
                    PSD(iM,iK) = PSD0 *&
                        exp( -( ((pperp(iM,iK)-pperp0)/(pperp0))**2 ) &
                            -( ((ppar(iM,iK)-ppar0)/(2*ppar0))**2 ) )
                enddo
            enddo
        end subroutine injection_PSD
        !-----------------------------------------------------------------
        subroutine single_value_tests
            real :: wp,wg
            wp = wplasma(N_dens)
            wg = wgyro(bo)
            ! Quick tests
            write(*,*) "======== begin single value tests ========"
            write(*,"(a,e10.3)") " test wp: ", wp
            write(*,"(a,e10.3)") " test wg: ", wg
            write(*,"(a,e10.3)") " test kr|wResW_I=0.1wg: ", &
                                                    calc_kr(0.1*wg,wp,wg)
            write(*,"(a,e10.3)") " test wResW_I/(kr c)|wResW_I=0.1wg: ", &
                                           0.1*wg/(calc_kr(0.1*wg,wp,wg)*cC)
            write(*,"(a,e10.3)") " test gammaL|wResW_I=0.1wg, pperpP_I=0: ", &
                                        calc_gamma_lorentz(0.1*wg,wp,wg,0.0)
            write(*,"(a,e10.3)") " test pR|wResW_I=0.1wg, pperpP_I=0: ", &
                                                   calc_pR(0.1*wg,wp,wg,0.0)
            write(*,"(a,e10.3)") " test delR|wResW_I=0.1wg, pperpP_I=0: ", &
                                                 calc_delR(0.1*wg,wp,wg,0.0)
            write(*,"(a,e10.3)") " test norm1|wResW_I=0.1wg, N=16: ", &
                                               calc_norm1(0.1*wg,wp,wg,16.0)
            write(*,"(a,e10.3)") " test vgroup/c|wResW_I=0.1wg: ", &
                                                calc_vgroup(0.1*wg,wp,wg)/cC
            write(*,*) "======== end single value tests ========"
        end subroutine single_value_tests
        !-----------------------------------------------------------------
        subroutine integration_tests
            ! Tests interpolating and integrating over frequency space
            real,dimension(400) :: pperp
            real :: wr,wp,wg,kr,dpperp,integral,p_res,del_R,dlnpperp,&
                    dfdpparR,dfdpperpR,eta_r,integral2,a_r,gamma_L
            real :: p0,integrand,pRtest,growth
            integer :: ip
            real :: f
    
            wp = wplasma(N_dens)
            wg = wgyro(bo)

            wr = 0.1*wgyro(bo)

            kr = calc_kr(wr,wp,wg)
            p_res = calc_pR(wr,wp,wg,0.)
            call logspace(1e-24,3000.*abs(p_res),size(pperp),pperp,dlnpperp)

            pRtest = abs(calc_pR(wr,wp,wg,0.))
            p_res = abs(calc_pR(wr,wp,wg,0.1*pRtest))
            call calc_derivativesPA(pp_MK,sina_K,PSD_MK/pp_MK**2,&
                                    p_res,pRtest,dlnpperp,dlnpperp,&
                                    f,dfdpparR,dfdpperpR)
            write(*,*) "======== begin wr integration tests ========"
            write(*,*) 'ppar, pperp, f: ', p_res,pRtest,f
            write(*,*) 'dfdpperp test: ', dfdpperpR

            p0 = 1.5e-17*1e-5
            p_res = calc_pR(wr,wp,wg,p0)
            del_R = calc_delR(wr,wp,wg,p0)
            call calc_derivativesPA(pp_MK,sina_K,PSD_MK/pp_MK**2,&
                                    abs(p_res),p0,dlnpperp,dlnpperp,&
                                    f,dfdpparR,dfdpperpR)
            integrand = ( (p0)**2 / del_R ) * dfdpperpR
            write(*,*) 'Single integrand: ', integrand

            integral=0.0
            write(*,*) 'Integral1 Result: i, pperp, f, dpperp, pres, &
                                                       delR, dfdpperpR, INT1 '
            do ip=1,size(pperp)-1
                dpperp = (pperp(ip+1)-pperp(ip))
                p_res = calc_pR(wr,wp,wg,pperp(ip))
                del_R = calc_delR(wr,wp,wg,pperp(ip))
                gamma_L = calc_gamma_lorentz(wr,wp,wg,pperp(ip))
                call calc_derivativesPA(pp_MK,sina_K,PSD_MK/pp_MK**2,&
                                    abs(p_res),pperp(ip),dlnpperp,dlnpperp,&
                                    f,dfdpparR,dfdpperpR)
                integral = integral + &
                            (pperp(ip)**2/del_R) * dpperp * dfdpperpR / &
                            (ccm2m*ckeV2Joule)
                integral2 = integral2 + &
                            (pperp(ip)**2/(gamma_L * del_R)) * dpperp *&
                            (pperp(ip)*dfdpparR - abs(p_res)*dfdpperpR)&
                            / (ccm2m*ckeV2Joule)
                write(*,"(i6,7e11.3)") &
                    ip,pperp(ip),f,dpperp,p_res,del_R,dfdpperpR,integral
            enddo
            eta_r = (cPi * cMe * (wr - wgyro(bo))) / (N_dens * kr) * integral
            a_r = kr/(cMe*(wr-wg)) * integral2 / integral
            growth = calc_growth_rate(wr,wp,wg,eta_r,a_r)
            write(*,*) 'eta: ', eta_r
            write(*,*) 'a: ', a_r
            write(*,*) 'growth: ', growth
            write(*,*) "======== end wr integration tests ========"
        end subroutine integration_tests
        !-----------------------------------------------------------------
        subroutine write_test_psd(nsteps)
            use ModIoUnit, ONLY: UnitTmp_
            integer,intent(in) :: nsteps
            ! Write out PSD to file
            if (IsFirstTestCall) then
                open(unit=UnitTmp_,&
                    file='IM/plots/localwave_PSD.psd',&
                    status='unknown')
                ! Header info
                write(UnitTmp_,*)"! ntimes, ndim1, ndim2"
                write(UnitTmp_,"(3i6)") nsteps,nM,nK
                write(UnitTmp_,*)"! ppar (2D)"
                write(UnitTmp_,"(7e10.3)") ppar_MK
                write(UnitTmp_,*)"! pperp (2D)"
                write(UnitTmp_,"(7e10.3)") pperp_MK
            else
                open(unit=UnitTmp_,&
                    file='IM/plots/localwave_PSD.psd',&
                    status='old',&
                    position='append')
            endif
            ! Entry info
            write(UnitTmp_,*)"! tSim"
            write(UnitTmp_,"(f6.1)") tSim
            do iK=1,nK
                write(UnitTmp_,"(8e9.2)") PSD_MK(:,iK)
            enddo
            close(UnitTmp_)
            if (IsFirstTestCall) write(*,*)"SAVED: IM/plots/localwave_PSD.psd"
        end subroutine
        !-----------------------------------------------------------------
        subroutine write_test_fin
            use ModIoUnit, ONLY: UnitTmp_
            ! Write out a .fin file for this synthetic distribution
            open(unit=UnitTmp_,&
                 file='IM/plots/synthetic.fin',&
               status='unknown')
            ! Header info
            write(UnitTmp_,"(2i6,a)")1,nM," ! il, ie"
            write(UnitTmp_,"(1i6,a)")1,&
                             " ! Synthetic flux; iunit: 1=(cm^2 s sr keV)^-1"
            write(UnitTmp_,"(1f7.2,a)") 4.50, " ! ro"
            write(UnitTmp_,"(7f11.3)") Ekev_MK(:,nK)
            write(UnitTmp_,"(7e10.3)") PSD_MK(:,nK)
            write(*,*)"SAVED: IM/plots/synthetic.fin"
            close(UnitTmp_)
        end subroutine write_test_fin
    end subroutine unit_test_localwave
    !-------------------------------------------------------------------------

end module CIMI_LocalWave

