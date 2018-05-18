subroutine CIMI_set_parameters(NameAction)

  use ModIoUnit,         ONLY: UnitTmp_, io_unit_new
  use ModReadParam
  use ModUtilities,	 ONLY: lower_case
  use ModCimiInitialize, ONLY: &
       IsEmptyInitial, IsDataInitial, IsRBSPData, IsGmInitial
  use ModCimiPlot
  use ModCimiTrace,	 ONLY: UseEllipse, UseSmooth, UseCorotation, &
       UsePotential, SmoothWindow, imod, iLatTest, iLonTest
  use ModCimi,           ONLY: UseMcLimiter, BetaLimiter, time, Pmin, &
       IsStandAlone, UseStrongDiff, UseDecay, DecayTimescale,&
       dt, dtmax, DoCalcPrecip, DtCalcPrecip
  use ModCimiRestart,    ONLY: IsRestart, DtSaveRestart
  use ModCimiPlanet,     ONLY: nspec, NameSpecies_I
  use ModImTime,         ONLY: iStartTime_I, TimeMax
  use ModCimiBoundary,   ONLY: &
       UseBoundaryEbihara, UseYoungEtAl, CIMIboundary, Outputboundary
  use ModIeCimi,         ONLY: UseWeimer
  use ModPrerunField,    ONLY: DoWritePrerun, UsePrerun, DtRead
  use ModGmCIMI,         ONLY: UseGm
  use ModWaveDiff,       ONLY: &
       UseWaveDiffusion, UseHiss, UseChorus, UseChorusUB, &
       DiffStartT, HissWavesD, ChorusWavesD, ChorusUpperBandD, &
       testDiff_aa, testDiff_EE, testDiff_aE, &
       NameAeFile, read_ae_wdc_kyoto
  use DensityTemp, 	 ONLY: densityP
  use ModImSat,          ONLY: DtSatOut, DoWritePrerunSat, UsePrerunSat, &
       DtReadSat, DoWriteSats, ReadRestartSat
  use ModCimiGrid
  
  implicit none

  logical 			:: DoEcho = .false.
  character (len=100)           :: NameCommand, StringCIMIPlot
  character (len=*), intent(in) :: NameAction
  character (len=7)             :: TypeBoundary
  character (len=3)             :: NameModel
  character (len=*), parameter  :: NameSub = 'CIMI_set_parameters'

  character (len=100) :: cTempLine
  character (len=100), dimension(100) :: cTempLines
  
  integer :: iError, iDate, iCIMIPlotType, nCIMIPlotType
  real :: DensitySW, VelSW, BxSW, BySW, BzSW, DtOutputCIMIPlot
  logical :: DoSaveSeparateFiles
  !\
  ! Description:
  ! This subroutine gets the inputs for CIMI
  !/

  !---------------------------------------------------------------------------
  
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE

     select case(NameCommand)

     case("#ECHO")
!        call check_stand_alone
        call read_var('DoEcho', DoEcho)
        if(iProc==0)call read_echo_set(DoEcho)

     case('#STOP')
        if(IsStandAlone)then
           call read_var('TimeMax',TimeMax)
        else
           write(*,*)'IM WARNING: #STOP command is ignored in the framework'
        end if

     case('#STARTTIME')
        if(IsStandAlone)then
           !read in iYear,iMonth,iDay,iHour,iMinute,iSecond into iStartTime
           do iDate=1,6
              call read_var('iStartTime',iStartTime_I(iDate))
           enddo
           
        else
           write(*,*)'IM WARNING:#STARTTIME command is ignored in the framework'
        end if

     case ("#NGDC_INDICES")
        cTempLines(1) = NameCommand
        call read_var('NameNgdcFile', cTempLine)
        cTempLines(2) = cTempLine
        cTempLines(3) = " "
        cTempLines(4) = "#END"
        
        call IO_set_inputs(cTempLines)
        call read_NGDC_Indices(iError)
        
        if (iError /= 0) then 
           write(*,*) "read indices was NOT successful (NOAA file)"
        endif

     case ("#MHD_INDICES")
        cTempLines(1) = NameCommand
        call read_var('UpstreamFile',cTempLine)
        cTempLines(2) = cTempLine
        cTempLines(3) = " "
        cTempLines(4) = "#END"
        
        call IO_set_inputs(cTempLines)
        call read_MHDIMF_Indices(iError)

     case ("#SOLARWIND")
        call read_var('n',DensitySW)
        call read_var('vx',VelSW)
        call read_var('bx',BxSW)
        call read_var('by',BySW)
        call read_var('bz',BzSW)
        
        call IO_set_imf_bx_single(BxSW)
        call IO_set_imf_by_single(BySW)
        call IO_set_imf_bz_single(BzSW)
        call IO_set_sw_v_single(abs(VelSW))
        call IO_set_sw_n_single(abs(DensitySW))

     case('#SMOOTH')
        !do you want to smooth the solar wind for Tsy model
        call read_var('UseSmooth',UseSmooth)
        call read_var('SmoothWindow',SmoothWindow)

     case('#BMODEL')
        call read_var('NameModel',NameModel)!t96,t04,MHD,Dip
        !call read_var('UseFixedB',UseFixedB)!T=fixed B config or 
                                            !F=changing B config
        if (NameModel == 'Dip') then
           iMod=0
           ! Checks the next two lines to see if cor and/or pot
           ! variables should be set to 0.
           call read_var('UseCorotation',UseCorotation)
           call read_var('UsePotential',UsePotential)
        elseif (NameModel == 't96') then
           iMod=1
        elseif(NameModel == 't04') then
           iMod=2
        elseif(NameModel == 'MHD')then
           iMod=3
           UseGm=.true.
        else
           call con_stop('Error: Model not found') 
        endif

     case('#IEMODEL')
        call read_var('UseWeimer',UseWeimer)

!!$     case('#SAVEPLOT')
!!$        call read_var('DtSavePlot',DtOutput)
!!$        call read_var('DoSaveFlux',DoSaveFlux)
!!$        call read_var('DoSaveDrifts',DoSaveDrifts)
!!$        call read_var('DoSavePSD',DoSavePSD)
!!$        ! If saving flux then decide if it should be just one file or many
!!$        if (DoSaveFlux .or. DoSaveDrifts .or. DoSavePSD) then
!!$           call read_var('UseSeparatePlotFiles',UseSeparatePlotFiles)
!!$        endif
!!$        
!!$        DoSavePlot = .true.

     case('#SAVELOG')
        call read_var('DtLogOut',DtLogOut)
        DoSaveLog = .true.

     case('#SAVEPLOT')

        call read_var( 'nCIMIPlotType', nCIMIPlotType )

        CIMI_PLOTTYPE: do iCIMIPlotType = 1, nCIMIPlotType

           call read_var( 'StringPlot', StringCIMIPlot )
           call lower_case( StringCIMIPlot )
           call read_var( 'DtOutput', DtOutputCIMIPlot )
           
           if ( index( StringCIMIPlot, 'fls'  ) > 0 .or. &
                index( StringCIMIPlot, 'flux' ) > 0 ) then

              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )

              if (iProc == 0 ) then
                 write(*,*) "DoSaveFlux: ", DoSaveFlux
                 write(*,*) "DtFluxOutput: ", DtFluxOutput
                 write(*,*) "DoSaveSeparateFluxFiles: ", &
                      DoSaveSeparateFluxFiles
              endif
              
              PLOT_FLUX_SPECIES: if &
                   ( index( StringCIMIPlot, 'all' ) > 0 ) then
                 
                 DoSaveFlux( 1 : nspec ) = .true.
                 DtFluxOutput( 1 : nspec ) = DtOutputCIMIPlot
                 DoSaveSeparateFluxFiles( 1 : nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'ions' ) > 0 ) then
                 
                 DoSaveFlux( 1 : nspec - 1 ) = .true.
                 DtFluxOutput( 1 : nspec - 1 ) = DtOutputCIMIPlot
                 DoSaveSeparateFluxFiles( 1 : nspec - 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'electrons' ) > 0 .or. &
                     index( StringCIMIPlot, 'e-' ) > 0 ) then
                 
                 DoSaveFlux( nspec ) = .true.
                 DtFluxOutput( nspec ) = DtOutputCIMIPlot
                 DoSaveSeparateFluxFiles( nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'h+' ) > 0 ) then
                 
                 DoSaveFlux( 1 ) = .true.
                 DtFluxOutput( 1 ) = DtOutputCIMIPlot
                 DoSaveSeparateFluxFiles( 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'o+' ) > 0 ) then
                 
                 if ( nspec .ge. 3 ) then
                    
                    DoSaveFlux( 2 ) = .true. 
                    DtFluxOutput( 2 ) = DtOutputCIMIPlot
                    DoSaveSeparateFluxFiles( 2 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_STOP( 'O+ not configured; Recompile CIMI '// &
                         'with ModEarthHO or ModEarthHOHe options.' )
                    
                 endif

              elseif &
                   ( index( StringCIMIPlot, 'he+' ) > 0 ) then

                 if ( nspec .gt. 3 ) then
                    
                    DoSaveFlux( 3 ) = .true. 
                    DtFluxOutput( 3 ) = DtOutputCIMIPlot
                    DoSaveSeparateFluxFiles( 3 ) = &
                         DoSaveSeparateFiles

                 else

                    call CON_STOP( 'He+ not configured; Recompile CIMI '//&
                         'with ModEarthHOHe option.' )

                 endif

              else

                 call CON_STOP( 'No flux species information; STOPPING' )

              endif PLOT_FLUX_SPECIES
              
              if (iProc == 0 ) then
                 write(*,*) "DoSaveFlux: ", DoSaveFlux
                 write(*,*) "DtFluxOutput: ", DtFluxOutput
                 write(*,*) "DoSaveSeparateFluxFiles: ", &
                      DoSaveSeparateFluxFiles
              endif
              
           elseif &
                ( index( StringCIMIPlot, 'psd' ) > 0 ) then

              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )

              if (iProc == 0 ) then
                 write(*,*) "DoSavePSD: ", DoSavePSD
                 write(*,*) "DtPSDOutput: ", DtPSDOutput
                 write(*,*) "DoSaveSeparatePSDFiles: ", &
                      DoSaveSeparatePSDFiles
              endif
              
              PLOT_PSD_SPECIES: if &
                   ( index( StringCIMIPlot, 'all' ) > 0 ) then

                 DoSavePSD( 1 : nspec ) = .true.
                 DtPSDOutput( 1 : nspec ) = DtOutputCIMIPlot
                 DoSaveSeparatePSDFiles( 1 : nspec ) = &
                      DoSaveSeparateFiles

              elseif &
                   ( index( StringCIMIPlot, 'ions' ) > 0 ) then

                 DoSavePSD( 1 : nspec - 1 ) = .true.
                 DtPSDOutput( 1 : nspec - 1 ) = DtOutputCIMIPlot
                 DoSaveSeparatePSDFiles( 1 : nspec - 1 ) = &
                      DoSaveSeparateFiles

              elseif &

                   ( index( StringCIMIPlot, 'electrons' ) > 0 .or. &
                     index( StringCIMIPlot, 'e-' ) > 0 ) then
                 DoSavePSD( nspec ) = .true.
                 DtPSDOutput( nspec ) = DtOutputCIMIPlot
                 DoSaveSeparatePSDFiles( nspec ) = &
                      DoSaveSeparateFiles

              elseif &
                   ( index( StringCIMIPlot, 'h+' ) > 0 ) then

                 DoSavePSD( 1 ) = .true.
                 DtPSDOutput( 1 ) = DtOutputCIMIPlot
                 DoSaveSeparatePSDFiles( 1 ) = &
                      DoSaveSeparateFiles

              elseif &
                   ( index( StringCIMIPlot, 'o+' ) > 0 ) then

                 if ( nspec .ge. 3 ) then
                    DoSavePSD( 2 ) = .true. 
                    DtPSDOutput( 2 ) = DtOutputCIMIPlot
                    DoSaveSeparatePSDFiles( 2 ) = &
                         DoSaveSeparateFiles
                 else
                    
                    call CON_STOP( 'O+ not configured; Recompile CIMI '// &
                         'with ModEarthHO or ModEarthHOHe options.' )
                    
                 endif
                 
              elseif &
                   ( index( StringCIMIPlot, 'he+' ) > 0 ) then
                 
                 if ( nspec .gt. 3 ) then
                    
                    DoSavePSD( 3 ) = .true. 
                    DtPSDOutput( 3 ) = DtOutputCIMIPlot
                    DoSaveSeparatePSDFiles( 3 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_STOP( 'He+ not configured; Recompile CIMI '//&
                         'with ModEarthHOHe option.' )
                    
                 endif
                 
              else
                 
                 call CON_STOP( 'No PSD species information; STOPPING' )
                 
              endif PLOT_PSD_SPECIES
              
              if (iProc == 0 ) then
                 write(*,*) "DoSavePSD: ", DoSavePSD
                 write(*,*) "DtPSDOutput: ", DtPSDOutput
                 write(*,*) "DoSaveSeparatePSDFiles: ", &
                      DoSaveSeparatePSDFiles
              endif
              
           elseif &
                ( index( StringCIMIPlot, 'vl'  ) > 0 .or. &
                  index( StringCIMIPlot, 'vldrift' ) > 0 ) then

              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )
              
              if (iProc == 0 ) then
                 write(*,*) "DoSaveVLDrift: ", DoSaveVLDrift
                 write(*,*) "DtVLDriftOutput: ", DtVLDriftOutput
                 write(*,*) "DoSaveSeparateVLDriftFiles: ", &
                      DoSaveSeparateVLDriftFiles
              endif
              
              PLOT_VLDRIFT_SPECIES: if &
                   ( index( StringCIMIPlot, 'all' ) > 0 ) then
                 
                 DoSaveVLDrift( 1 : nspec ) = .true.
                 DtVLDriftOutput( 1 : nspec ) = DtOutputCIMIPlot
                 DoSaveSeparateVLDriftFiles( 1 : nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'ions' ) > 0 ) then
                 
                 DoSaveVLDrift( 1 : nspec - 1 ) = .true.
                 DtVLDriftOutput( 1 : nspec - 1 ) = DtOutputCIMIPlot
                 DoSaveSeparateVLDriftFiles( 1 : nspec - 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'electrons' ) > 0 .or. &
                     index( StringCIMIPlot, 'e-' ) > 0 ) then
                 
                 DoSaveVLDrift( nspec ) = .true.
                 DtVLDriftOutput( nspec ) = DtOutputCIMIPlot
                 DoSaveSeparateVLDriftFiles( nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'h+' ) > 0 ) then
                 
                 DoSaveVLDrift( 1 ) = .true.
                 DtVLDriftOutput( 1 ) = DtOutputCIMIPlot
                 DoSaveSeparateVLDriftFiles( 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'o+' ) > 0 ) then
                 
                 if ( nspec .ge. 3 ) then
                    
                    DoSaveVLDrift( 2 ) = .true. 
                    DtVLDriftOutput( 2 ) = DtOutputCIMIPlot
                    DoSaveSeparateVLDriftFiles( 2 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_STOP( 'O+ not configured; Recompile CIMI '// &
                         'with ModEarthHO or ModEarthHOHe options.' )
                    
                 endif
              elseif &
                   ( index( StringCIMIPlot, 'he+' ) > 0 ) then
                 
                 if ( nspec .gt. 3 ) then
                    
                    DoSaveVLDrift( 3 ) = .true. 
                    DtVLDriftOutput( 3 ) = DtOutputCIMIPlot
                    DoSaveSeparateVLDriftFiles( 3 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_STOP( 'He+ not configured; Recompile CIMI '//&
                         'with ModEarthHOHe option.' )
                    
                 endif
                 
              else
                 
                 call CON_STOP( 'No VLDrift species information; STOPPING' )
                 
              endif PLOT_VLDRIFT_SPECIES
              
              if (iProc == 0 ) then
                 write(*,*) "DoSaveVLDrift: ", DoSaveVLDrift
                 write(*,*) "DtVLDriftOutput: ", DtVLDriftOutput
                 write(*,*) "DoSaveSeparateVLDriftFiles: ", &
                      DoSaveSeparateVLDriftFiles
              endif
              
           elseif &
                ( index( StringCIMIPlot, 'vp'  ) > 0 .or. &
                  index( StringCIMIPlot, 'vpdrift' ) > 0 ) then

              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )

              if (iProc == 0 ) then
                 write(*,*) "DoSaveVPDrift: ", DoSaveVPDrift
                 write(*,*) "DtVPDriftOutput: ", DtVPDriftOutput
                 write(*,*) "DoSaveSeparateVPDriftFiles: ", &
                      DoSaveSeparateVPDriftFiles
              endif
              
              PLOT_VPDRIFT_SPECIES: if &
                   ( index( StringCIMIPlot, 'all' ) > 0 ) then
                 
                 DoSaveVPDrift( 1 : nspec ) = .true.
                 DtVPDriftOutput( 1 : nspec ) = DtOutputCIMIPlot
                 DoSaveSeparateVPDriftFiles( 1 : nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'ions' ) > 0 ) then
                 
                 DoSaveVPDrift( 1 : nspec - 1 ) = .true.
                 DtVPDriftOutput( 1 : nspec - 1 ) = DtOutputCIMIPlot
                 DoSaveSeparateVPDriftFiles( 1 : nspec - 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'electrons' ) > 0 .or. &
                     index( StringCIMIPlot, 'e-' ) > 0 ) then

                 DoSaveVPDrift( nspec ) = .true.
                 DtVPDriftOutput( nspec ) = DtOutputCIMIPlot
                 DoSaveSeparateVPDriftFiles( nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'h+' ) > 0 ) then
                 
                 DoSaveVPDrift( 1 ) = .true.
                 DtVPDriftOutput( 1 ) = DtOutputCIMIPlot
                 DoSaveSeparateVPDriftFiles( 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'o+' ) > 0 ) then
                 
                 if ( nspec .ge. 3 ) then
                    
                    DoSaveVPDrift( 2 ) = .true. 
                    DtVPDriftOutput( 2 ) = DtOutputCIMIPlot
                    DoSaveSeparateVPDriftFiles( 2 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_STOP( 'O+ not configured; Recompile CIMI '// &
                         'with ModEarthHO or ModEarthHOHe options.' )
                    
                 endif
                 
              elseif &
                   ( index( StringCIMIPlot, 'he+' ) > 0 ) then
                 
                 if ( nspec .gt. 3 ) then
                    
                    DoSaveVPDrift( 3 ) = .true. 
                    DtVPDriftOutput( 3 ) = DtOutputCIMIPlot
                    DoSaveSeparateVPDriftFiles( 3 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_STOP( 'He+ not configured; Recompile CIMI '//&
                         'with ModEarthHOHe option.' )
                    
                 endif
                 
              else
                 
                 call CON_STOP( 'No VPDrift species information; STOPPING' )
                 
              endif PLOT_VPDRIFT_SPECIES
              
              if (iProc == 0 ) then
                 write(*,*) "DoSaveVPDrift: ", DoSaveVPDrift
                 write(*,*) "DtVPDriftOutput: ", DtVPDriftOutput
                 write(*,*) "DoSaveSeparateVPDriftFiles: ", &
                      DoSaveSeparateVPDriftFiles
              endif
              
           elseif &
                ( index( StringCIMIPlot, 'precipitation'  ) > 0 .or. &
                  index( StringCIMIPlot, 'precip' ) > 0 .or. &
                  index( StringCIMIPlot, 'preci' ) > 0 ) then

              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )
              
              if (iProc == 0 ) then
                 write(*,*) "DoSavePreci: ", DoSavePreci
                 write(*,*) "DtPreciOutput: ", DtPreciOutput
                 write(*,*) "DoSaveSeparatePreciFiles: ", &
                      DoSaveSeparatePreciFiles
              endif
              
              PLOT_PRECI_SPECIES: if &
                   ( index( StringCIMIPlot, 'all' ) > 0 ) then
                 
                 DoSavePreci( 1 : nspec ) = .true.
                 DtPreciOutput( 1 : nspec ) = DtOutputCIMIPlot
                 DoSaveSeparatePreciFiles( 1 : nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'ions' ) > 0 ) then
                 
                 DoSavePreci( 1 : nspec - 1 ) = .true.
                 DtPreciOutput( 1 : nspec - 1 ) = DtOutputCIMIPlot
                 DoSaveSeparatePreciFiles( 1 : nspec - 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'electrons' ) > 0 .or. &
                     index( StringCIMIPlot, 'e-' ) > 0 ) then

                 DoSavePreci( nspec ) = .true.
                 DtPreciOutput( nspec ) = DtOutputCIMIPlot
                 DoSaveSeparatePreciFiles( nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'h+' ) > 0 ) then
                 
                 DoSavePreci( 1 ) = .true.
                 DtPreciOutput( 1 ) = DtOutputCIMIPlot
                 DoSaveSeparatePreciFiles( 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'o+' ) > 0 ) then
                 
                 if ( nspec .ge. 3 ) then
                    
                    DoSavePreci( 2 ) = .true. 
                    DtPreciOutput( 2 ) = DtOutputCIMIPlot
                    DoSaveSeparatePreciFiles( 2 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_STOP( 'O+ not configured; Recompile CIMI '// &
                         'with ModEarthHO or ModEarthHOHe options.' )
                    
                 endif
                 
              elseif &
                   ( index( StringCIMIPlot, 'he+' ) > 0 ) then
                 
                 if ( nspec .gt. 3 ) then
                    
                    DoSavePreci( 3 ) = .true. 
                    DtPreciOutput( 3 ) = DtOutputCIMIPlot
                    DoSaveSeparatePreciFiles( 3 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_STOP( 'He+ not configured; Recompile CIMI '//&
                         'with ModEarthHOHe option.' )
                    
                 endif
                 
              else
                 
                 call CON_STOP( 'No Preci species information; STOPPING' )
                 
              endif PLOT_PRECI_SPECIES
              
              if (iProc == 0 ) then
                 write(*,*) "DoSavePreci: ", DoSavePreci
                 write(*,*) "DtPreciOutput: ", DtPreciOutput
                 write(*,*) "DoSaveSeparatePreciFiles: ", &
                      DoSaveSeparatePreciFiles
              endif

           elseif &
                ( index( StringCIMIPlot, '2d'  ) > 0 ) then

              if (iProc == 0 ) then
                 write(*,*) "DoSaveEq: ", DoSaveEq
                 write(*,*) "DoSaveIono: ", DoSaveIono
                 write(*,*) "DtOutput: ", DtOutput
              endif

              PLOT_2D: if &
                   ( index( StringCIMIPlot, 'all'  ) > 0 .or. &
                     index( StringCIMIPlot, 'both'  ) > 0 ) then

                 DoSaveEq = .true.
                 DoSaveIono = .true.
                 DtOutput = DtOutputCIMIPlot
              elseif &
                   ( index( StringCIMIPlot, 'equator' ) > 0 .or. &
                     index( StringCIMIPlot, 'eq' ) > 0 ) then

                 DoSaveEq = .true.
                 DtOutput = DtOutputCIMIPlot

              elseif &
                   ( index( StringCIMIPlot, 'ionosphere' ) > 0 .or. &
                     index( StringCIMIPlot, 'iono' ) > 0 ) then

                 DoSaveIono = .true.
                 DtOutput = DtOutputCIMIPlot
                 
              else

                 call CON_STOP( 'No CIMI 2D Plot information; STOPPING' )
                 
              endif PLOT_2D

              if (iProc == 0 ) then
                 write(*,*) "DoSaveEq: ", DoSaveEq
                 write(*,*) "DoSaveIono: ", DoSaveIono
                 write(*,*) "DtOutput: ", DtOutput
              endif
              
           elseif &
                ( index( StringCIMIPlot, 'lstar'  ) > 0 .or. &
                  index( StringCIMIPlot, 'l*'  ) > 0 ) then

              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )

              if (iProc == 0 ) then
                 write(*,*) "DoSaveLstar: ", DoSaveLstar
                 write(*,*) "DtLstarOutput: ", DtLstarOutput
                 write(*,*) "DoSaveSeparateLstarFiles: ", &
                      DoSaveSeparateLstarFiles
              endif
              
              DoSaveLstar = .true.
              DtLstarOutput = DtOutputCIMIPlot
              DoSaveSeparateLstarFiles = DoSaveSeparateFiles

              if (iProc == 0 ) then
                 write(*,*) "DoSaveLstar: ", DoSaveLstar
                 write(*,*) "DtLstarOutput: ", DtLstarOutput
                 write(*,*) "DoSaveSeparateLstarFiles: ", &
                      DoSaveSeparateLstarFiles
              endif
              
           else
              call CON_STOP( 'No plot type specified.' )
           endif

        end do CIMI_PLOTTYPE

        DoSavePlot = .true.
     
     case('#INITIALF2')
        call read_var('IsEmptyInitial',IsEmptyInitial)
        call read_var('IsGmInitial',   IsGmInitial)
        call read_var('IsDataInitial', IsDataInitial)
        if (IsDataInitial) &
             call read_var('IsRBSPData', IsRBSPData)
        
        !IsDataInitial only works with EarthHO or EarthH configurations
        if (nspec > 3 .and. IsDataInitial) &
             call CON_STOP('IsDataInitial only works with EarthHO or EarthH')  

     case('#TYPEBOUNDARY')
        call read_var('TypeBoundary',TypeBoundary)
        if(TypeBoundary == 'Ellipse') then
           UseEllipse = .true.
        else
           UseEllipse = .false.
        endif
        
     case('#PLASMASHEET')
        call read_var('UseYoungEtAl',UseYoungEtAl)
        if(IsStandAlone)then
           ! When Standalone specify boundary type
           call read_var('UseBoundaryEbihara',UseBoundaryEbihara)
        endif
        
     case('#RESTART')
        call read_var('IsRestart',IsRestart) !T:Continuous run, F: Initial Run
        if(IsRestart) call read_var('ReadRestartSat',ReadRestartSat)
        
     case('#SAVERESTART')
        ! when in standalone mode read restart save frequency
        if(IsStandAlone) then
           call read_var('DtSaveRestart',DtSaveRestart)
        else
           call CON_STOP('ERROR: Restart save frequency only set in standalone')
        end if

     case('#LIMITER')
        call read_var('UseMcLimiter', UseMcLimiter)
        if(UseMcLimiter) call read_var('BetaLimiter', BetaLimiter)

     ! minimum pressure in nPa passed to GM
     case('#MINIMUMPRESSURETOGM')   
        call read_var('MinimumPressureToGM', Pmin)
        
     case('#TIMESIMULATION')
        call read_var('TimeSimulation',time)

     case('#PRERUNSAT')
        call read_var('DoWritePrerunSat',DoWritePrerunSat)
        if(.not.DoWritePrerunSat) call read_var('UsePrerunSat',UsePrerunSat)
        if(UsePrerunSat) then
           call read_var('DtReadSat',   DtReadSat)
           DtSatout = DtReadSat
           DoWriteSats = .true.
        end if
        
     case('#PRERUNFIELD')
        call read_var('DoWritePrerun',DoWritePrerun)
        if(.not.DoWritePrerun) call read_var('UsePrerun',UsePrerun)
        if(UsePrerun)          call read_var('DtRead',   DtRead)
        if(UsePrerun .or. DoWritePrerun) UseGm=.true.

     case('#DECAY')
        call read_var('UseDecay',UseDecay)
        if ( UseDecay ) &
             call read_var('DecayTimescale in seconds', DecayTimescale)
        
     case('#STRONGDIFFUSION')
        call read_var('UseStrongDiff',UseStrongDiff)
        
     case('#WAVEDIFFUSION')
        call read_var('UseWaveDiffusion',UseWaveDiffusion)
        if(UseWaveDiffusion) then
           call read_var('DiffStartT',  DiffStartT)
           call read_var('UseHiss',  UseHiss)
           call read_var('UseChorus',  UseChorus)
           call read_var('UseChorusUB',  UseChorusUB)
           call read_var('HissWavesD', HissWavesD)
           call read_var('ChorusWavesD', ChorusWavesD)
           call read_var('ChorusUpperBandD', ChorusUpperBandD)
           call read_var('NameAeFile',NameAeFile)
           call read_ae_wdc_kyoto(iError)

           if (iError /= 0) then
              write(*,*) "read AE index was not successful "//&
                   "(WDC Kyoto file; IAGA2002 Format)"
           end if

        end if

     case('#ENERGYGRID')
        call read_var('MinIonEnergy (in keV), MinElectronEnergy x10',&
             MinIonEnergy)
        call read_var('MaxIonEnergy (in keV), MaxElectronEnergy x10',&
             MaxIonEnergy)
        
     case('#SETENERGYGRID')
        call read_var('neng', neng)
        if (neng .lt. 2) &
             call CON_STOP('ERROR: neng must be greater than 1.')
        call read_var('UseLogEGrid', UseLogEGrid)
        call read_var('MinIonEnergy (in keV), MinElectronEnergy x10',&
             MinIonEnergy)
        call read_var('MaxIonEnergy (in keV), MaxElectronEnergy x10',&
             MaxIonEnergy)

     case('#RBSPENERGYGRID')
        call read_var('UseRBSPGrid', UseRBSPGrid)
        if (UseRBSPGrid) neng = 30
        
     case('#IMTIMESTEP')
        call read_var('IMDeltaT [s]',dt)
        call read_var('IMDeltaTMax [s]',dtmax)

     case('#TESTFIELDLINE')
        call read_var('iLatTest',iLatTest)
        call read_var('iLonTest',iLonTest)
        
     case('#PLASMAPAUSEDENSITY')
        call read_var('DensityP [m^3]',densityP)
        
     case('#SETRB')
        call read_var('rb [R_E]', rb)

     case('#BOUNDARYCHECK')
        call read_var('CIMIboundary',CIMIboundary)
        if (CIMIboundary) call read_var('Outputboundary',Outputboundary)

     case('#PRECIPITATION')
        call read_var('DoCalcPrecip',DoCalcPrecip)
        if (DoCalcPrecip) call read_var('DtCalcPrecip',DtCalcPrecip)
!!$        if (DoCalcPrecip) call read_var('PrecipOutput',PrecipOutput)
!!$        if (PrecipOutput) call read_var('DtPreOut',DtPreOut)

     end select
     
  enddo

end subroutine CIMI_set_parameters
