module CIMI_SetParams
    use ModReadParam, ONLY: read_line, read_command, read_var, read_echo_set

    implicit none

    private ! except

    public:: CIMI_set_parameters
    character (len=*), parameter  :: NameSub = 'CIMI_set_parameters'


contains
subroutine CIMI_set_parameters(NameAction)

  ! Description:
  ! This subroutine gets the inputs for CIMI
  
  character (len=*), intent(in) :: NameAction
  character (len=100)           :: NameCommand

  !---------------------------------------------------------------------------
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE

     select case(NameCommand)

     case("#ECHO")
        call ECHO_command()
     case("#STOP")
        call STOP_command()
     case("#STARTTIME")
        call STARTTIME_command()
     case ("#NGDC_INDICES")
        call NGDC_INDICES_command()
     case ("#MHD_INDICES")
        call MHD_INDICES_command()
     case ("#KYOTO_DST")
        call KYOTO_DST_command()
     case ("#KYOTO_AE")
        call KYOTO_AE_command()
     case ("#POTSDAM_KP_AP_F107")
        call POTSDAM_KP_AP_F107_command()
     case ("#SOLARWIND")
        call SOLARWIND_command()
     case("#SINGLEGRID")
        call SINGLEGRID_command()
     case("#SMOOTH")
        call SMOOTH_command()
     case("#BMODEL")
        call BMODEL_command()
     case("#IEMODEL")
        call IEMODEL_command()
     case("#GIMME")
        call GIMME_command() 
     case("#SAVELOG")
        call SAVELOG_command()
     case("#SAVEPLOT")
        call SAVEPLOT_command()
     case('#DTSATOUT')
        call DTSATOUT_command()
     case('#VERBOSELSTAR')
        call VERBOSELSTAR_command()
     case('#INITIALF2')
        call INITIALF2_command()
     case('#INITIALLSTAR')
        call INITIALLSTAR_command()
     case('#TYPEBOUNDARY')
        call TYPEBOUNDARY_command()
     case('#PLASMASHEET')
        call PLASMASHEET_command()
     case('#BOUNDARYKAPPA')
        call BOUNDARYKAPPA_command()
     case('#RESTART')
        call RESTART_command()
     case('#SAVERESTART')
        call SAVERESTART_command()
     case('#LIMITER')
        call LIMITER_command()
     case('#MINIMUMPRESSURETOGM')   
        call MINIMUMPRESSURETOGM_command()
     case('#ELECTRONTEMPERATURERATIO')
        call ELECTRONTEMPERATURERATIO_command()
     case('#GMTEMPERATUREFACTORS')
        call GMTEMPERATUREFACTORS_command()
     case('#TIMESIMULATION')
        call TIMESIMULATION_command()
     case('#PRERUNSAT')
        call PRERUNSAT_command()
     case('#PRERUNFIELD')
        call PRERUNFIELD_command()
     case('#DECAY')
        call DECAY_command()
     case('#ELECTRONDECAY')
        call ELECTRONDECAY_command()
     case('#FLC')
        call FLC_command()
     case('#COMPOSITION')
        call COMPOSITION_command()
     case('#STRONGDIFFUSION')
        call STRONGDIFFUSION_command()
     case('#WAVES')
        call WAVES_command()
     case('#DIAGONALIZEDDIFFUSION')
        call DIAGONALIZEDDIFFUSION_command()
     case('#DIAGDIFFUSIONTEST')
        call DIAGDIFFUSIONTEST_command()
     case('#ENERGYGRID')
        call ENERGYGRID_command()
     case('#SETENERGYGRID')
        call SETENERGYGRID_command()
     case('#RBSPENERGYGRID')
        call RBSPENERGYGRID_command()
     case('#ELECTRONGRIDRATIO')
        call ELECTRONGRIDRATIO_command()
     case('#IMTIMESTEP')
        call IMTIMESTEP_command()
     case('#TESTFIELDLINE')
        call TESTFIELDLINE_command()
     case('#PLASMAPAUSE')
        call PLASMAPAUSE_command()
     case('#PLASMAPAUSEDENSITY')
        call PLASMAPAUSEDENSITY_command()
     case('#COREPLASMASPHERE')
        call COREPLASMASPHERE_command()
     case('#SETRB')
        call SETRB_command()
     case('#SETBOUNDARYPARAMS')
        call SETBOUNDARYPARAMS_command()
     case('#PRECIPITATION')
        call PRECIPITATION_command()
     case('#STRICTDRIFT')
        call STRICTDRIFT_command()
     case('#LATITUDINALGRID')
        call LATITUDINALGRID_command()
     case('#VERBOSELATGRID')
        call VERBOSELATGRID_command()
     case('#HIGHERORDERDRIFT')   ! use higher order drift
        call HIGHERORDERDRIFT_command()
     case('#DRIFTSCHEME')   ! use higher order drift
        call DRIFTSCHEME_command()
     end select
     
  enddo

    contains
    !-------------------------------------------------------------------------
    subroutine ECHO_command
        use ModCimiGrid, ONLY: iProc

        logical :: DoEcho = .false.

        ! call check_stand_alone
        call read_var('DoEcho', DoEcho)
        if(iProc==0)call read_echo_set(DoEcho)
    end subroutine ECHO_command
    !-------------------------------------------------------------------------
    subroutine STOP_command
        use ModCimi,		 ONLY: IsStandAlone
        use ModImTime,	     ONLY: TimeMax

        if(IsStandAlone)then
            call read_var('TimeMax',TimeMax)
        else
            write(*,*)'IM WARNING: #STOP command is ignored in the framework'
        end if
    end subroutine STOP_command
    !-------------------------------------------------------------------------
    subroutine STARTTIME_command
        use ModTimeConvert,  ONLY: time_int_to_real
        use ModCimi,		 ONLY: IsStandAlone
        use ModImTime,	     ONLY: iStartTime_I, StartTime
        use CON_axes,        ONLY: init_axes
        use ModPlanetConst,  ONLY: iPlanet

        integer :: iDate

        if(IsStandAlone)then
            !read in iYear,iMonth,iDay,iHour,iMinute,iSecond into iStartTime
            do iDate=1,6
                call read_var('iStartTime',iStartTime_I(iDate))
            enddo
            call time_int_to_real(iStartTime_I, StartTime)
            call init_axes(StartTime)

        else
            write(*,*)'IM WARNING:#STARTTIME command is ignored in the &
                        framework'
        end if
    end subroutine STARTTIME_command
    !-------------------------------------------------------------------------
    subroutine NGDC_INDICES_command
        use ModUtilities, ONLY: CON_stop

        character (len=100) :: cTempLine
        character (len=100), dimension(100) :: cTempLines
        integer:: iError

        cTempLines(1) = '#NGDC_INDICES'
        call read_var('NameNGDCFile', cTempLine)
        cTempLines(2) = cTempLine
        cTempLines(3) = " "
        cTempLines(4) = "#END"

        call IO_set_inputs(cTempLines)
        call read_NGDC_Indices(iError)

        if (iError /= 0) then
            write(*,*) 'IM: read_NGDC_Indices iError=', iError
            call CON_stop('IM: Could not read '//trim(cTempLine))
        endif
    end subroutine NGDC_INDICES_command
    !-------------------------------------------------------------------------
    subroutine MHD_INDICES_command
        use ModUtilities, ONLY: CON_stop

        character (len=100) :: cTempLine
        character (len=100), dimension(100) :: cTempLines
        integer:: iError

        cTempLines(1) = '#MHD_INDICES'
        call read_var('UpstreamFile',cTempLine)
        cTempLines(2) = cTempLine
        cTempLines(3) = " "
        cTempLines(4) = "#END"

        call IO_set_inputs(cTempLines)
        call read_MHDIMF_Indices(iError)
        if(iError /= 0)then
            write(*,*)'IM: read_MHDIMF_Indices iError=', iError
            call CON_stop('IM: Could not read '//trim(cTempLine))
        end if
    end subroutine MHD_INDICES_command
    !-------------------------------------------------------------------------
    subroutine KYOTO_DST_command
        use ModImIndices, ONLY: UseDstKyoto, NameDstFile, &
                                     read_dst_wdc_kyoto
        use ModUtilities, ONLY: CON_stop

        integer:: iError

        call read_var('UseDstKyoto',UseDstKyoto)
        if(UseDstKyoto)then
           call read_var('NameDstFile',NameDstFile)

           call read_dst_wdc_kyoto(iError)
                
           if (iError /= 0) then
              write(*,*)'IM: read_dst_wdc_kyotom iError=', iError
              call CON_stop('IM: Could not read '//trim(NameDstFile))
           end if
        end if
    end subroutine KYOTO_DST_command
    !-------------------------------------------------------------------------
    subroutine KYOTO_AE_command
        use ModImIndices,      ONLY: UseAeKyoto,NameAeFile,read_ae_wdc_kyoto
        use ModUtilities, ONLY: CON_stop

        integer:: iError

        call read_var('UseAeKyoto',UseAeKyoto)
        if(UseAeKyoto)then
           call read_var('NameAeFile',NameAeFile)
           call read_ae_wdc_kyoto(iError)
                
           if (iError /= 0) then
              write(*,*)'IM: read_ae_wdc_kyoto iError=', iError
              call CON_stop('IM: Could not read '//trim(NameAeFile))
           end if
        end if
    end subroutine KYOTO_AE_command
    !-------------------------------------------------------------------------
    subroutine POTSDAM_KP_AP_F107_command
        use ModImIndices,      ONLY: UseKpApF107IndicesFile, &
                                     read_kpapf107_indices_file

        !use Kp, Ap, and F107 from Geomagnetic Observatory Niemegk, &
        !GFZ German Research Centre for Geosciences
        !note this file has all data from 1937, but does need to be
        ! periodically
        !updated from https://www.gfz-potsdam.de/en/section/geomagnetism/data-products-services/geomagnetic-kp-index
        call read_var('UseKpApF107IndicesFile',UseKpApF107IndicesFile)
        call read_kpapf107_indices_file
    end subroutine POTSDAM_KP_AP_F107_command
    !-------------------------------------------------------------------------
    subroutine SOLARWIND_command
        real :: DensitySW, VelSW, BxSW, BySW, BzSW

        call read_var('DensitySW',DensitySW)
        call read_var('VelSW',VelSW)
        call read_var('BxSW',BxSW)
        call read_var('BySW',BySW)
        call read_var('BzSW',BzSW)
        
        call IO_set_imf_bx_single(BxSW)
        call IO_set_imf_by_single(BySW)
        call IO_set_imf_bz_single(BzSW)
        call IO_set_sw_v_single(abs(VelSW))
        call IO_set_sw_n_single(abs(DensitySW))
    end subroutine SOLARWIND_command
    !-------------------------------------------------------------------------
    subroutine SINGLEGRID_command
        use ModCimiGrid, ONLY: np,DoReadParamGrid, &
                               LatitudeDeg, dLatitudeDeg, &
                               LongitudeRad, dLongitudeRad

        call read_var('DoReadParamGrid',DoReadParamGrid)
        if (np.gt.1 .and. DoReadParamGrid) then
            write(*,*) 'WARNING: #SINGLEGRID command used, '
            write(*,*) '   CIMI is configured with more than one grid point! '
            write(*,*) '   ( setting doReadParamGrid=F )'
            DoReadParamGrid=.false.
        else
            if (DoReadParamGrid) call read_var('LatitudeDeg',LatitudeDeg)
            if (DoReadParamGrid) call read_var('dLatitudeDeg',dLatitudeDeg)
            if (DoReadParamGrid) call read_var('LongitudeRad',LongitudeRad)
            if (DoReadParamGrid) call read_var('dLongitudeRad',dLongitudeRad)
        endif
    end subroutine SINGLEGRID_command
    !-------------------------------------------------------------------------
    subroutine SMOOTH_command
        use ModCimiTrace,	 ONLY: UseSmooth, SmoothWindow

        !do you want to smooth the solar wind for Tsy model
        call read_var('UseSmooth',UseSmooth)
        if (UseSmooth) call read_var('SmoothWindow',SmoothWindow)
    end subroutine SMOOTH_command
    !-------------------------------------------------------------------------
    subroutine BMODEL_command
        use ModUtilities, ONLY: lower_case, CON_stop
        use ModCimiTrace, ONLY: iMod, UseCorotation, UsePotential
        use ModGmCimi,    ONLY: UseGm

        character (len=3) :: NameModel

        call read_var('NameModel',NameModel)!t96,t04,MHD,Dip
        !call read_var('UseFixedB',UseFixedB)!T=fixed B config or 
                                            !F=changing B config
        call lower_case( NameModel )
        if (NameModel == 'dip') then
           iMod=0
           ! Checks the next two lines to see if cor and/or pot
           ! variables should be set to 0.
           call read_var('UseCorotation',UseCorotation)
           call read_var('UsePotential',UsePotential)
        elseif (NameModel == 't96') then
           iMod=1
        elseif(NameModel == 't04') then
           iMod=2
        elseif(NameModel == 'mhd')then
           iMod=3
           UseGm=.true.
        else
           call CON_stop('Error: Model not found') 
        endif
    end subroutine BMODEL_command
    !-------------------------------------------------------------------------
    subroutine IEMODEL_command
        use ModIeCimi, ONLY: UseWeimer

        call read_var('UseWeimer',UseWeimer)
    end subroutine IEMODEL_command
    !-------------------------------------------------------------------------
    subroutine GIMME_command
        use GIMME_cimi_interface, ONLY: DtGimmePlot, UseGimme

        call read_var('UseGimme',UseGimme)
        call read_var('DtGimmePlot',DtGimmePlot)
    end subroutine GIMME_command
    !-------------------------------------------------------------------------
    subroutine SAVELOG_command
        use ModCimiPlot, ONLY: DtLogOut, DoSaveLog

        call read_var('DtLogOut',DtLogOut)
        DoSaveLog = .true.
    end subroutine SAVELOG_command
    !-------------------------------------------------------------------------
    subroutine SAVEPLOT_command
        use ModUtilities, ONLY: lower_case
        use ModCimiPlot,     ONLY: DtOutput, DoSavePlot,&
                  DoSaveFlux, DtFluxOutput, DoSaveSeparateFluxFiles,&
                  DoSavePSD, DtPSDOutput, DoSaveSeparatePSDFiles,&
                  DoSaveVLDrift, DtVLDriftOutput, DoSaveSeparateVLDriftFiles,&
                  DoSaveVPDrift, DtVPDriftOutput, DoSaveSeparateVPDriftFiles,&
                  DoSavePreci, DtPreciOutput, DoSaveSeparatePreciFiles,&
                  DoSaveEq, DoSaveIono,&
                  DoSaveLstar, DtLstarOutput, DoSaveSeparateLstarFiles
        use ModCimiPlanet,     ONLY: nspec
        use ModUtilities,    ONLY: CON_stop
        use ModCimi,         ONLY: DoCalcPrecip
        use ModPlasmasphere, ONLY: DoSavePlas, DtPlasOutput

        integer :: iCIMIPlotType, nCIMIPlotType
        character (len=100) :: StringCIMIPlot
        real :: DtOutputCIMIPlot
        logical :: DoSaveSeparateFiles

        call read_var( 'nCIMIPlotType', nCIMIPlotType )

        CIMI_PLOTTYPE: do iCIMIPlotType = 1, nCIMIPlotType
           
           call read_var( 'StringPlot', StringCIMIPlot )
           call lower_case( StringCIMIPlot )
           call read_var( 'DtOutput', DtOutputCIMIPlot )
           
           if ( index( StringCIMIPlot, 'fls'  ) > 0 .or. &
                index( StringCIMIPlot, 'flux' ) > 0 ) then
              
              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )
              
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
                     index( StringCIMIPlot, 'e' ) > 0 ) then
                 
                 DoSaveFlux( nspec ) = .true.
                 DtFluxOutput( nspec ) = DtOutputCIMIPlot
                 DoSaveSeparateFluxFiles( nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'h' ) > 0 ) then
                 
                 DoSaveFlux( 1 ) = .true.
                 DtFluxOutput( 1 ) = DtOutputCIMIPlot
                 DoSaveSeparateFluxFiles( 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'o' ) > 0 ) then
                 
                 if ( nspec .ge. 3 ) then
                    
                    DoSaveFlux( 2 ) = .true. 
                    DtFluxOutput( 2 ) = DtOutputCIMIPlot
                    DoSaveSeparateFluxFiles( 2 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_stop( 'O+ not configured; Recompile CIMI '// &
                         'with ModEarthHO or ModEarthHOHe options.' )
                    
                 endif
                 
              elseif &
                   ( index( StringCIMIPlot, 'he' ) > 0 ) then
                 
                 if ( nspec .gt. 3 ) then
                    
                    DoSaveFlux( 3 ) = .true. 
                    DtFluxOutput( 3 ) = DtOutputCIMIPlot
                    DoSaveSeparateFluxFiles( 3 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_stop( 'He+ not configured; Recompile CIMI '//&
                         'with ModEarthHOHe option.' )
                    
                 endif
                 
              else
                 
                 call CON_stop( 'No flux species information; STOPPING' )
                 
              endif PLOT_FLUX_SPECIES
              
           elseif &
                ( index( StringCIMIPlot, 'psd' ) > 0 ) then
              
              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )
              
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
                     index( StringCIMIPlot, 'e' ) > 0 ) then
                 DoSavePSD( nspec ) = .true.
                 DtPSDOutput( nspec ) = DtOutputCIMIPlot
                 DoSaveSeparatePSDFiles( nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'h' ) > 0 ) then
                 
                 DoSavePSD( 1 ) = .true.
                 DtPSDOutput( 1 ) = DtOutputCIMIPlot
                 DoSaveSeparatePSDFiles( 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'o' ) > 0 ) then

                 if ( nspec .ge. 3 ) then
                    DoSavePSD( 2 ) = .true. 
                    DtPSDOutput( 2 ) = DtOutputCIMIPlot
                    DoSaveSeparatePSDFiles( 2 ) = &
                         DoSaveSeparateFiles
                 else
                    
                    call CON_stop( 'O+ not configured; Recompile CIMI '// &
                         'with ModEarthHO or ModEarthHOHe options.' )
                    
                 endif
                 
              elseif &
                   ( index( StringCIMIPlot, 'he' ) > 0 ) then
                 
                 if ( nspec .gt. 3 ) then
                    
                    DoSavePSD( 3 ) = .true. 
                    DtPSDOutput( 3 ) = DtOutputCIMIPlot
                    DoSaveSeparatePSDFiles( 3 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_stop( 'He+ not configured; Recompile CIMI '//&
                         'with ModEarthHOHe option.' )
                    
                 endif
                 
              else
                 
                 call CON_stop( 'No PSD species information; STOPPING' )
                 
              endif PLOT_PSD_SPECIES
              
           elseif &
                ( index( StringCIMIPlot, 'vl'  ) > 0 .or. &
                  index( StringCIMIPlot, 'vldrift' ) > 0 ) then
              
              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )
              
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
                     index( StringCIMIPlot, 'e' ) > 0 ) then
                 
                 DoSaveVLDrift( nspec ) = .true.
                 DtVLDriftOutput( nspec ) = DtOutputCIMIPlot
                 DoSaveSeparateVLDriftFiles( nspec ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'h' ) > 0 ) then
                 
                 DoSaveVLDrift( 1 ) = .true.
                 DtVLDriftOutput( 1 ) = DtOutputCIMIPlot
                 DoSaveSeparateVLDriftFiles( 1 ) = &
                      DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'o' ) > 0 ) then
                 
                 if ( nspec .ge. 3 ) then
                    
                    DoSaveVLDrift( 2 ) = .true. 
                    DtVLDriftOutput( 2 ) = DtOutputCIMIPlot
                    DoSaveSeparateVLDriftFiles( 2 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_stop( 'O+ not configured; Recompile CIMI '// &
                         'with ModEarthHO or ModEarthHOHe options.' )
                    
                 endif
              elseif &
                   ( index( StringCIMIPlot, 'he' ) > 0 ) then
                 
                 if ( nspec .gt. 3 ) then
                    
                    DoSaveVLDrift( 3 ) = .true. 
                    DtVLDriftOutput( 3 ) = DtOutputCIMIPlot
                    DoSaveSeparateVLDriftFiles( 3 ) = &
                         DoSaveSeparateFiles
                    
                 else
                    
                    call CON_stop( 'He+ not configured; Recompile CIMI '//&
                         'with ModEarthHOHe option.' )
                    
                 endif
                 
              else
                 
                 call CON_stop( 'No VLDrift species information; STOPPING' )
                 
              endif PLOT_VLDRIFT_SPECIES
              
           elseif &
                ( index( StringCIMIPlot, 'vp'  ) > 0 .or. &
                  index( StringCIMIPlot, 'vpdrift' ) > 0 ) then

              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )

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
                     index( StringCIMIPlot, 'e' ) > 0 ) then

                 DoSaveVPDrift( nspec ) = .true.
                 DtVPDriftOutput( nspec ) = DtOutputCIMIPlot
                 DoSaveSeparateVPDriftFiles( nspec ) = &
                      DoSaveSeparateFiles

              elseif &
                   ( index( StringCIMIPlot, 'h' ) > 0 ) then

                 DoSaveVPDrift( 1 ) = .true.
                 DtVPDriftOutput( 1 ) = DtOutputCIMIPlot
                 DoSaveSeparateVPDriftFiles( 1 ) = &
                      DoSaveSeparateFiles

              elseif &
                   ( index( StringCIMIPlot, 'o' ) > 0 ) then

                 if ( nspec .ge. 3 ) then

                    DoSaveVPDrift( 2 ) = .true.
                    DtVPDriftOutput( 2 ) = DtOutputCIMIPlot
                    DoSaveSeparateVPDriftFiles( 2 ) = &
                         DoSaveSeparateFiles

                 else

                    call CON_stop( 'O+ not configured; Recompile CIMI '// &
                         'with ModEarthHO or ModEarthHOHe options.' )

                 endif

              elseif &
                   ( index( StringCIMIPlot, 'he' ) > 0 ) then

                 if ( nspec .gt. 3 ) then

                    DoSaveVPDrift( 3 ) = .true.
                    DtVPDriftOutput( 3 ) = DtOutputCIMIPlot
                    DoSaveSeparateVPDriftFiles( 3 ) = &
                         DoSaveSeparateFiles

                 else

                    call CON_stop( 'He+ not configured; Recompile CIMI '//&
                         'with ModEarthHOHe option.' )

                 endif

              else

                 call CON_stop( 'No VPDrift species information; STOPPING' )

              endif PLOT_VPDRIFT_SPECIES

           elseif &
                ( index( StringCIMIPlot, 'precipitation'  ) > 0 .or. &
                  index( StringCIMIPlot, 'precip' ) > 0 .or. &
                  index( StringCIMIPlot, 'preci' ) > 0 ) then

              DoCalcPrecip = .true.
              call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )

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
                     index( StringCIMIPlot, 'e' ) > 0 ) then

                 DoSavePreci( nspec ) = .true.
                 DtPreciOutput( nspec ) = DtOutputCIMIPlot
                 DoSaveSeparatePreciFiles( nspec ) = &
                      DoSaveSeparateFiles

              elseif &
                   ( index( StringCIMIPlot, 'h' ) > 0 ) then

                 DoSavePreci( 1 ) = .true.
                 DtPreciOutput( 1 ) = DtOutputCIMIPlot
                 DoSaveSeparatePreciFiles( 1 ) = &
                      DoSaveSeparateFiles

              elseif &
                   ( index( StringCIMIPlot, 'o' ) > 0 ) then

                 if ( nspec .ge. 3 ) then

                    DoSavePreci( 2 ) = .true.
                    DtPreciOutput( 2 ) = DtOutputCIMIPlot
                    DoSaveSeparatePreciFiles( 2 ) = &
                         DoSaveSeparateFiles

                 else

                    call CON_stop( 'O+ not configured; Recompile CIMI '// &
                         'with ModEarthHO or ModEarthHOHe options.' )

                 endif

              elseif &
                   ( index( StringCIMIPlot, 'he' ) > 0 ) then

                 if ( nspec .gt. 3 ) then

                    DoSavePreci( 3 ) = .true.
                    DtPreciOutput( 3 ) = DtOutputCIMIPlot
                    DoSaveSeparatePreciFiles( 3 ) = &
                         DoSaveSeparateFiles

                 else

                    call CON_stop( 'He+ not configured; Recompile CIMI '//&
                         'with ModEarthHOHe option.' )

                 endif

              else

                 call CON_stop( 'No Preci species information; STOPPING' )

              endif PLOT_PRECI_SPECIES

           elseif &
                ( index( StringCIMIPlot, '2d'  ) > 0 ) then

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
                 
              elseif &
                   ( index( StringCIMIPlot, 'lstar'  ) > 0 .or. &
                     index( StringCIMIPlot, 'l*'  ) > 0 ) then
                 
                 call read_var( 'DoSaveSeparateFiles', DoSaveSeparateFiles )
                 
                 DoSaveLstar = .true.
                 DtLstarOutput = DtOutputCIMIPlot
                 DoSaveSeparateLstarFiles = DoSaveSeparateFiles
                 
              elseif &
                   ( index( StringCIMIPlot, 'coreplas' ) > 0 ) then
                 
                 DoSavePlas = .true.
                 DtPlasOutput = DtOutputCIMIPlot
                 
              else

                 call CON_stop( 'No CIMI 2D Plot information; STOPPING' )
                 
              endif PLOT_2D

           else
              
              call CON_stop( 'No plot type specified.' )
              
           endif

        end do CIMI_PLOTTYPE

        DoSavePlot = .true.

!!  END OF SAVEPLOT ROUTINE        

    end subroutine SAVEPLOT_command
    !-------------------------------------------------------------------------
    subroutine DTSATOUT_command
        use ModImSat, ONLY: DtSatOut

        call read_var('DtSatOut', DtSatOut)
    end subroutine DTSATOUT_command
    !-------------------------------------------------------------------------
    subroutine VERBOSELSTAR_command
        use ModLstar, ONLY: DoVerboseLstar

        call read_var('DoVerboseLstar',DoVerboseLstar)
    end subroutine VERBOSELSTAR_command
    !-------------------------------------------------------------------------
    subroutine INITIALF2_command
        use ModCimiInitialize, ONLY: IsEmptyInitial, IsGmInitial,&
                                     IsDataInitial, IsRBSPData
        use ModCimiPlanet,     ONLY: nspec
        use ModUtilities,      ONLY: CON_stop

        call read_var('IsEmptyInitial',IsEmptyInitial)
        call read_var('IsGmInitial',   IsGmInitial)
        call read_var('IsDataInitial', IsDataInitial)
        if (IsDataInitial) &
             call read_var('IsRBSPData', IsRBSPData)
        
        !IsDataInitial only works with EarthHO or EarthH configurations
        if (nspec > 3 .and. IsDataInitial) &
             call CON_stop('IsDataInitial only works with EarthHO or EarthH')  
    end subroutine INITIALF2_command
    !-------------------------------------------------------------------------
    subroutine INITIALLSTAR_command
        use ModCimiInitialize, ONLY: DoLstarInitialization

        call read_var('DoLstarInitialization',DoLstarInitialization)
    end subroutine INITIALLSTAR_command
    !-------------------------------------------------------------------------
    subroutine TYPEBOUNDARY_command
        use ModCimiTrace,	 ONLY: UseEllipse

        character (len=7) :: TypeBoundary

        call read_var('TypeBoundary',TypeBoundary)
        if(TypeBoundary == 'Ellipse') then
           UseEllipse = .true.
        else
           UseEllipse = .false.
        endif
    end subroutine TYPEBOUNDARY_command
    !-------------------------------------------------------------------------
    subroutine PLASMASHEET_command
        use ModCimiBoundary, ONLY: UseYoungEtAl, UseBoundaryEbihara
        use ModCimi,	     ONLY: IsStandAlone

        call read_var('UseYoungEtAl',UseYoungEtAl)
        if(IsStandAlone)then
           ! When Standalone specify boundary type
           call read_var('UseBoundaryEbihara',UseBoundaryEbihara)
        endif
    end subroutine PLASMASHEET_command
    !-------------------------------------------------------------------------
    subroutine BOUNDARYKAPPA_command
        use ModUtilities, ONLY: CON_stop
        use ModCimiPlanet,	     ONLY: nspec
        use ModCimiBoundary,	 ONLY: &
                                UseDynamicKappa,UseKappaTotalFit, &
                                KappaMin_I, KappaMax_I, KappaConst_I, &
                                DoRescaleKappaTemp
        integer :: ispec

        call read_var('UseDynamicKappa',UseDynamicKappa)
        if(UseDynamicKappa)then
           call read_var('UseKappaTotalFit',UseKappaTotalFit)
           do iSpec = 1, nspec
              call read_var('KappaMin',KappaMin_I(iSpec))
              call read_var('KappaMax',KappaMax_I(iSpec))
              if(KappaMin_I(iSpec) <= 1.5) &
                  call CON_stop('IM: KappaMin must be >= 1.5')
           end do
        else
           do iSpec = 1, nspec
              call read_var('KappaConst',KappaConst_I(iSpec))
              if(KappaConst_I(iSpec) <= 1.5 .and. KappaConst_I(iSpec) > 0.0) &
                  call CON_stop('IM: KappaConst must be >= 1.5')
           end do
        end if
        call read_var('DoRescaleKappaTemp',DoRescaleKappaTemp)
    end subroutine BOUNDARYKAPPA_command
    !-------------------------------------------------------------------------
    subroutine RESTART_command
        use ModCimiRestart,	ONLY: IsRestart
        use ModImSat,		ONLY: ReadRestartSat

        call read_var('IsRestart',IsRestart) !T:Continuous run, F: Initial Run
        if(IsRestart) call read_var('ReadRestartSat',ReadRestartSat)
    end subroutine RESTART_command
    !-------------------------------------------------------------------------
    subroutine SAVERESTART_command
        use ModCimi,	    ONLY: IsStandAlone
        use ModCimiRestart,	ONLY: DtSaveRestart
        use ModUtilities,   ONLY: CON_stop

        ! when in standalone mode read restart save frequency
        if(IsStandAlone) then
           call read_var('DtSaveRestart',DtSaveRestart)
        else
           call CON_stop(&
                       'ERROR: Restart save frequency only set in standalone')
        end if
    end subroutine SAVERESTART_command
    !-------------------------------------------------------------------------
    subroutine LIMITER_command
        use ModCimi,		 ONLY: UseMcLimiter, BetaLimiter

        call read_var('UseMcLimiter', UseMcLimiter)
        if(UseMcLimiter) call read_var('BetaLimiter', BetaLimiter)
    end subroutine LIMITER_command
    !-------------------------------------------------------------------------
    subroutine MINIMUMPRESSURETOGM_command
        use ModCimi, ONLY: Pmin

        call read_var('MinimumPressureToGM', Pmin)
    end subroutine MINIMUMPRESSURETOGM_command
    !-------------------------------------------------------------------------
    subroutine ELECTRONTEMPERATURERATIO_command
        use ModCimiPlanet, ONLY: nspec, tFactor_I

        real :: ElectronTempRatio
        integer :: iSpec

        call read_var('ElectronTempRatio', ElectronTempRatio)
        do iSpec = 1, nspec - 1
            tFactor_I(iSpec) = 1/(1 + ElectronTempRatio)
        end do
        tFactor_I(nspec) = ElectronTempRatio/(1 + ElectronTempRatio)
    end subroutine ELECTRONTEMPERATURERATIO_command
    !-------------------------------------------------------------------------
    subroutine GMTEMPERATUREFACTORS_command
        use ModCimiPlanet,   ONLY: nspec, tFactor_I

        integer :: iSpec

         do iSpec = 1, nspec
            call read_var('tFactor', tFactor_I(iSpec))
         end do
    end subroutine GMTEMPERATUREFACTORS_command
    !-------------------------------------------------------------------------
    subroutine TIMESIMULATION_command
        use ModCimi,		 ONLY: time

        call read_var('TimeSimulation',time)
    end subroutine TIMESIMULATION_command
    !-------------------------------------------------------------------------
    subroutine PRERUNSAT_command
        use ModImSat, ONLY: DtSatOut, DoWritePrerunSat,&
                            UsePrerunSat, DtReadSat, DoWriteSats

        call read_var('DoWritePrerunSat',DoWritePrerunSat)
        if(.not.DoWritePrerunSat) call read_var('UsePrerunSat',UsePrerunSat)
        if(UsePrerunSat) then
           call read_var('DtReadSat',   DtReadSat)
           DtSatout = DtReadSat
           DoWriteSats = .true.
        end if
    end subroutine PRERUNSAT_command
    !-------------------------------------------------------------------------
    subroutine PRERUNFIELD_command
        use ModPrerunField,	ONLY: DoWritePrerun, UsePrerun, DtRead
        use ModGmCIMI,	    ONLY: UseGm

        call read_var('DoWritePrerun',DoWritePrerun)
        if(.not.DoWritePrerun) call read_var('UsePrerun',UsePrerun)
        if(UsePrerun)          call read_var('DtRead',   DtRead)
        if(UsePrerun .or. DoWritePrerun) UseGm=.true.
    end subroutine PRERUNFIELD_command
    !-------------------------------------------------------------------------
    subroutine DECAY_command
        use ModCimi, ONLY: UseDecay, DecayTimescale

        call read_var('UseDecay',UseDecay)
        if ( UseDecay ) &
             call read_var('DecayTimescale in seconds', DecayTimescale)
    end subroutine DECAY_command
    !-------------------------------------------------------------------------
    subroutine ELECTRONDECAY_command
        use ModCimi, ONLY: UseElectronDecay

        call read_var('UseElectronDecay', UseElectronDecay)
    end subroutine ELECTRONDECAY_command
    !-------------------------------------------------------------------------
    subroutine FLC_command
        use ModCimi, ONLY: UseFLC

        call read_var('UseFLC',UseFLC)
    end subroutine FLC_command
    !-------------------------------------------------------------------------
    subroutine COMPOSITION_command
        use ModCimiPlanet, ONLY: nspec, dFactor_I
        use ModUtilities,  ONLY: CON_stop

        integer :: iSpec
        character (len=5) :: TypeComposition='FIXED'

        call read_var('TypeComposition', TypeComposition, IsUpperCase=.true.)
        select case(TypeComposition)
        case('FIXED')
           do iSpec = 1, nspec
              call read_var('DensityFraction', dFactor_I(iSpec))
           end do
        case default
           call CON_stop(&
                       NameSub//': unknown TypeComposition='//TypeComposition)
        end select
    end subroutine COMPOSITION_command
    !-------------------------------------------------------------------------
    subroutine STRONGDIFFUSION_command
        use ModCimi, ONLY: UseStrongDiff

        call read_var('UseStrongDiff',UseStrongDiff)
    end subroutine STRONGDIFFUSION_command
    !-------------------------------------------------------------------------
    subroutine WAVES_command
        use CIMI_waves,	ONLY: UseWaves, UseHiss, UseChorus,&
                              NameHissFile, NameChorusFile, UseKpIndex

        call read_var('UseWaves',UseWaves)
        if(UseWaves) then
           call read_var('UseHiss',  UseHiss)
           call read_var('UseChorus',  UseChorus)
           call read_var('NameHissFile', NameHissFile)
           call read_var('NameChorusFile', NameChorusFile)
           call read_var('UseKpIndex',UseKpIndex)
        end if
    end subroutine WAVES_command
    !-------------------------------------------------------------------------
    subroutine DIAGONALIZEDDIFFUSION_command
        use CIMI_waves, ONLY: UseDiagDiffusion

        call read_var('UseDiagDiffusion',UseDiagDiffusion)
    end subroutine DIAGONALIZEDDIFFUSION_command
    !-------------------------------------------------------------------------
    subroutine DIAGDIFFUSIONTEST_command
        use ModDiagDiff, ONLY: UsePitchAngleDiffusionTest,&
                               UseEnergyDiffusionTest

        call read_var('UsePitchAngleDiffusionTest',&
                       UsePitchAngleDiffusionTest)
        call read_var('UseEnergyDiffusionTest',&
                       UseEnergyDiffusionTest)
    end subroutine DIAGDIFFUSIONTEST_command
    !-------------------------------------------------------------------------
    subroutine ENERGYGRID_command
        use ModCimiGrid, ONLY: MinIonEnergy, MaxIonEnergy

        call read_var('MinIonEnergy (in keV), MinElectronEnergy x10',&
             MinIonEnergy)
        call read_var('MaxIonEnergy (in keV), MaxElectronEnergy x10',&
             MaxIonEnergy)
    end subroutine ENERGYGRID_command
    !-------------------------------------------------------------------------
    subroutine SETENERGYGRID_command
        use ModCimiGrid,  ONLY: neng, UseLogEGrid, MinIonEnergy, MaxIonEnergy
        use ModUtilities, ONLY: CON_stop

        call read_var('neng', neng)
        if (neng .lt. 2) &
             call CON_stop('ERROR: neng must be greater than 1.')
        call read_var('UseLogEGrid', UseLogEGrid)
        call read_var('MinIonEnergy (in keV), MinElectronEnergy x10',&
             MinIonEnergy)
        call read_var('MaxIonEnergy (in keV), MaxElectronEnergy x10',&
             MaxIonEnergy)
    end subroutine SETENERGYGRID_command
    !-------------------------------------------------------------------------
    subroutine RBSPENERGYGRID_command
        use ModCimiGrid, ONLY: neng, UseRBSPGrid, nRBSPEnergy

        call read_var('UseRBSPGrid', UseRBSPGrid)
        if (UseRBSPGrid) neng = nRBSPEnergy
    end subroutine RBSPENERGYGRID_command
    !-------------------------------------------------------------------------
    subroutine ELECTRONGRIDRATIO_command
        use ModCimiGrid, ONLY: ElectronGridRatio

        call read_var('ElectronGridRatio', ElectronGridRatio)
    end subroutine ELECTRONGRIDRATIO_command
    !-------------------------------------------------------------------------
    subroutine IMTIMESTEP_command
        use ModCimi, ONLY: dt,dtmax

        call read_var('IMDeltaT [s]',dt)
        call read_var('IMDeltaTMax [s]',dtmax)
    end subroutine IMTIMESTEP_command
    !-------------------------------------------------------------------------
    subroutine TESTFIELDLINE_command
        use ModCimiTrace, ONLY: iLatTest, iLonTest

        call read_var('iLatTest',iLatTest)
        call read_var('iLonTest',iLonTest)
    end subroutine TESTFIELDLINE_command
    !-------------------------------------------------------------------------
    subroutine PLASMAPAUSE_command
        use ModUtilities, ONLY: CON_stop
        use ModPlasmasphere, ONLY: PlasmaPauseDensity, TypePlasmapause

        call read_var('TypePlasmapause',TypePlasmapause, IsLowerCase=.true.)
        select case(TypePlasmapause)
         case('relative')
         case('absolute')
            call read_var('PlasmaPauseDensity [m^-3]',PlasmaPauseDensity)
         case default
            call CON_stop(NameSub//': unknown TypePlasmapause='// &
                TypePlasmapause)
         end select
    end subroutine PLASMAPAUSE_command
    !-------------------------------------------------------------------------
    subroutine PLASMAPAUSEDENSITY_command
        use ModPlasmasphere, ONLY: PlasmaPauseDensity, TypePlasmapause

        TypePlasmapause = 'absolute'
        call read_var('DensityP [m^3]',PlasmaPauseDensity)
    end subroutine PLASMAPAUSEDENSITY_command
    !-------------------------------------------------------------------------
    subroutine COREPLASMASPHERE_command
        use ModPlasmasphere, ONLY: UseCorePsModel, PlasMinDensity

        call read_var('UseCorePsModel',UseCorePsModel)
        call read_var('PlasMinDensity',PlasMinDensity)
    end subroutine COREPLASMASPHERE_command
    !-------------------------------------------------------------------------
    subroutine SETRB_command
        use ModCimiGrid, ONLY: rb

        call read_var('rb [R_E]', rb)
    end subroutine SETRB_command
    !-------------------------------------------------------------------------
    subroutine SETBOUNDARYPARAMS_command
        use ModCimiTrace, ONLY: DeltaRMax, xmltlim

        call read_var('DeltaRMax', DeltaRMax)
        call read_var('DeltaMLTmax', xmltlim)
    end subroutine SETBOUNDARYPARAMS_command
    !-------------------------------------------------------------------------
    subroutine PRECIPITATION_command
        use ModCimiTrace, ONLY: UseAltitudePrecip, UsePrecipEnergyLoss

       call read_var('UseAltitudePrecip', UseAltitudePrecip)
       if(UseAltitudePrecip) &
          call read_var('UsePrecipEnergyLoss', UsePrecipEnergyLoss)
    end subroutine PRECIPITATION_command
    !-------------------------------------------------------------------------
    subroutine STRICTDRIFT_command
        use ModCimi, ONLY: IsStrictDrift

        call read_var('IsStrictDrift',IsStrictDrift) ! .T : STOP when f2 < 0
    end subroutine STRICTDRIFT_command
    !-------------------------------------------------------------------------
    subroutine LATITUDINALGRID_command
        use ModCimiGrid, ONLY: DoDefineVarNpower, varNpower, xlatmax

        call read_var('DoDefineVarNpower',DoDefineVarNpower) 
        call read_var('varNpower',varNpower)   ! n in L = 1/cos(xlat)**n
        call read_var('xlatmax',xlatmax)       ! upper boundary latitude
    end subroutine LATITUDINALGRID_command
    !-------------------------------------------------------------------------
    subroutine VERBOSELATGRID_command
        use ModCimiGrid, ONLY: DoVerboseLatGrid

        call read_var( 'DoVerboseLatGrid', DoVerboseLatGrid )
    end subroutine VERBOSELATGRID_command
    !-------------------------------------------------------------------------
    subroutine HIGHERORDERDRIFT_command
        use ModInterFlux, ONLY: UseHigherOrder, iOrderLat, iOrderLon
        use ModUtilities, ONLY: CON_stop

        call read_var('UseHigherOrder',UseHigherOrder)
        if ( UseHigherOrder ) then
           
           call read_var( 'iOrderLat', iOrderLat )   ! order in latitude
           call read_var( 'iOrderLon', iOrderLon )   ! order in longitude
           if ( iOrderLat .eq. iOrderLon ) then
              if ( ( iOrderLat .ne. 2 ) .and. ( iOrderLat .ne. 7 ) ) &
                   call CON_stop('IM: ERROR: iOrderLat = iOrderLon = 2 '//&
                   	'or iOrderLat = iOrderLon = 7')
           else
              call CON_stop('IM: ERROR: iOrderLat = iOrderLon = 2 '//&
                   'or iOrderLat = iOrderLon = 7')
           endif
           
        endif
    end subroutine HIGHERORDERDRIFT_command
    !-------------------------------------------------------------------------
    subroutine DRIFTSCHEME_command
        use ModInterFlux, ONLY: UseHigherOrder, iOrderLat, iOrderLon
        use ModUtilities, ONLY: CON_stop

        call read_var( 'iOrderLat', iOrderLat )   ! order in latitude
        call read_var( 'iOrderLon', iOrderLon )   ! order in longitude
        if ( iOrderLat .eq. iOrderLon ) then
           if ( iOrderLat .eq. 7 ) UseHigherOrder=.true.
           if ( ( iOrderLat .ne. 2 ) .and. ( iOrderLat .ne. 7 ) ) &
                call CON_stop('IM: ERROR: iOrderLat = iOrderLon = 2 '//&
                'or iOrderLat = iOrderLon = 7')
        else
           call CON_stop('IM: ERROR: iOrderLat = iOrderLon = 2 '//&
                'or iOrderLat = iOrderLon = 7')
        endif
        write(*,*) "UseHigherOrder: ", UseHigherOrder
    end subroutine DRIFTSCHEME_command
    !-------------------------------------------------------------------------
end subroutine CIMI_set_parameters
end module CIMI_SetParams
