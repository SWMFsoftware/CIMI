Module ModCimiRestart
  implicit none
  
  SAVE

  logical :: IsRestart = .false.
  real    :: DtSaveRestart=-1.0

  character(len=100), public:: NameRestartInDir="IM/restartIN/"
  character(len=100), public:: NameRestartOutDir="IM/restartOUT/"
contains
  !============================================================================

  subroutine cimi_read_restart
    use ModCimiPlanet,ONLY: nspec
    use ModCimiGrid,  ONLY: np,nt,nm,nk,neng,d4Element_C
    use ModCimi,      ONLY: f2, phot, Pressure_IC, PressurePar_IC, FAC_C, &
         Ppar_IC, Bmin_C, &
         eTimeAccumult_ICI, eChangeOperator_VICI, eChangeGlobal, &
         pTimeAccumult_ICI, pChangeOperator_VICI, &
         driftin, driftout, rbsumGlobal, rcsumGlobal, nOperator, &
         PreF, PreP, Eje1, nTotalCell
    use ModCimiTrace, ONLY: iba, iw2
    use ModGmCimi,    ONLY: Den_IC
    use ModIoUnit,    ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file
    use ModCimiGrid,  ONLY: iProc,nProc,iComm
    use ModMpi

    integer :: iError
    character(len=4) :: iProcStr

    character(len=*), parameter:: NameSub = 'cimi_read_restart'    
    !--------------------------------------------------------------------------
    !When nProc>1, proc0 reads and then bcasts restart infor
    !when only 1 proc is used then just read restart info
    
     write(iProcStr,'(I4.4)') iProc

     call open_file(file=trim(NameRestartInDir)//'data.restart'//trim(iProcStr),&
          status='old',form='unformatted', NameCaller=NameSub)

     read(UnitTmp_) f2  
     read(UnitTmp_) Den_IC
     read(UnitTmp_) phot
     read(UnitTmp_) Pressure_IC           
     read(UnitTmp_) PressurePar_IC           
     read(UnitTmp_) FAC_C             
     read(UnitTmp_) iba
     read(UnitTmp_) Ppar_IC
     read(UnitTmp_) Bmin_C
     read(UnitTmp_) eTimeAccumult_ICI
     read(UnitTmp_) eChangeOperator_VICI
     read(UnitTmp_) eChangeGlobal
     read(UnitTmp_) pTimeAccumult_ICI
     read(UnitTmp_) pChangeOperator_VICI
     read(UnitTmp_) rbsumGlobal
     read(UnitTmp_) rcsumGlobal
     read(UnitTmp_) driftin
     read(UnitTmp_) driftout
     read(UnitTmp_) iw2
     read(UnitTmp_) PreF
     read(UnitTmp_) PreP
     read(UnitTmp_) Eje1
     call close_file
    
  end subroutine cimi_read_restart
  !============================================================================
  subroutine cimi_write_restart

    use ModCimiPlanet,ONLY: nspec
    use ModCimiGrid,  ONLY: np,nt,nm,nk,neng,MinLonPar,MaxLonPar
    use ModCimi,      ONLY: f2,time, phot, Pressure_IC, PressurePar_IC, &
         FAC_C, Ppar_IC, Bmin_C, &
         eTimeAccumult_ICI, eChangeOperator_VICI, eChangeGlobal, &
         pTimeAccumult_ICI, pChangeOperator_VICI, &
         driftin, driftout, rbsumGlobal, rcsumGlobal, nOperator, &
         PreF, PreP, Eje1
    use ModCimiTrace,ONLY: iba,iw2    
    use ModGmCimi,    ONLY: Den_IC
    use ModIoUnit,    ONLY: UnitTmp_
    use ModUtilities, ONLY: open_file, close_file, cTab
    use ModCimiGrid,  ONLY: iProc,nProc,iComm,nLonPar,nLonPar_P,nLonBefore_P
    use ModImSat,     ONLY: DoWriteSats,nImSats,NameSat_I,SatLoc_3I
    use ModMpi

    integer ::iSendCount, iSpecies, iM, iK, iError, iOperator, iEnergy
    real :: BufferSend_C(np,nt),BufferRecv_C(np,nt)
    integer :: BufferSend_I(nt),BufferRecv_I(nt)
    integer,allocatable :: iReceiveCount_P(:),iDisplacement_P(:)

    ! Initialize necessary variables for saving restart.sat file
    integer :: iSat,iRow
    character(len=4) :: iProcStr

    character(len=*), parameter:: NameSub = 'cimi_write_restart'
    !--------------------------------------------------------------------------

     write(iProcStr,'(I4.4)') iProc
     call open_file(file=trim(NameRestartOutDir)//'data.restart'//trim(iProcStr), &
          form='unformatted', NameCaller=NameSub)
     write(UnitTmp_) f2
     write(UnitTmp_) Den_IC
     write(UnitTmp_) phot
     write(UnitTmp_) Pressure_IC
     write(UnitTmp_) PressurePar_IC
     write(UnitTmp_) FAC_C
     write(UnitTmp_) iba                
     write(UnitTmp_) Ppar_IC
     write(UnitTmp_) Bmin_C
     write(UnitTmp_) eTimeAccumult_ICI
     write(UnitTmp_) eChangeOperator_VICI
     write(UnitTmp_) eChangeGlobal
     write(UnitTmp_) pTimeAccumult_ICI
     write(UnitTmp_) pChangeOperator_VICI
     write(UnitTmp_) rbsumGlobal
     write(UnitTmp_) rcsumGlobal
     write(UnitTmp_) driftin
     write(UnitTmp_) driftout
     write(UnitTmp_) iw2
     write(UnitTmp_) PreF
     write(UnitTmp_) PreP
     write(UnitTmp_) Eje1
     call close_file

     if(iProc == 0) then
       call open_file(file=trim(NameRestartOutDir)//'restart.H', &
            NameCaller=NameSub)
       write(UnitTmp_,'(a)') '#TIMESIMULATION'
       write(UnitTmp_,'(es15.8,a)') time, cTab//cTab//'tSimulation'
       write(UnitTmp_,*)
       write(UnitTmp_,'(a)') '#RESTART'
       write(UnitTmp_,'(a)') 'T'//cTab//cTab//cTab//'IsRestart'
       write(UnitTmp_,'(l1,a)') DoWriteSats, cTab//cTab//cTab//'DoReadRestartSatellite'
       call close_file

       if (DoWriteSats) then
          call open_file(file=trim(NameRestartOutDir)//'restart.sat', &
               form="unformatted", NameCaller=NameSub) 
          write(UnitTmp_) nImSats
       
          do iSat=1,nImSats
             
             write(UnitTmp_) NameSat_I(iSat)
             
          end do
       
          do iSat=1,nImSats
          
             do iRow=1,2
             
                write(UnitTmp_) SatLoc_3I(1:4,iRow,iSat)
                
             end do
             
          end do
          
          call close_file

       end if
       
    endif
  end subroutine cimi_write_restart

end Module ModCimiRestart
