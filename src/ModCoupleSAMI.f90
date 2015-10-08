Module ModCoupleSami
    implicit none 
    private !except

    logical, public :: DoCoupleSami = .false.
    
    integer, public :: iProc0CIMI,iProc0SAMI, iCommGlobal
    
  contains
    !========================================================================
    subroutine cimi_get_init_for_sami(iStartTime_I)
      use ModCrcmGrid, ONLY: iComm,iProc
      integer,intent(in) ::  iStartTime_I(7)
      integer :: iError
      !----------------------------------------------------------------------
      
      ! Set the coupling to be true
      DoCoupleSami=.true.
      
      ! recieve and then bcast the starttime from cimi
      if (iProc = 0) call MPI_send(iStartTime_I,7,MPI_INTEGER,iProc0SAMI,&
           1,iCommGlobal,iStatus_I,iError)
      
      
      
      
    end subroutine cimi_get_init_for_sami
    
  end Module ModCoupleSami
  
