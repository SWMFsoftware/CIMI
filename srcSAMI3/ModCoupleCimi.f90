Module ModCoupleCimi
    implicit none 
    private !except

    logical, public :: DoCoupleCimi = .false.
    
    integer, public ::  iStartTime_I(7)  =(/1976,6,28,0,0,0,0/)
    
    integer, public :: iProc0CIMI,iProc0SAMI, iCommGlobal
    
  contains
    !========================================================================
    subroutine sami_put_init_from_cimi
      use ModSAMI, ONLY: iComm,iProc
      integer ::  iStartTimeCIMI_I(7),iError
      !----------------------------------------------------------------------
      
      ! Set the coupling to be true
      DoCoupleCimi=.true.
      
      ! recieve and then bcast the starttime from cimi
      if (iProc = 0) call MPI_recv(iStartTimeCIMI_I,7,MPI_INTEGER,iProc0CIMI,&
           1,iCommGlobal,iStatus_I,iError)
      
      ! Set the start time from cimi start time
      iStartTime_I = iStartTimeCIMI_I 
      
      call MPI_bcast(iStartTime_I,7,MPI_LOGICAL,0,iComm,iError)
      
      
      
      
      
      
    end subroutine sami_put_init_from_cimi
    
  end Module ModCoupleCimi
  
