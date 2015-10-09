Module ModCoupleCimi
    implicit none 
    private !except

    logical, public :: DoCoupleCimi = .false.
    
    integer, public ::  iStartTime_I(7)  =(/1976,6,28,0,0,0,0/)
    
    integer, public :: iProc0CIMI,iProc0SAMI, iCommGlobal
    
    ! cimi Latitude (radians) and phi (radians) grid
    real, allocatable :: LatCimi_C(:),PhiCimi_C(:)
    integer :: nLatCIMI, nLonCIMI
    
    !public methods
    public :: sami_set_global_mpi
    public :: sami_put_init_from_cimi
  contains
    !==========================================================================
    ! call by all procs
    subroutine sami_set_global_mpi(iProc0CimiIN,iProc0SamiIN, iCommGlobalIN)
      integer, intent(in) ::iProc0CimiIN,iProc0SamiIN, iCommGlobalIN
      !------------------------------------------------------------------------
      iProc0CIMI=iProc0CimiIN
      iProc0Sami=iProc0SamiIN
      iCommGlobal=iCommGlobalIN

      ! Set the coupling to be true
      DoCoupleCimi=.true.
    end subroutine sami_set_global_mpi
    !========================================================================
    subroutine sami_put_init_from_cimi
      use ModSAMI, ONLY: iComm,iProc
      use ModMpi
      integer ::  iStartTimeCIMI_I(7),iError
      integer :: iStatus_I(MPI_STATUS_SIZE)
      !----------------------------------------------------------------------

      ! recieve and then bcast the starttime from cimi
      if (iProc == 0) call MPI_recv(iStartTime_I,7,MPI_INTEGER,iProc0CIMI,&
           1,iCommGlobal,iStatus_I,iError)
      call MPI_bcast(iStartTime_I,7,MPI_LOGICAL,0,iComm,iError)
      
      ! recieve grid size
      if (iProc == 0) call MPI_recv(nLatCIMI,1,MPI_INTEGER,iProc0CIMI,&
           2,iCommGlobal,iStatus_I,iError)

      if (iProc == 0) call MPI_recv(nLonCIMI,1,MPI_INTEGER,iProc0CIMI,&
           3,iCommGlobal,iStatus_I,iError)
      
      ! bcast grid size
      call MPI_bcast(nLatCIMI,1,MPI_LOGICAL,0,iComm,iError)
      call MPI_bcast(nLonCIMI,1,MPI_LOGICAL,0,iComm,iError)
      
      ! allocate CIMI grid
      allocate(LatCimi_C(nLatCIMI))
      allocate(PhiCimi_C(nLonCIMI))

      if (iProc == 0) call MPI_recv(LatCimi_C,nLatCIMI,MPI_REAL,iProc0CIMI,&
           4,iCommGlobal,iStatus_I,iError)

      if (iProc == 0) call MPI_recv(PhiCimi_C,nLonCIMI,MPI_REAL,iProc0CIMI,&
           5,iCommGlobal,iStatus_I,iError)

    end subroutine sami_put_init_from_cimi
    
  end Module ModCoupleCimi
  
