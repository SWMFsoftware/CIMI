Module ModCoupleSami
    implicit none 
    private !except

    logical, public :: DoCoupleSami = .false.
    
    integer, public :: iProc0CIMI,iProc0SAMI, iCommGlobal

    !public methods
    public :: cimi_set_global_mpi
    public :: cimi_get_init_for_sami
    public :: cimi_send_to_sami
  contains
    !==========================================================================
    subroutine cimi_set_global_mpi(iProc0CimiIN,iProc0SamiIN, iCommGlobalIN)
      integer, intent(in) ::iProc0CimiIN,iProc0SamiIN, iCommGlobalIN
      !------------------------------------------------------------------------
      iProc0CIMI=iProc0CimiIN
      iProc0Sami=iProc0SamiIN
      iCommGlobal=iCommGlobalIN

      ! Set the coupling to be true
      DoCoupleSami=.true.
    end subroutine cimi_set_global_mpi
    !==========================================================================
    ! only call by proc 0
    subroutine cimi_get_init_for_sami
      use ModCrcmGrid, ONLY: xlatr, nLat=>np, phi, nLon=>nt
      use ModImTime,   ONLY: iStartTime_I
      use ModMPI
      integer :: iError
      integer :: iStatus_I(MPI_STATUS_SIZE)
      !------------------------------------------------------------------------
      
      ! Send the starttime from cimi to sami
      call MPI_send(iStartTime_I,7,MPI_INTEGER,iProc0SAMI,&
           1,iCommGlobal,iError)

      ! Send the grid parameters from cimi to sami
      call MPI_send(nLat,1,MPI_INTEGER,iProc0SAMI,&
           2,iCommGlobal,iError)
      
      call MPI_send(nLon,1,MPI_INTEGER,iProc0SAMI,&
           3,iCommGlobal,iError)

      call MPI_send(xlatr,nLat,MPI_REAL,iProc0SAMI,&
           4,iCommGlobal,iError)

      call MPI_send(phi,nLon,MPI_REAL,iProc0SAMI,&
           5,iCommGlobal,iError)
      
    end subroutine cimi_get_init_for_sami
    !==========================================================================
    ! 
    subroutine cimi_send_to_sami
      use ModCrcmGrid, ONLY: nLat=>np, phi, nLon=>nt, iProc
      use ModIeCrcm,   ONLY: pot
      use ModMPI
      integer :: iError
      integer :: iStatus_I(MPI_STATUS_SIZE)
      !------------------------------------------------------------------------
      
      !if(iProc==0) write(*,*) 'cimi_send_to_sami',nLat,nLon
      ! Send the starttime from cimi to sami
      if(iProc ==0) call MPI_send(pot,nLat*nLon,MPI_REAL,iProc0SAMI,&
           1,iCommGlobal,iError)

    end subroutine cimi_send_to_sami
    
  end Module ModCoupleSami
  
