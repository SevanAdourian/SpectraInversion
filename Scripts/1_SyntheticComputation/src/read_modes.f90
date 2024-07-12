MODULE read_cluster

! cleaned up module just for cluster 
! information
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_modesetup


CONTAINS

  SUBROUTINE read_modesetup(setup_file,tcluster)

    ! HL: NOV 2018
    ! simplified version for tides (just reads clusters)
    ! we consider just one cluster (full coupling)
    !
    ! reads the file 'setup' in ./bin/ where line 1
    ! is number of clusters which should be 1
    ! assumes full coupling so only one cluster

    include "coupling.h"  !ONLY: CCLUSTMAX

    ! INPUTS
    character(len=*), INTENT(IN) :: setup_file
    ! OUTPUTS
    character(len=CCLUSTMAX), INTENT(OUT) :: tcluster
    ! OTHER
    integer :: ncpl,fid,iread,ict
    character(len=10) :: modestring
    character(len=CCLUSTMAX) :: string    

    fid=22
    ict=1
    open(unit=fid,file=setup_file,status='old',action='read')
    read(fid,*) ncpl
    do iread=1,ncpl
       read(fid,*) modestring
       string(ict:ict+LEN_TRIM(modestring))=trim(modestring)
       ict=ict+LEN_TRIM(modestring)
    enddo
    close(fid)
    
    tcluster=trim(string(1:ict))
    write(*,*) 'Cluster modes:'
    write(*,*) tcluster(1:ict)

  END SUBROUTINE read_modesetup


END MODULE read_cluster 


