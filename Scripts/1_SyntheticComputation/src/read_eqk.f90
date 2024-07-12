MODULE read_eqk

  ! ================================================================ !
  ! Utility scripts for reading source-receiver information          !
  ! ================================================================ !

  IMPLICIT NONE

  TYPE multisource_type
    real*8 :: elat,elon,depth ! event lat, event lon, depth
    real*8 :: moment_tensor(6)
  END TYPE multisource_type

  TYPE polarization_type
    real*8 :: pola1,pola2,pola3
  END TYPE polarization_type

  TYPE hdr_type
    character(len=8) :: ntwk, name
    real*8 :: stla, stlo, stel
  END TYPE hdr_type

  TYPE(multisource_type), allocatable :: multisource(:)
  TYPE(polarization_type), allocatable :: polarization(:)
  TYPE(hdr_type), allocatable :: hdr(:)

  PRIVATE
  PUBLIC :: read_stations,&
            read_cmt,&
            multisource,&
            hdr,&
            polarization,&
            receiver_polarization,&
            polarization_type,&
            multisource_type

CONTAINS
  
  SUBROUTINE read_stations(stn_file,num_station)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! read station info from station file               !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! input/output:
    character(len=*), intent(in) :: stn_file
    integer, intent(out) :: num_station
    ! internal:
    integer :: fid=21
    integer :: ierr
    integer :: istn=0

    ! begin:
    open(unit=fid,file=stn_file,status='old',action='read')
    read(unit=fid,fmt=*) ! 1st line is not used
    do
      read(unit=fid,fmt=*,iostat=ierr)
      if(ierr<0) exit
      if(ierr>0) stop 'reading station problem'
      istn = istn + 1
    end do

    num_station = istn
    write(*,*) 'Num station = ', num_station

    allocate(hdr(num_station))
    allocate(polarization(num_station))

    rewind(unit=fid)
    read(unit=fid,fmt=*)
    do istn=1,num_station
      read(unit=fid,fmt=*) hdr(istn)%ntwk,&
                           hdr(istn)%name,&
                           hdr(istn)%stla,&
                           hdr(istn)%stlo,&
                           hdr(istn)%stel,&
                           polarization(istn)%pola1,&
                           polarization(istn)%pola2,&
                           polarization(istn)%pola3 
    end do
    close(fid)

  END SUBROUTINE read_stations

  ! ================================================================ !

  SUBROUTINE read_cmt(cmt_file)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! read source info from cmt file                    !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! input/output:
    character(len=*), intent(in) :: cmt_file
    ! internal: 
    integer :: fid = 12,i
    real*8 :: input(9) ! 3 event params + 6 Moment tensor

    allocate(multisource(1))
 
    open(unit=fid,file=cmt_file,status='old',action='read')
    call get_cmt_info(fid,'latitude',multisource(1)%elat)
    call get_cmt_info(fid,'longitude',multisource(1)%elon)
    call get_cmt_info(fid,'depth',multisource(1)%depth)
    call get_cmt_info(fid,'Mrr',multisource(1)%moment_tensor(1))
    call get_cmt_info(fid,'Mtt',multisource(1)%moment_tensor(2))
    call get_cmt_info(fid,'Mpp',multisource(1)%moment_tensor(3))
    call get_cmt_info(fid,'Mrt',multisource(1)%moment_tensor(4))
    call get_cmt_info(fid,'Mrp',multisource(1)%moment_tensor(5))
    call get_cmt_info(fid,'Mtp',multisource(1)%moment_tensor(6))
    close(unit=fid)

  CONTAINS

    SUBROUTINE get_cmt_info(fid,cheader,cmtinfo)
      ! input/output:
      integer, intent(in) :: fid
      character(len=*), intent(in) :: cheader
      real*8, intent(out) :: cmtinfo
      ! internal:
      character(len=128) :: string
      integer :: itmp1, itmp2, ierr

      read(fid,FMT='(A)') string
      itmp1=index(string,':')-1
      itmp2=len_trim(string)
      if (string(1:itmp1) == cheader) then
        read(string(itmp1+2:itmp2),*,iostat=ierr) cmtinfo
        if (ierr/=0) stop 'reading cmt file error1, STOP!'
      else
        write(*,*) ' reading cmt file error2: ', trim(cheader)  
        stop
      end if

    END SUBROUTINE get_cmt_info
    
  END SUBROUTINE read_cmt

  ! ================================================================ !


  FUNCTION receiver_polarization(a,b,c)
    real*8, intent(in) :: a,b,c
    real*8 :: receiver_polarization(3)
    receiver_polarization(1)= a
    receiver_polarization(2)= b
    receiver_polarization(3)= c

  END FUNCTION receiver_polarization

  ! ================================================================ !

END MODULE read_eqk
