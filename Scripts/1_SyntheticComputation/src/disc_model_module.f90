MODULE disc_model_module
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: moho_topography,&
            cmb_topography

CONTAINS

  SUBROUTINE moho_topography(idisc)

    !  read fin_moho_variation w.r.t. 1-D moho (in km, positive value represents
    !  thicker crust).
    !  The file contains "modified" coefficients of real spherical harmonics. 
    !
    !  Moho variation is always stored in ddis(1,:,:) and idisc(1)
    !
    !  USE NS_D_MOHO, NS,ND, RA, fin_moho_variation

    USE modellayer_module, ONLY : NMOHO
    IMPLICIT NONE
    INCLUDE "coupling.h" 
    
    INTEGER, INTENT(INOUT) :: idisc(:) !ND
 
    INTEGER :: ierr,k,l,m,id
    REAL    :: rtmp
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: Ac,Bc
    
    COMPLEX dB2(0:NK,0:NS,0:NS),dS2(0:NK,0:NS,0:NS),drho(0:NK,0:NS,0:NS)
    COMPLEX ddisc(ND,0:NS,0:NS)
    COMMON/mantlemodel/dB2,dS2,drho,ddisc
  
    id=1 

    ALLOCATE(Ac(id,0:NS_D_MOHO,0:NS_D_MOHO),Bc(id,0:NS_D_MOHO,0:NS_D_MOHO), &
         STAT=ierr)
    
    IF (ierr/=0) STOP 'allocate error 2 in SUB moho_topography STOP'
    IF ( ND > 0 ) THEN
       Ac=0.0
       Bc=0.0
       ddisc(id,:,:)=CMPLX(0.0,0.0)    
       OPEN(10,file=fin_moho_variation,STATUS='OLD',ACTION='READ',IOSTAT=ierr)
       IF(ierr/=0) STOP 'open error in reading MOHO topography, STOP'
       DO k=id,id
          DO l=0,NS_D_MOHO
             READ(10,*,IOSTAT=ierr) Ac(k,l,0),(Ac(k,l,m),Bc(k,l,m),m=1,l)
          END DO
          READ(10,*,IOSTAT=ierr) rtmp
          IF (ierr>=0) STOP 'reading crust.dat error, STOP'
          
          DO l=0,NS_D_MOHO
             m=0
             ddisc(k,l,m)=-cmplx(Ac(k,l,m),0.0)*1000.0/RA
             DO m=1,l
                ddisc(k,l,m)=-0.5*cmplx(Ac(k,l,m),-Bc(k,l,m))*1000.0/RA
             END DO
          END DO
       END DO
       CLOSE(10)
       idisc(1)=NMOHO
       WRITE(*,*) 'Moho variation has been considered...'
    END IF
    
    DEALLOCATE(Ac,Bc)

  END SUBROUTINE moho_topography


  SUBROUTINE cmb_topography(idisc, cmb_topo_mod)

    USE modellayer_module, ONLY : NCMB
    IMPLICIT NONE
    INCLUDE "coupling.h"

    CHARACTER(len=128), INTENT(IN) :: cmb_topo_mod
    INTEGER, INTENT(INOUT) :: idisc(:)

    INTEGER :: ierr,k,l,m,id
    REAL    :: rtmp
    REAL, DIMENSION(:,:), ALLOCATABLE :: Rs,Is
    character*300 :: fcmbr,fcmbi

    COMPLEX dB2(0:NK,0:NS,0:NS),dS2(0:NK,0:NS,0:NS),drho(0:NK,0:NS,0:NS)
    COMPLEX ddisc(ND,0:NS,0:NS)
    COMMON/mantlemodel/dB2,dS2,drho,ddisc
      id=2   
    ALLOCATE(Rs(0:NS,0:NS),Is(0:NS,0:NS), &
             STAT=ierr)

    fcmbr = trim(cmb_topo_mod)//'dcmb_re.dat'
    fcmbi = trim(cmb_topo_mod)//'dcmb_im.dat'

    print *, 'CMB topography set to:  ', trim(cmb_topo_mod)

    IF (ierr/=0) STOP 'allocate error 2 in SUB cmb_topography STOP'
    IF ( ND > 0 ) THEN
       ddisc(id,:,:)=CMPLX(0.0,0.0)
       write(*,*) 'considering CMB:',fcmbr

       OPEN(10,file=fcmbr,STATUS='OLD',ACTION='READ',IOSTAT=ierr)
       IF(ierr/=0) STOP 'open error in reading CMB_r topography, STOP'
       OPEN(11,file=fcmbi,STATUS='OLD',ACTION='READ',IOSTAT=ierr)
       IF(ierr/=0) STOP 'open error in reading CMB_i topography, STOP'
       DO k=id,id
          DO l=0,6
             read(10,*,iostat=ierr) (Rs(l,m),m=0,l)
             !IF (ierr>=0) STOP 'reading cmb_r.dat error, STOP'
             read(11,*,iostat=ierr) (Is(l,m),m=0,l)
             !IF (ierr>=0) STOP 'reading cmb_i.dat error, STOP'
          END DO
          DO l=0,6
             m=0
             ddisc(k,l,m)=cmplx(Rs(l,m),Is(l,m))*1000.0/RA
             DO m=1,l
                ddisc(k,l,m)=cmplx(Rs(l,m),Is(l,m))*1000.0/RA
             END DO
          END DO
       END DO
       CLOSE(10)
       CLOSE(11)

       idisc(2)=NCMB

       WRITE(*,*) 'CMB variation has been considered...'

    ENDIF

    DEALLOCATE(Rs,Is)

  END SUBROUTINE cmb_topography


END MODULE disc_model_module
