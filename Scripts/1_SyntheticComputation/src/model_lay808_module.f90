MODULE model_lay808_module
  ! MAY 2022
  ! Makes model for every layer in the background PREM 808 model

  PRIVATE
  PUBLIC :: lay808_mantle_model, &
             lay808_mantle_basis
CONTAINS

  SUBROUTINE lay808_mantle_model(is_corr_crust,is_corr_cmb,idisc,&
                                 lay808_model,cmb_topo_mod)

  ! ================================================================ !
  ! --> common block: /mantlemodel/dB2,dS2,drho,ddisc
  !
  !       B2=kappa/rho is the squared bulk velocity
  !       S2=mu/rho is the squared shear velocity
  !       dB2 is the relative squared bulk velocity perturbation
  !       delta(B2)/B2 and dS2 is the relative squared shear
  !       velocity perturbation delta(S2)/S2
  ! ================================================================ !

  USE disc_model_module, ONLY : moho_topography,cmb_topography
  USE modellayer_module, ONLY : NCMB
 
  IMPLICIT NONE

  include "coupling.h"

  CHARACTER(len=128), intent(in) :: lay808_model, cmb_topo_mod

  LOGICAL, INTENT(IN) :: is_corr_crust
  LOGICAL, INTENT(IN) :: is_corr_cmb
  INTEGER, INTENT(OUT):: idisc(ND)

  REAL    :: rtmp
  INTEGER :: ierr
  INTEGER ::k,l,m
  COMPLEX*16 dB2(0:NK,0:NS,0:NS),dS2(0:NK,0:NS,0:NS),drho(0:NK,0:NS,0:NS)
  COMPLEX*16 ddisc(ND,0:NS,0:NS)
  REAL*8 :: Rs(0:NK,0:NS,0:NS),Is(0:NK,0:NS,0:NS)
  REAL*8 :: R0
  INTEGER :: kk,j
  
  character(len=100) :: vsfile,vbfile,rofile
  integer :: fout1,fout2,fout3,lmax_mod
  real*8 :: dXr1,dXi1,dXr2,dXi2,dXr3,dXi3
  real :: fjunk

  ! BM23:
  real*8 :: dpar

  common/mantlemodel/dB2,dS2,drho,ddisc

  write(*,*) 'Model used: ', lay808_model

  lmax_mod = 6

  dS2(:,:,:)=cmplx(0.0,0.0)
  dB2(:,:,:)=cmplx(0.0,0.0)
  drho(:,:,:)=cmplx(0.0,0.0)

  fout1=325
  fout2=326
  fout3=327
  ! the files obtained from HL were renamed so that
  ! the code doesn't change
  ! now reading model degree 12 
  vsfile = trim(lay808_model)//'dlnvs_csh.dat' 
  vbfile = trim(lay808_model)//'dlnvb_csh.dat'
  rofile = trim(lay808_model)//'dlnro_csh.dat'
 
  open(fout1,file=trim(vsfile),status='unknown')
  open(fout2,file=trim(vbfile),status='unknown')
  open(fout3,file=trim(rofile),status='unknown')

  do k=0,NK
    ! do l=0, m=0
    l=0
    m=0
    read(fout1,*) fjunk,fjunk,fjunk,dXr1,dXi1
    read(fout2,*) fjunk,fjunk,fjunk,dXr2,dXi2
    read(fout3,*) fjunk,fjunk,fjunk,dXr3,dXi3
    dS2(k,l,m)=cmplx(0.0,0.0)
    dB2(k,l,m)=cmplx(0.0,0.0)
    drho(k,l,m)=cmplx(0.0,0.0)
    do l=1,lmax_mod!NS
      m=0
      read(fout1,*) fjunk,fjunk,fjunk,dXr1,dXi1
      read(fout2,*) fjunk,fjunk,fjunk,dXr2,dXi2
      read(fout3,*) fjunk,fjunk,fjunk,dXr3,dXi3
      dS2(k,l,m) = 2.d0*cmplx(dXr1,0.0)
      dB2(k,l,m) = 2.d0*cmplx(dXr2,0.0)
      drho(k,l,m)= cmplx(dXr3,0.0)
      do m=1,l
        read(fout1,*) fjunk,fjunk,fjunk,dXr1,dXi1
        read(fout2,*) fjunk,fjunk,fjunk,dXr2,dXi2
        read(fout3,*) fjunk,fjunk,fjunk,dXr3,dXi3
        dS2(k,l,m)  = 2.0*cmplx(dXr1,dXi1)
        dB2(k,l,m)  = 2.0*cmplx(dXr2,dXi2)
        drho(k,l,m) = cmplx(dXr3,dXi3)
      enddo
    enddo
  enddo
  close(fout1)
  close(fout2)
  close(fout3)

  ! =================================
  ! You can also put other discontinuity layer here
  ! !->idisc,ddisc
  !==================================

  ! crust correction
  IF (is_corr_crust) THEN
    CALL moho_topography(idisc) !-> idisc(s1), ddisc
  END IF

  IF (is_corr_cmb) THEN
    CALL cmb_topography(idisc, cmb_topo_mod)
  END IF

  return

END SUBROUTINE lay808_mantle_model

SUBROUTINE lay808_mantle_basis(r,t)

  USE modellayer_module, ONLY : NR,NCMB,NMOHO
  IMPLICIT NONE
  include "coupling.h"
  REAL*8, POINTER :: r(:)
  REAL*8, POINTER :: t(:,:)

  INTEGER :: nlay,i,j

  IF (.NOT.ASSOCIATED(t) .OR. .NOT.ASSOCIATED(r)) &
    STOP 't and r have to be associated, STOP in lay25_mantle_basis'

  t(:,:)=0.0
  do i=0,NK
    j = NCMB + i
    t(j,i)=1.
  enddo
  
  return

END SUBROUTINE lay808_mantle_basis

END MODULE model_lay808_module
