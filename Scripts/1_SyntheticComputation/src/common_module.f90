MODULE modellayer_module

  ! Vs_moho: velocity between V_moho- and V_moho+ (in m/s)
  ! NR: number of radial knots in eigenfunctions and its correpsonding model
  ! NOC: layer index of ocean floor in the 1D model (r-)
  !      if ocean floor does not exist, NOC=NR 
  ! NMOHO: layer index of moho in the 1D model (r-)
  ! NCMB: layer index of cmb in the 1D model  ***** (r+) ******
  ! NICB: layer index of icb in the 1D model (r-)
  ! kdis: index of discontinuity layer in the 1D model (r-)
  !
  ! These variables are defined by read_prem_fun, which reads
  ! model stored in binayr file of "funSph".   

  IMPLICIT NONE
  REAL*8, PARAMETER :: Vs_moho=4000.0  
  INTEGER :: NR,NOC,NMOHO,NCMB,NICB
  INTEGER, DIMENSION(:), POINTER :: kdis
  
  PRIVATE
  PUBLIC ::  NR,NOC,NMOHO,NCMB,NICB,Vs_moho,kdis

END MODULE modellayer_module

!----------------------------------
!----------------------------------
!----------------------------------

MODULE prem_module
  ! rho: density
  ! epsilon: ellipticity
  ! eta: used by the process of calculating "epsilon"
  ! g: gravity
  ! c01: C
  ! c02: 2N
  ! c03: A-N
  ! c04: -F
  ! c05: -L
  ! Qkappa: quality factor of kappa
  ! Qmu:    quality facto of mu
  !
  ! "non-dimensional" 1D model defined by read_prem_fun
  
  IMPLICIT NONE
  ! DIMENSION(NR)
  REAL*8, DIMENSION(:), POINTER :: rho,epsilon,eta,g, &
                                 c01,c02,c03,c04,c05,Qkappa,Qmu
END MODULE prem_module

!----------------------------------
!----------------------------------
!----------------------------------

MODULE eigfun_module
  ! defined in SUBROUTINE decipher_cluster 
  ! used in initialize_rot_ell

  ! eigfun_radius: radius of the eigenfunctions (non-dimensional) 
  ! eigfun_(u,du,v,dv,w,dw,p,dp) : eigenfunctions 
  !      In fact, v=V/k; dv=dV/k
  !               w=W/k; dw=dW/k
  ! eiginfo: information of mode
  !        w: eigenfrequency (in rad/s)
  !        q: quality factor
  !        n: radial order of mode
  !        l: angular degree of mode
  !    itype: =1 for Spheroidal and Radial modes
  !           =0 for Toridial modes
  ! eigac: eigenfuncions which consider gravit effect on seismometer
  !           at the top of the solid earth (_tsod)
  !        contains
  !             u=U*_tsod, 
  !             v=V*_tsod/k
  !             w=W_tsod/k 
  !             p=P
  !        see Chap 10.4 in Dahlen & Tromp, 1998 for details   
  !

  IMPLICIT NONE
  TYPE modeinfo_type
     INTEGER :: n, l, itype 
     REAL    :: w, q
  END TYPE modeinfo_type
  
  TYPE modeac_type
     REAL :: u,v,w,p
  END TYPE modeac_type

  REAL*8, DIMENSION(:), ALLOCATABLE, TARGET   :: eigfun_radius         
  REAL*8, DIMENSION(:,:), ALLOCATABLE, TARGET :: &
                                       eigfun_u, eigfun_du,&
                                       eigfun_v, eigfun_dv,&
                                       eigfun_w, eigfun_dw,&  
                                       eigfun_p, eigfun_dp
  TYPE(modeinfo_type), DIMENSION(:), ALLOCATABLE, TARGET :: eiginfo
  TYPE(modeac_type), DIMENSION(:), ALLOCATABLE, TARGET :: eigac
  PRIVATE
  PUBLIC :: eigfun_radius, &         
            eigfun_u, eigfun_du,&
            eigfun_v, eigfun_dv,&
            eigfun_w, eigfun_dw,&  
            eigfun_p, eigfun_dp, &
            modeinfo_type, &
            eiginfo, &
            eigac, &
            deallocate_eigfun, &
            allocate_eigfun
  
  CONTAINS

    SUBROUTINE allocate_eigfun(nr,nm)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nr, nm
      INTEGER :: ierr
      CALL deallocate_eigfun
      ALLOCATE(eigfun_radius(nr), &
               eigfun_u(nr,nm), &
               eigfun_du(nr,nm), &
               eigfun_v(nr,nm), &
               eigfun_dv(nr,nm), &
               eigfun_w(nr,nm), &
               eigfun_dw(nr,nm), &
               eigfun_p(nr,nm), &
               eigfun_dp(nr,nm), &
               eiginfo(nm),&
               eigac(nm), &
               STAT=ierr)
      IF (ierr/=0) STOP 'allocate error of eigfun. STOP'

    END SUBROUTINE allocate_eigfun

    SUBROUTINE deallocate_eigfun

      IF (ALLOCATED(eigfun_radius)) DEALLOCATE(eigfun_radius)
      IF (ALLOCATED(eigfun_u)) DEALLOCATE(eigfun_u)
      IF (ALLOCATED(eigfun_v)) DEALLOCATE(eigfun_v)
      IF (ALLOCATED(eigfun_w)) DEALLOCATE(eigfun_w)
      IF (ALLOCATED(eigfun_p)) DEALLOCATE(eigfun_p)
      IF (ALLOCATED(eigfun_du)) DEALLOCATE(eigfun_du)
      IF (ALLOCATED(eigfun_dv)) DEALLOCATE(eigfun_dv)
      IF (ALLOCATED(eigfun_dw)) DEALLOCATE(eigfun_dw)
      IF (ALLOCATED(eigfun_dp)) DEALLOCATE(eigfun_dp)          
      IF (ALLOCATED(eiginfo)) DEALLOCATE(eiginfo) 
      IF (ALLOCATED(eigac)) DEALLOCATE(eigac)
 
    END  SUBROUTINE deallocate_eigfun

END MODULE eigfun_module

