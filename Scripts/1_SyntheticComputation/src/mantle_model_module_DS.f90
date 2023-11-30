MODULE mantle_model_module_DS

  ! ================================================================ !

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_prem_fun, &
            find_mdis, &
            initialize_mantle_model, &
            allocate_3D_structure, &
            deallocate_3D_structure, &
            set_3D_structure

  ! ================================================================ !
  ! variables associated with the 3D model                           !
  ! ================================================================ !
  integer, save, public :: sst
  integer, save, public :: nbnd
  integer, dimension(:), allocatable, save, public :: bnd
  real*8, dimension(:), allocatable, save, public :: rho_el,d_el
  real*8, dimension(:), allocatable, save, public :: mu_el,kap_el
  complex*16, dimension(:), allocatable, save, public :: c1,c2,c3
  complex*16, dimension(:), allocatable, save, public :: c4,c5
  complex*16, dimension(:,:,:), allocatable, save, public :: rho_st
  complex*16, dimension(:,:,:), allocatable, save, public :: mu_st
  complex*16, dimension(:,:,:), allocatable, save, public :: dmu_st
  complex*16, dimension(:,:,:), allocatable, save, public :: kap_st
  complex*16, dimension(:,:,:), allocatable, save, public :: dkap_st
  complex*16, dimension(:,:,:), allocatable, save, public :: d_st
  complex*16, dimension(:), allocatable, save, public :: dkap, dmu

CONTAINS

  SUBROUTINE read_prem_fun(r)

    ! HL MAR 2018
    ! edited mainly how dispersion is treated
    ! we do not treat dispersion here.

    ! we read the reference model here from MINEOS
    ! outputs.
    ! the reference model is read in at wref (in coupling.h)
    ! the eigenfrequencies are read in at wf or omega 
    ! (i.e., the forcing frequency of the system)

    ! INTENT(OUT) :: r
    ! r (dimensionaless) is radius of model with dimension 
    ! of n
    ! 
    ! write the follwing variables into prem_module 
    !   rho,epsilon,eta,g,c01,c02,c03,c04,c05,Qkappa,Qmu
    !   These density and elastic parameters are non-dimensional.
    !     c01 = C  =  C{0000}
    !     c02 = 2N =  C{++--}
    !     c03 = A-N = C{+-+-}
    !     c04 = -F  = C{+-00}
    !     c05 = -L  = C{+0-0}
    !   The index in C{} represents the components in genearlized system.
    !   Because prem is model of transversely isotropic, the order is 
    !   zero, i.e. C{abcd}, where a+b+c+d=0
    !   Please refer to Mochizuki, Geophys. J. R. astr. Soc. (1986).
    !   
    !   epsilon: ellipticity as a function of radius (eq. 14.20)
    !   g: graivity
    !   for earth model, epsilon(ra)~1/300
    !     g(ra)~1.33 (non-dimensional)
    !
    !   write NR,NMOMO-,NOC-,NCMB+,NICB-, kdis into modellayer_module
    
    USE modellayer_module, ONLY : NR,NMOHO,NCMB,NICB,kdis
    USE prem_module
    USE futil_module, ONLY : intgrl,deriv
    IMPLICIT NONE     
    include "coupling.h" 

    REAL*8, POINTER :: r(:)
    REAL*8 :: z,bom,exponent1,i_rho,i_radau
    REAL*8, DIMENSION(:),ALLOCATABLE :: k,radau,s1,s2,s3,&
         i_rho_r,i_radau_r
 
    INTEGER :: n
    REAL*8 :: eps = 1.e-7
    REAL*8,ALLOCATABLE :: model_store(:,:)  !model_store(NRMAX,9)
    CHARACTER(LEN=128) :: bin_file

    INTEGER :: ierr, i, j, ndis

    !!! MINEOS NOTES: 
    !!! THESE ARE FROM MINEOS'S UNFORMATTED READ - MUST STAY
    !!! IN THIS FORMAT:
    !!! VALUES TO BE TRANSFERED WILL HAVE "_ms"
    character(len=8) :: kfmt
    integer :: itmp,nicfun,nocfun,ifanis
    real :: rtmp,trefl
    real, dimension(:,:), allocatable :: mstore_ms
    !!!
     

    ALLOCATE(mstore_ms(NRMAX,9),STAT=ierr)
    IF(ierr/=0) STOP 'allocate error 1 in SUB read_prem_fun, STOP!'
    ALLOCATE(model_store(NRMAX,9),STAT=ierr)
    IF(ierr/=0) STOP 'allocate error 1 in SUB read_prem_fun, STOP!'

    !!! in this part we do not read the entire eigenfunction file
    !!! we just extract PREM

    ! read binary model file
    bin_file=TRIM(funSph)
      
    OPEN(UNIT=1,FILE=bin_file,STATUS='OLD',FORM='UNFORMATTED', &
         ACCESS='SEQUENTIAL',IOSTAT=ierr) 
    IF (ierr>0) STOP 'reading funfile error, STOP!'
    READ(UNIT=1) kfmt
    IF ( kfmt == 'VFUN1.02' ) THEN
       READ(UNIT=1) itmp,rtmp,rtmp,itmp,itmp,rtmp
    ELSE IF ( kfmt == 'VFUN1.04') THEN
       READ(UNIT=1) itmp,rtmp,rtmp,itmp,itmp,itmp,itmp,itmp,rtmp
    ELSE
       STOP 'reading eigfunction error, STOP'
    END IF
    ! model_store: r, rho, vpv, vsv,
    ! vph, vsh, Qkappa, Qmu, eta
    READ(UNIT=1) n, nicfun, nocfun, ifanis,trefl,((mstore_ms(i,j),i=1,n),j=1,9)
    CLOSE(UNIT=1)

    !! REFILL INTO REAL*8:
    do i=1,n
      do j=1,9
        model_store(i,j) = dble(mstore_ms(i,j))
      enddo
    enddo

    ALLOCATE(r(n),rho(n),epsilon(n),eta(n),g(n), &
         c01(n),c02(n),c03(n),c04(n),c05(n), &
         Qkappa(n),Qmu(n),STAT=ierr)
    IF (n<=0 .OR. ierr/=0) STOP &
         'allocate error 1 in read_prem_fun, STOP'
    r=model_store(1:n,1)/l_norm
    rho=model_store(1:n,2)/d_norm
    c01=model_store(1:n,2)*model_store(1:n,3)**2/c_norm  ! C
    c02=model_store(1:n,2)*model_store(1:n,6)**2/c_norm  ! N
    c03=model_store(1:n,2)*model_store(1:n,5)**2/c_norm  ! A
    c05=model_store(1:n,2)*model_store(1:n,4)**2/c_norm  ! L
    c04=model_store(1:n,9)*(c03-2.0*c05)             ! F
    c03=c03-c02                                      ! A-N
    c02=2.0*c02                                      ! 2N
    c04=-c04                                         ! -F
    c05=-c05                                         ! -L
    Qkappa=model_store(1:n,7)                        ! Qkappa
    Qmu=model_store(1:n,8)                           ! Qmu
    
    ! assign n to global variable NR
    NR=n

    ! find major discontinuities --> (NMOHO-,NCMB+,NICB-,NOC)
    CALL find_mdis(r,model_store(1:n,6))
    
    ! calculating eplision (hydrostatic ellipticity)
    ! Consistent with (14.21) and (14.20)
    ! Note that the variables here are non-dimensional

    ALLOCATE(k(n),radau(n),s1(n),s2(n),s3(n),i_rho_r(n),i_radau_r(n),STAT=ierr)
    IF (ierr/=0) STOP 'allocate error 2 in read_prem_fun, STOP'
    
    radau=rho*r**2   
    eta(1)=0.0
    k(1)=0.0
    
    CALL intgrl(i_rho,r,1,n,rho,s1,s2,s3,i_rho_r)       ! integral of rho*r^2
    CALL intgrl(i_radau,r,1,n,radau,s1,s2,s3,i_radau_r) ! integral of rho*r^4
    WHERE (r>0.0)
       g=4.0*i_rho_r/r**2                               ! gravity
       eta=(25.0/4.0)*((1.0-(i_radau_r/(i_rho_r*r**2)))**2)-1.0  ! (14.19)
       k=eta/r**3
    END WHERE
    
    i_rho=i_rho_r(n)
    bom=PI2/(24.0*3600.0)/f_norm
    epsilon(n)=15.0*(bom**2.0)/(24.0*i_rho*(eta(n)+2.0))
    DO i=1,n-1
       ! integrand of k*r^2=(eta/r^3)*(r^2)=eta/r 
       call intgrl(exponent1,r,i,n,k,s1,s2,s3) 
       epsilon(i)=epsilon(n)*exp(-exponent1)            ! ellipticity (14.20)
    END DO
    WRITE(*,*) '  1D model, At suface g, epsilon=', g(n), epsilon(n)
    DEALLOCATE(k,radau,s1,s2,s3,model_store,i_rho_r,i_radau_r)
    RETURN           

  END SUBROUTINE read_prem_fun

  SUBROUTINE find_mdis(r,vs)

    ! find major discontinuties
    ! --> NMOHO,NCMB,NICB,NOC
    ! --> kdis
  
    USE modellayer_module, ONLY: NOC,NMOHO,NCMB,NICB,NR,Vs_moho,kdis
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: r(:),vs(:)

    INTEGER :: nn
    INTEGER :: i, i1, ihit, jhit
    REAL*8    :: eps,zero,sum1
    INTEGER :: kdis0(30)

    zero=0.0
    eps=1E-7
    kdis0=0
  
    nn=SIZE(r)
    IF (nn/=NR .OR. SIZE(vs)/=NR) STOP 'error in SUB find_mdis, STOP'
    ! check if vs is in m/s
  
    sum1=0.0
    nn=NR
    ihit=0
    jhit=0
    nn=nn-1
    DO i =1,nn
       sum1=sum1+vs(i)
       i1=i+1
       IF (real(r(i1)-r(i)) <= eps) THEN
         ! discontinuties
         ihit=ihit+1
         IF (ihit>SIZE(kdis0)) STOP 'please enlarge the size of kdis0, &
             in SUB find_mdis, STOP'
         kdis0(ihit)=i
         IF ((vs(i1)==zero).AND.(vs(i)/=zero)) THEN
           ! ICB- or NOC- (sea floor)
           jhit=jhit+1
           IF(jhit==1) THEN
             NICB=i
           ELSE IF (jhit==2) THEN
             NOC=i
           ELSE
           STOP 'more than two liquid layers in the earth model,STOP'
         ENDIF
           
       ELSE IF ((vs(i1)/=0).AND.(vs(i)==0)) THEN
         ! CMB+
         NCMB=i1
           
       ELSE
         ! other discontinuity layers
         IF ((vs(i1)<=Vs_moho) .AND. (vs(i)>=Vs_moho)) THEN
            NMOHO=i
         ENDIF
       END IF
     END IF
   END DO
   IF (sum1/NR <=1000.0) STOP &
       'vs in find_mdis has to be in unit of m/s, STOP'
  
   IF (jhit==1) NOC=NR
   WRITE(*,'(A/,4(A,1X,I5,/))') '  major discontinuties...>', &
       '  NOC-=', NOC, &
       '  NIC-=', NMOHO, &
       '  NCMB+=', NCMB, &
       '  NICB-=', NICB
  
   IF (ASSOCIATED(kdis)) THEN
     STOP 'kdis has be defined, STOP!'
   ELSE
     ALLOCATE(kdis(ihit),STAT=i1)
     IF(i1/=0) STOP 'allocate error of kdis in SUB find_mdis, STOP'
     kdis=kdis0(1:ihit)
   END IF

  END SUBROUTINE find_mdis
    
  SUBROUTINE initialize_mantle_model(is_3d,is_3dQ,model_3d,&
       is_crust,is_cmb,model_3dQ,mantle_basis,idisc,r)

    ! HL Mar 2018

    ! including 1D spherical symmetric model and heterogeneous model
    ! write out common/mantlemodel/dB2,dS2,drho,ddisc
    ! assign NR, NMOHO,NCMB,NICB in foward_module
    ! assign all variables in prem_module

    USE modellayer_module,      ONLY : NR,NMOHO,NCMB,NICB
    USE model_lay808_module,    ONLY : lay808_mantle_model,lay808_mantle_basis
!    USE mantle_qmodel_partialsf

    IMPLICIT NONE
    include "coupling.h"

    ! inputs 
    LOGICAL, INTENT(IN) :: is_3d,is_3dQ,is_crust,is_cmb
    CHARACTER(LEN=*), INTENT(IN) :: model_3d,model_3dQ
    CHARACTER(LEN=128) :: lay808_model,cmb_topo_mod 

    ! work
    INTEGER :: k,l,m,ierr
    COMPLEX*16 dB2(0:NK,0:NS,0:NS),dS2(0:NK,0:NS,0:NS),drho(0:NK,0:NS,0:NS)
    COMPLEX*16 ddisc(ND,0:NS,0:NS)
    COMPLEX*16 dqinvk(0:NK,0:NS,0:NS),dqinvm(0:NK,0:NS,0:NS)
    REAL*8, POINTER :: r(:)
    REAL*8, POINTER :: mantle_basis(:,:)
    
    ! outputs
    INTEGER, INTENT(OUT) :: idisc(ND)

    common/mantlemodel/dB2,dS2,drho,ddisc
    common/attnmodel/dqinvk,dqinvm

    ! B2=kappa/rho is the squared bulk velocity
    ! S2=mu/rho is the squared shear velocity
    ! dB2 is the relative squared bulk velocity perturbation delta(B2)/B2
    ! dS2 is the relative squared shear velocity perturbation delta(S2)/S2
    ! dqinvk is the relative bulk modulus qinv (1/Q) model
    ! dqinvm is the relative shear modulus qinv (1/Q) model

    ! read 1d model and write NR and discontinuities 
    call read_prem_fun(r)
    allocate(mantle_basis(NR,0:NK),stat=ierr)
    if (NR<=0.or.ierr/=0) STOP &
         "allocate error in SUB initialize_mantle_module, STOP"

    ! begin reading 3D structure
    dS2=CMPLX(0.0,0.0)
    dB2=CMPLX(0.0,0.0)
    drho=CMPLX(0.0,0.0)
    dqinvk=cmplx(0.0,0.0)
    dqinvm=cmplx(0.0,0.0)

    if (is_3d) then
       ! 3d elastic model
       IF (model_3D=='lay808') THEN
          lay808_model=trim(mdir)//'MODEL_LAY808_SP12RTS/'
          cmb_topo_mod=trim(mdir)//'cmb/'
          CALL lay808_mantle_model(is_crust,is_cmb,idisc,lay808_model,cmb_topo_mod)
          CALL lay808_mantle_basis(r,mantle_basis)
          WRITE(*,*) ' mantle model: lay808 cases'
       else
          STOP 'model_3d is not defined in SUB initial_mantle_model, STOP!'
       endif
       
      ! 3d anelastic model
      if (is_3dQ) then
        if (model_3dQ=='partialsQ') then
          ! same radial basis
         ! call partialsf_mantle_qmodel
          write(*,*) ' Q model: S20RTSQ -- still to be completed'
        else
          STOP 'model_3dQ is not defined in SUB initial_mantle_model, STOP!'   
        endif
     else
       ! just 1D anelasticity
       dqinvk=cmplx(0.0,0.0)
       dqinvm=cmplx(0.0,0.0)
    endif

    if (ND>0 .and. .not.(is_crust)) then
       ddisc(1,:,:)=CMPLX(0.0,0.0)
       idisc(1)=NMOHO
    endif
    if (ND>0 .and. .not.(is_cmb)) THEN
      ddisc(2,:,:)=cmplx(0.0,0.0)
       idisc(2)=NCMB
    endif

    ! check discontinutiy layers
     if (ND>0 .and. (any(idisc<=0).or. any(idisc>NR))) then
        STOP 'idisc has to be defined in 3d model, STOP'
     endif
  else
       dS2=cmplx(0.0,0.0)
       dB2=cmplx(0.0,0.0)
       drho=cmplx(0.0,0.0)
       dqinvk=cmplx(0.0,0.0)
       dqinvm=cmplx(0.0,0.0)
    endif

  END SUBROUTINE initialize_mantle_model

  SUBROUTINE allocate_3D_structure(is_3D,is_rotate,smax)

    ! this routine allocates the arrays used to store the model parameters
    ! MD: smax is now an output : for allocating kernels
  
    USE modellayer_module, ONLY : NR,NOC,kdis
    implicit none
    include 'coupling.h'

    ! inputs
    logical, intent(in) :: is_rotate
    logical, intent(in) :: is_3d

    !output
    integer, intent(out) :: smax

    ! local variables
    integer :: ndis

    ! set the degree of the model
    smax = 0
    if(is_rotate) smax = 2
    if(is_3D) smax = NS
  
    ! set the boundary arrays
    nbnd = size(kdis)+1
    if(allocated(bnd)) deallocate(bnd)
    allocate(bnd(nbnd))
    bnd(1:nbnd-1) = kdis(:)
    bnd(nbnd) = NR
   
    ! allocate the coefficient arrays
    if(allocated(c1)) deallocate(c1)
    allocate(c1(NR))
    if(allocated (c2)) deallocate(c2)
    allocate(c2(NR))
    if(allocated(c3)) deallocate(c3)
    allocate(c3(NR))
    if(allocated(c4)) deallocate(c4)
    allocate(c4(NR))
    if(allocated(c5)) deallocate(c5)
    allocate(c5(NR))
    if(allocated(rho_st)) deallocate(rho_st)
    allocate(rho_st(NR,0:smax,-smax:smax))
    if(allocated(mu_st)) deallocate(mu_st)
    allocate(mu_st(NR,0:smax,-smax:smax))
    if(allocated(kap_st)) deallocate(kap_st)
    allocate(kap_st(NR,0:smax,-smax:smax))
    if(allocated(dmu_st)) deallocate(dmu_st)
    allocate(dmu_st(NR,0:smax,-smax:smax))
    if(allocated(dkap_st)) deallocate(dkap_st)
    allocate(dkap_st(NR,0:smax,-smax:smax))
    if(allocated(d_st)) deallocate(d_st)
    allocate(d_st(nbnd,0:smax,-smax:smax))

    rho_st = 0.d0
    mu_st = 0.d0
    kap_st = 0.d0
    dmu_st = 0.d0
    dkap_st = 0.d0
    d_st = 0.d0
  
    if(allocated(rho_el)) deallocate(rho_el)
    allocate(rho_el(NR))
    if(allocated(kap_el)) deallocate(kap_el)
    allocate(kap_el(NR))
    if(allocated( mu_el)) deallocate( mu_el)
    allocate( mu_el(NR))
    if(allocated(  d_el)) deallocate(  d_el)
    allocate(  d_el(nbnd))

    rho_el = 0.d0
    kap_el = 0.d0
    mu_el = 0.d0
    d_el = 0.d0

    ! Added for dkap, dmu
    if(allocated(dkap)) deallocate(dkap)
    if(allocated(dmu)) deallocate(dmu)
    allocate(dmu(NR))
    allocate(dkap(NR))
  
    return

  END SUBROUTINE allocate_3D_structure

  SUBROUTINE deallocate_3D_structure()

    ! this routine allocates the arrays used to store the model parameters
    ! MD: smax is now an output : for allocating kernels
  
    USE modellayer_module, ONLY : NR,NOC,kdis
    implicit none
    include 'coupling.h'

    ! local variables
    integer :: ndis

    deallocate(rho_st)
    deallocate(mu_st)
    deallocate(kap_st)
    deallocate(dmu_st)
    deallocate(dkap_st)
    deallocate(d_st)
  
    return

  END SUBROUTINE deallocate_3D_structure
  

  SUBROUTINE set_3D_structure(is_3d,is_3dQ,is_rotate,&
                              is_crust,is_CMB,mantle_basis,idisc)

  ! this routine sets the values for the material parameters of the model
  ! at the given frequency

  USE futil_module, ONLY : intgrl,thrj,deriv
  USE prem_module
  USE modellayer_module, ONLY : NR,NOC,kdis
  USE eigfun_module, ONLY : r=>eigfun_radius,eiginfo
  USE anelastic_module
  
  IMPLICIT NONE
  include 'coupling.h'
  ! inputs
  logical, intent(in) :: is_rotate
  logical, intent(in) :: is_3d
  logical, intent(in) :: is_3dQ 
  logical, intent(in) :: is_crust
  logical, intent(in) :: is_cmb
  REAL*8, POINTER :: mantle_basis(:,:)
  integer, dimension(ND), intent(in) :: idisc

  ! local variables
  integer :: smax,i,j,s,t,ndis,ibnd,id
  real*8, dimension(NR) :: fctk,fctm
  real*8, dimension(NR) :: f,fp,s1,s2,s3
  complex*16, parameter :: ii = cmplx(0.,1.)
  complex*16, dimension(NR) :: mu,kap
  
  ! common data
  COMPLEX*16 :: dB2(0:NK,0:NS,0:NS),dS2(0:NK,0:NS,0:NS),drho(0:NK,0:NS,0:NS)
  COMPLEX*16 :: ddisc(ND,0:NS,0:NS)
  common/mantlemodel/dB2,dS2,drho,ddisc
  
  ! set the degree of the model
  smax = 0
  if(is_rotate) smax = 2
  if(is_3D) smax = NS

  ! store the model degree
  sst = smax

  ! initialise the coefficients
  rho_st  = 0.d0
  mu_st   = 0.d0
  kap_st  = 0.d0
  dmu_st  = 0.d0
  dkap_st = 0.d0
  d_st    = 0.d0

  
  !====================================================================!
  !                               1D structure                         !
  !====================================================================!

  !------------------------------------------------!
  !       add in the 1D (visco)elastic modulii     !
  !------------------------------------------------!
 
  ! set the refererence modulii
  kap = (4.0*c03+c01-4.0*c04)/9.0
  mu  = (c03+3.0*c02+c01+2.0*c04-6.0*c05)/15.0
  c1 = c01
  c2 = c02
  c3 = c03
  c4 = c04
  c5 = c05

 
  ! HL AA 
  ! make anelastic corrections if needed
  ! add factors as prefactors to the frequency dependent term
  ! see eq. 6.109
  fctk = 2./(pi*Qkappa)
  fctm = 2./(pi*Qmu)
  dkap = kap*fctk
  dmu = mu*fctm

  !====================================================================!
  !                         ellipsoidal structure                      !
  !====================================================================!
  
  if(is_rotate) then
     
    ! density perturbation
    ndis = size(kdis)
    f = rho
    call deriv(f,fp,NR,r,ndis,kdis,s1,s2,s3)
    rho_st(:,2,0) = rho_st(:,2,0) + (2./3.)*sqrt(4.*pi/5.)*r*epsilon*fp
     
    ! mu perturbation
    f = real(mu)
    call deriv(f,fp,NR,r,ndis,kdis,s1,s2,s3)
    mu_st(:,2,0) = mu_st(:,2,0) + (2./3.)*sqrt(4.*pi/5.)*r*epsilon*fp
    f = imag(mu)
    call deriv(f,fp,NR,r,ndis,kdis,s1,s2,s3)
    mu_st(:,2,0) = mu_st(:,2,0) + ii*(2./3.)*sqrt(4.*pi/5.)*r*epsilon*fp
     
    ! kappa perturbation
    f = real(kap)
    call deriv(f,fp,NR,r,ndis,kdis,s1,s2,s3)
    kap_st(:,2,0) = kap_st(:,2,0) + (2./3.)*sqrt(4.*pi/5.)*r*epsilon*fp
    f = imag(kap)
    call deriv(f,fp,NR,r,ndis,kdis,s1,s2,s3)
    kap_st(:,2,0) = kap_st(:,2,0) + ii*(2./3.)*sqrt(4.*pi/5.)*r*epsilon*fp

    ! boundary perturbations
    do ibnd = 1,nbnd
      ! get lower index of the discontinuity
      i = bnd(ibnd)
      ! set the degree two topography
      d_st(ibnd,2,0) = -(2./3.)*sqrt(4*pi/5.)*r(i)*epsilon(i)
    end do

    ! store the elliptical terms for reference
    rho_el = rho_st(:,2,0)
    kap_el = kap_st(:,2,0)
    mu_el =  mu_st(:,2,0)
    d_el =   d_st(:,2,0)
     
  end if

  !====================================================================!
  !                             3D structure                           !
  !====================================================================!

  if(is_3d) then
    !-------------------------------------------------!
    !         deal with volumetic perturbations       !
    !-------------------------------------------------!
    
    ! loop over radial basis functions
    
    do j = 0,NK
      ! loop over radial nodes
       do i = 1,nr
        do s = 0,smax
          do t = 0,s
            rho_st(i,s,t)  = rho_st(i,s,t) + rho(i)*drho(j,s,t)*mantle_basis(i,j)
            mu_st(i,s,t)   = mu_st(i,s,t)  +  mu(i)*(dS2(j,s,t) + &
                 drho(j,s,t))*mantle_basis(i,j)
            kap_st(i,s,t)  = kap_st(i,s,t) + kap(i)*(dB2(j,s,t) + &
                             drho(j,s,t))*mantle_basis(i,j)
            if(t > 0) then
              rho_st(i,s,-t)  = rho_st(i,s,-t) + (-1)**t*rho(i)*conjg(drho(j,s,t))*&
                                mantle_basis(i,j)
              mu_st(i,s,-t)   = mu_st(i,s,-t)  +  (-1)**t*mu(i)*conjg(dS2(j,s,t)+&
                                drho(j,s,t))*mantle_basis(i,j)
              kap_st(i,s,-t)  = kap_st(i,s,-t) + (-1)**t*kap(i)*conjg(dB2(j,s,t)+&
                                drho(j,s,t))*mantle_basis(i,j)
            end if
              
          end do
        end do
      end do
      ! end loop over radial nodes
    end do
    ! end loop over radial basis functions
  end if

    !-------------------------------------------------!
    !         deal with boundary perturbations        !
    !-------------------------------------------------!

    if(is_crust .or. is_cmb) then
      ! deal with the crust
      do id = 1,ND
        i = idisc(id)
        do j = 1,nbnd
          if(id == 1) then
            if(i == bnd(j)) then
              ibnd = j
              exit
            end if
          else
            if(i-1 == bnd(j)) then
              ibnd = j
              exit
            end if
          end if
        end do
        do s = 0,smax
          do t = 0,smax
            d_st(ibnd,s,t) = d_st(ibnd,s,t) + ddisc(id,s,t)
            if(t > 0) then
              d_st(ibnd,s,-t) = d_st(ibnd,s,-t) + (-1)**t*conjg(ddisc(id,s,t))
            end if
          end do
        end do
      end do
    end if

    !-------------------------------------------------!
    !             deal with 3d anelasticity           !
    !-------------------------------------------------!

    if (is_3dQ) then
      do j = 0,NK
        ! loop over radial nodes
          do i = 1,nr
            do s = 0,smax
              do t = -s,s
                dmu_st(i,s,t)  = mu_st(i,s,t)*2./(pi*Qmu(i))
                dkap_st(i,s,t) = kap_st(i,s,t)*2./(pi*Qkappa(i))
              end do
            end do
          end do
          ! end loop over radial nodes
        end do
      else
    dmu_st(:,:,:) = cmplx(0.0,0.0)
    dkap_st(:,:,:) = cmplx(0.0,0.0)
  end if

  return

END SUBROUTINE set_3D_structure
  
END MODULE mantle_model_module_DS
