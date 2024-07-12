!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONSTANTS                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real*8 PI,PI2
       parameter(PI=3.14159265358979,PI2=6.28318530717958)
       real*8 w2mhz
       parameter(w2mhz=159.15494309189549983131)
       real*8 RA,GRAV,RHOAV,SCALE_FAC
       parameter(RA=6371000.0,GRAV=6.6723e-11,RHOAV=5515.0,SCALE_FAC=1.0E+10)

       ! RA = radius of Earth (m)
       ! GRAV = gravitational constant (N m^2 / kg^2)
       ! RHO_AV = average density of the Earth (kg / m^3)
       ! SCALE_FAC: scaling factor for seismograms; one needs to divide
       !            data and synthetics by SCALE_FAC to get true
       !            ground acceleration


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NORMALIZATION PARAMETERS                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real*8, parameter :: l_norm = RA                        ! length
       real*8, parameter :: d_norm = RHOAV                     ! density
       real*8, parameter :: f_norm = dsqrt(PI*GRAV*RHOAV)      ! frequency
       real*8, parameter :: m_norm = d_norm*l_norm**3          ! mass
       real*8, parameter :: t_norm = 1.0/f_norm                ! time
       real*8, parameter :: v_norm = RA/t_norm                 ! velocity
       real*8, parameter :: c_norm = d_norm*v_norm**2          ! stress
       real*8, parameter :: a_norm = v_norm/t_norm             ! acceleration
       real*8, parameter :: p_norm = a_norm*l_norm             ! gravitation potential   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODE CLUSTER PARAMETERS                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! CCLUSTMAX: maximum length of characters of cluster
      ! MMAX: maximum number of modes to be coupled
      integer CCLUSTMAX,MMAX
      parameter(CCLUSTMAX=1250,MMAX=200)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! REFERENCE MODEL                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real*8 wref
       parameter(wref=1.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D MANTLE MODEL                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! NK: number of radial degrees for mantle model
       ! NS: number of angular degrees = MAX(NS_D_MOHO,NS_[S,P]_[mantle_model]) 
       ! ND: number of discontinuities for mantle model
       integer NK,NS,ND
       parameter(NK=476,NS=20,ND=2) ! lay 808
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KERNEL OPTIONS                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! skmax: maximum structural degree to calculate
       integer skmax
       parameter(skmax=8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FILE PATHS                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! Top directory
       character(len=132), parameter :: main_dir='/data/SA/SpectraInversion/Scripts/1_SyntheticComputation/'
       ! Directory to store matrices
       character(len=132), parameter :: sdir=trim(main_dir)//'/matrices/'
       ! Output directory for results 
       character(len=132), parameter :: odir=trim(main_dir)//'/results/'
       ! Directory for models
       character(len=132), parameter :: mdir=trim(main_dir)//'/models/'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EIGENFUNCTION FILES                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! NHEAD  : number of headers in eigfunctions (see MINEOS)
       ! NRMAX  : number of maximal stored knots for eigenfunctions
       ! funXXX : filename of binary eigenfunction files
       integer, parameter :: NHEAD=6, NRMAX=1000
       character(len=128), parameter :: funTor=trim(main_dir)//&
                           '/mineos/isoprem808_elastic/toroidal_isoprem808.bin'
       character(len=128), parameter :: funRad=trim(main_dir)//&
                           '/mineos/isoprem808_elastic/radial_isoprem808.bin'
       character(len=128), parameter :: funSph=trim(main_dir)//&
                           '/mineos/isoprem808_elastic/spheroidal_isoprem808.bin'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D MODELS                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! CRUST MODELS
       integer, parameter :: NS_D_MOHO=24
       character(len=128), parameter :: fin_moho_variation=trim(main_dir)//&
                           '/MODELS/crust/crust1.0_l24_d24.4.dat'
       integer, parameter :: NS_D_CMB=6 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MISC                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, parameter :: LMAX=55 ! Largest angular degree of mode included (c.f receiver function)
      ! freq domain (fourier transform related to avoid ringing)
      integer, parameter :: dpcc=kind((1.0d00,1.0d0))
      integer, parameter :: nfmax = 10000000 ! largest frequency storage
      integer, parameter :: mpi1 = 4 ! used only with mpirun
      integer, parameter :: mpi2 = 4



