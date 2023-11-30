program driver_forward_only

  ! ================================================================ !
  ! May 2023                                                         !
  ! ================================================================ !
  !
  ! In this test, we will perform the forward
  ! calculations only.
  ! Outputting only the raw accelerations.

  ! We do 1D (spherically symmstric +  rotation)
  ! We do 3D (SP12RTS at 0.1x,... 2.0x + rotation)
  !
  ! ================================================================ !
  ! PREAMBLE  
  ! ================================================================ !


  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
  ! External modules:                                                 
  USE nrtype
  USE read_cluster,            ONLY : read_modesetup
  USE decipher_cluster_module, ONLY : decipher_cluster
  USE read_eqk,                ONLY : read_stations, &
                                      read_cmt, &
                                      multisource, &
                                      hdr, &
                                      polarization, &
                                      receiver_polarization, &
                                      polarization_type, &
                                      multisource_type
  USE receiver_source_module,  ONLY : build_r, &
                                      build_s
  USE futil_module, ONLY : intgrl, dirac_delta ! HL-BM
  use module_util
  USE mantle_model_module_DS
  USE splitting_seismo_module
  USE disol_module
  USE freq_proc_module
  USE modellayer_module, ONLY:NR,NCMB,NMOHO,kdis
  USE model_lay808_module

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
  ! Declarations:
  IMPLICIT NONE
  include "coupling.h"
 
  ! user model and run options
  integer :: mat_op, sbuild0, sbuild1
  logical :: is_3d, is_3dQ, is_rotate
  logical :: is_cmb, is_crust
  character(len=128) :: model_3d, model_3dQ
  character(len=128) :: flabel
  
  ! mode parameters
  real*8 :: om(MMAX), gma(MMAX), Qval(MMAX)
  ! cluster parameter
  character(len=CCLUSTMAX) :: tcluster
  integer :: nmode, ndim, type(MMAX), la(MMAX), na(MMAX)

  ! forcing and frequency parameters
  real*8 :: w_re, w_re1, w_re2, w_im, dw_re
  real*8 :: f1,f2,dt,tout,df0,t1,t2,ep 
  real*8 :: wtb,df,i11,i22
  real*8 :: wr1,wr2,wr,wi,dwr,f, dwp
  real*8, dimension(:), allocatable :: wpro
  complex*16 :: w,wfor,wadj
  integer :: ntpro,nfpro

  ! output frequency dependent data
  complex*16, dimension(:), allocatable :: rec,sce
  complex*16, dimension(:,:), allocatable :: dat
  !complex*16, dimension(:,:), allocatable :: acl_raw --> HL-BM
  complex*16, dimension(:), allocatable :: wf,acl_iw
  complex*16, dimension(:), allocatable :: aclfpro
  real*8, dimension(:), allocatable :: acltpro,tpro,fpro

  ! source/station information
  integer :: station,nstation,nsources
  real*8 :: rlat,rlon,elat,elon,evdp,nu(3)
  character(len=128) :: station_file,cmt_file, filename
  character(len=2) :: file_id
  
  ! woodhouse kernels (pre=calculated)
  integer :: iwd,nmwd,swd0,swd1
  integer :: twd,nwd,lwd
  integer, dimension(:,:,:), allocatable :: iwdhs
  complex*16, dimension(:,:,:,:), allocatable :: &
                                  K_mu_wh,K_kp_wh, &
                                  K_Tro_wh,K_Vro_wh,&
                                  K_Td_wh,K_Vd_wh

  ! coupling matrices:
  complex*16, dimension(:,:), allocatable :: &
                              V_mat,W_mat,T_mat,dV_mat
  complex*16, dimension(:,:), allocatable :: a0,a1,a2,a3

  ! adjoint kernels:
  integer :: iadj,ss,tt
  integer :: process_bool, ind_freq, freq_start_ind,freq_end_ind
  real*8 :: dfk
  real*8 :: integration_band, central_frequency, thresh
  ! HL-BM: add extra dimension for adjoint
  !        kernels (frequency points -- don't integrate f)
  ! complex*16, dimension(:,:,:,:,:), allocatable :: &
  !                                 kmu,kkp,&
  !                                 kro
  complex*16, dimension(:,:,:,:), allocatable :: &
                                  kmutest,kkptest,&
                                  krotest, kds
  complex*16, dimension(:,:), allocatable :: afor,aadj
  complex*16, dimension(:,:,:,:), allocatable :: aclf_raw
  complex*16, dimension(:,:,:,:), allocatable :: acla_raw
  complex*16, dimension(:,:), allocatable :: aclf_iw,&
                                             acla_iw
  complex*16, dimension(:,:), allocatable :: vs_adj
  complex*16, dimension(:), allocatable :: aclffpro
  complex*16, dimension(:), allocatable :: aclafpro
  complex*16, dimension(:,:), allocatable :: ukf
  complex*16, dimension(:,:), allocatable :: uka
  complex*16, dimension(:,:), allocatable :: ukatest
  ! complex*16, dimension(:,:,:), allocatable :: uka
  ! complex*16, dimension(:,:,:), allocatable :: ukatest
  complex*16, dimension(:), allocatable :: acl_sol
  real*8, dimension(:), allocatable :: for_wr,adj_wr
  complex*16, dimension(:,:,:), allocatable :: acl_raw

  ! solver parameters:
  complex*16, dimension(:,:), allocatable :: a
  complex*16, dimension(:), allocatable :: wk

  ! paths:
  character(len=32), parameter :: fmtt = '(2E15.7)'
  character(len=32), parameter :: fmtf = '(4E15.7,I5)'
  character(len=32), parameter :: fmtk = '(4I5,2E18.9)'
  character(len=32) :: freefmt1,freefmt2
  character(len=128) :: f_aclt,f_aclf
  character(len=128) :: f_kmu,f_kkp,f_kro,f_kds

  ! misc:
  integer, parameter :: io1 = 7,io2=8,io3=9,io4=10
  integer :: im,i,j,info,istat,nindex,ifgot,nw,iw,ir1,&
             nt,i1,i2,mex,qex,nt0,ne,ntb,md,iter,mindex,ib,&
             jb,tbd,count,ntime,ijob,ifcor,im2,nbyts,&
             lmod,lpath,imode,s,t,ir, iwtest, ifreq
  integer, dimension(:), allocatable :: nn,ll,ity
  complex*16, dimension(:), allocatable :: vs,vsw ! HL-mat
  complex*16, dimension(:,:), allocatable :: vr
  logical :: exists
  integer :: it
  logical :: ltmp
  integer :: nt_orig
  real*8, pointer :: mantle_basis(:,:), r(:)
  integer, dimension(nd) :: idisc
  integer :: smax
  real*8 :: valr,vali
  character(len=128) :: sjunk

  ! BM23:
  integer, parameter :: nper = 2,nftest = 1
  real*8, dimension(nper) :: dper
  character*10 :: dperc
  integer, dimension(nftest) :: iftest
  real*8, dimension(:), allocatable :: tmp_array
  integer :: iper,sper,tper,ik0,ik1,ierr
  character(len=128) :: laymod,fout1
  character(len=128),dimension(nper) :: fout2
  real*8 :: ff,per
  complex*16, dimension(:,:), allocatable :: aclf_amp,adjforce

  COMPLEX*16 :: dB2(0:NK,0:NS,0:NS),dS2(0:NK,0:NS,0:NS),drho(0:NK,0:NS,0:NS)
  COMPLEX*16, dimension(0:NK,0:NS,0:NS) :: dB2_0,dS2_0,drho_0
  COMPLEX*16 :: ddisc(ND,0:NS,0:NS)
  complex*16, dimension(:,:), allocatable :: acl_raw0
  real*8, dimension(:), allocatable :: fr,fi,s1,s2,s3,rr
  real*8 :: sumr,sumi
  real*8 :: freq_for_delta
  real*8, dimension(:),allocatable :: du_sum
  complex*16, dimension(:,:),allocatable :: du
  complex*16, dimension(:), allocatable :: dfull,ratio
  complex*16, dimension(:,:), allocatable :: dfullout
 
  complex*16, dimension(:), allocatable :: spro1d
  complex*16, dimension(:,:), allocatable :: spro3d
  complex*16, dimension(:,:,:), allocatable :: hforce
  complex*16, dimension(:), allocatable :: hforw
  complex*16, dimension(:), allocatable :: htest
  complex*16, dimension(:,:), allocatable :: delta_s
  complex*16, dimension(:), allocatable :: dirac_tmp, dirac_pro_tmp
  complex*16, dimension(:,:), allocatable :: dirac, dirac_pro
  real*8 :: real_dir, imag_dir
  
  NAMELIST/forward_inp/&
       is_3d,&
       is_3dQ,&
       is_rotate,&
       is_crust,&
       is_cmb,&
       model_3d,&
       model_3dQ,&
       mat_op,&
       sbuild0,&
       sbuild1

  NAMELIST/freq_param/&
       f1,f2,dt,df0,&
       wtb,t1,t2,central_frequency,integration_band,&
       process_bool,&
       station_file,&
       cmt_file


  common/mantlemodel/dB2,dS2,drho,ddisc

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
  ! end of preamble                                                  !
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

  ! ================================================================ !
  ! BEGIN PROGRAM
  ! ================================================================ !

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
  ! Preliminaries
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

  print *, '===================================================='
  print *, '               BEGIN CALCULATIONS                   '
  print *, '===================================================='
  print *, ''
  print *, ''

  write(*,*) 'INITIAL SET UP'

    
  ! (1) some directories here:
  laymod    = trim(mdir)//'MODEL_LAY808_SP12RTS/'
  
  ! fout1     = 'for_acc_1d.dat'
  ! fout2     = (/'for_acc_3d_x1.0.dat'/)
       
       


  
  call getarg(1,flabel)
  open(unit=15,file='forward_inp.namlist'//trim(flabel),&
       action='read',status='old')
  read(unit=15,nml=forward_inp)
  close(15)

  open(unit=15,file='freq_param.namelist',&
       action='read',status='old')
  read(unit=15,nml=freq_param)
  close(15)

  !     print these out:
  print *, '===================================================='
  print *, '               Perturbations considered             '
  print *, '===================================================='
  if (is_3d)     write(*,*) '.... considering 3D structure'
  if (is_3dQ)    write(*,*) '.... considering 3D-Q structure'
  if (is_rotate) write(*,*) '.... considering rotation'
  if (is_crust)  write(*,*) '.... considering crustal structure'
  if (is_cmb)    write(*,*) '.... considering CMB topogrpahy'

  ! (2) set up a few bits:
  call read_prem_fun(r)
  ! HL-BM: set from above
  call read_modesetup('setup',tcluster)
  call decipher_cluster(tcluster,nmode,ndim,na,type,la,om,gma,Qval)
  nbnd = 13 ! hard code

  ! (3) read in precalculated Woodhouse kernels:
  open(unit=15,file=trim(sdir)//'wdhsinfo',action='read',status='old')
  read(15,*) sjunk
  read(15,*) sjunk
  read(15,*) nmwd, swd0, swd1
  print*, nmwd, swd0, swd1
  allocate(K_mu_wh(1:nmwd,1:nmwd,1:nr,swd0:swd1),&
       K_kp_wh(1:nmwd,1:nmwd,1:nr,swd0:swd1),&
       K_Tro_wh(1:nmwd,1:nmwd,1:nr,swd0:swd1),&
       K_Vro_wh(1:nmwd,1:nmwd,1:nr,swd0:swd1),&
       K_Td_wh(1:nmwd,1:nmwd,1:nbnd,swd0:swd1),&
       K_Vd_wh(1:nmwd,1:nmwd,1:nbnd,swd0:swd1))
  allocate(iwdhs(0:1,0:120,0:120))
  iwdhs(:,:,:) = 0
  do iwd=1,nmwd
     read(15,*) twd,nwd,lwd
     iwdhs(twd,nwd,lwd) = iwd
  enddo
  close(15)
  open(15,file=trim(sdir)//'wdhskernels',form='unformatted')
  read(15) K_mu_wh
  read(15) K_kp_wh  
  read(15) K_Tro_wh 
  read(15) K_Vro_wh 
  read(15) K_Td_wh 
  read(15) K_Vd_wh  
  close(15)
    
  ! (4) get source-receiver data:
  write(*,*) 'Calculate source-receiver'
  station_file = 'station.info'
  cmt_file     = 'cmt.info'
  
  call read_stations(station_file,nstation)
  write(*,*) '.... Station read, num station = ', nstation
  call read_cmt(cmt_file)
  write(*,*) '.... Moment tensor read'

  ! Scale the moment tensor
  multisource(1)%moment_tensor = multisource(1)%moment_tensor &
       * dble(0.95179e-30)
  allocate(rec(ndim),sce(ndim))
  allocate(vr(ndim,nstation),vs(ndim))
  allocate(vsw(ndim)) ! HL-BM
  
  !     compute source/receiver vectors
  do station=1,nstation
     ! parameters for station considered
     rlat = hdr(station)%stla
     rlon = hdr(station)%stlo
     nu   = receiver_polarization(polarization(station)%pola1,&
          polarization(station)%pola2,&
          polarization(station)%pola3)
     ! build receiver vector
     call build_r(rec,nmode,na,type,la,rlat,rlon,nu)
     ! store in container
     vr(:,station) = rec
  enddo ! end loop build receiver

  !     prepare the source elements
  elat = multisource(1)%elat
  elon = multisource(1)%elon
  evdp = multisource(1)%depth
  !     build the source vector
  call build_s(sce,nmode,na,type,la,elat,elon,evdp,&
       multisource(1)%moment_tensor)  
  vs(:) = sce(:)

  ! (5) non-dimensionalize frequency/time parameters:
  f1  = f1/1000.d0  /f_norm
  f2  = f2/1000.d0  /f_norm
  df0 = df0/1000.d0 /f_norm
  wtb = 2.*pi*wtb/1000.d0 /f_norm
  t1  = t1*3600.d0  /t_norm 
  t2  = t2*3600.d0  /t_norm 
  ! dt  = dt /t_norm
  dt = 1.d0 / (2.d0*f2)
  print*, dt
  
  ! non dimensionalisation necessary for dirac delta
  ! fixed values of mex and qex from David's thesis:
  mex = 5
  qex = 4
  call fcal_new(f1,f2,dt,t2,mex,qex,df,ep,nt,i1,i2)
  wr1     = pi2*f1
  wr2     = pi2*f2
  dwr     = pi2*df
  wi      = ep
  nt_orig = nt      ! store full number of time steps
  nt      = i2-i1+1 ! the number of steps to be iterated
                    ! considering the band  
  
  print*, f1, f2, df*f_norm*1000.d0, nt, dwr, dt

  ! (6) set radial basis:
  allocate(mantle_basis(NR,0:NK),stat=ierr)
  mantle_basis(:,:) = 0.0
  do i=0,NK
     j = NCMB + i
     mantle_basis(j,i) = 1.
  enddo
   
  ! (7) Set perturbations and *final* frequencies to compare
  ! dper   = (/0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0/)
  dper   = (/0.0,1.0/)

  ! (8) Get mode info:
  allocate(wk(nmode),ll(nmode))
  do iw=1,nmode
     wk(iw) = dcmplx(om(iw),gma(iw))
     ll(iw) = la(iw)
     write(*,*) ll(iw),wk(iw)*1000.*f_norm/pi2
     ! write(*,*) realpart(wk(iw))*1000.d0*f_norm/pi2
  enddo

  ! stop
  ! ==================================================== !
  ! BEGIN LOOPING OVER                                   !
  ! ==================================================== !

  ! ================================= !
  ! STEP 1 (in preamble)              !
  ! ================================= !
  ! Calculate all forward spectra:
  
  ! allocate and initialize:
  allocate(a(ndim,ndim),afor(ndim,ndim),aadj(ndim,ndim))
  allocate(a0(ndim,ndim),a1(ndim,ndim))
  allocate(a2(ndim,ndim),a3(ndim,ndim))
  allocate(acl_raw0(nt,nstation),acl_raw(nper,nt,nstation))
  allocate(acl_iw(nstation))
  allocate(wf(nt),for_wr(nt))
    
  acl_raw0(:,:)  = dcmplx(0.d0,0.d0)
  acl_raw(:,:,:) = dcmplx(0.d0,0.d0)
  wf(:)          = dcmplx(0.d0,0.d0)
  for_wr(:)      = 0.d0

  ! 1D problem:
  ! call allocate_3D_structure(.false.,is_rotate,smax)
  ! call set_3D_structure(.false.,is_3dQ,is_rotate,&
  !      is_crust,is_CMB,mantle_basis,idisc)
  ! call coupling_matrix(0,2,is_rotate,ndim,nmode,&
  !      swd0,swd1,nmwd,iwdhs,&
  !      K_mu_wh,K_kp_wh,K_Tro_wh,K_Vro_wh,&
  !      K_Td_wh,K_Vd_wh,&
  !      a0,a1,a2,a3)
  ! call deallocate_3D_structure()

  ! write(*,*) 'BEGIN CALCULATING ALL FORWARD SPECTRA'
  ! do iw=1,nt
  !    write(*,*)"Processing frequency ",iw,"/",nt
  !    if (iw .eq. 1)then
  !       acl_raw0(iw,:) = dcmplx(0.d0,0.d0)
  !    else
  !       a(:,:)     = dcmplx(0.,0.)
  !       wr         = wr1 + (iw-1)*dwr
  !       wf(iw)     = wr - ii*wi
  !       for_wr(iw) = wr
  !       do i=1,ndim
  !          do j=1,ndim
  !             a(i,j) = a0(i,j) + wf(iw)*a1(i,j) &
  !                  + wf(iw)*wf(iw)*a2(i,j)
  !          enddo
  !          a(i,i) = a(i,i) - wf(iw)*wf(iw)
  !       enddo

  !       ! call solver
  !       vsw(:) = vs(:)/(ii*wf(iw))
  !       call solve(wf(iw),wtb,iw,wk,nmode,ndim,nstation,ll,a,vsw,vr,acl_iw)
  !       acl_raw0(iw,:) = -wf(iw)*wf(iw)*acl_iw
  !    endif
  ! enddo

  ! ! write(dperc,'(F4.4)') dper(iper)
  ! filename = trim(hdr(1)%name)//'.dat'
  ! open(211,file=trim(filename),action="write", status='replace')
  ! do i=1,nt
  !    write(211,*) for_wr(i)*1000/(2.*pi)*f_norm, &
  !         dble(acl_raw0(i,1)), aimag(acl_raw0(i,1)), &
  !         abs(acl_raw0(i,1))
  ! end do
  ! close(211)

  ! stop
  
  ! 3D problems:     
  ! set SP12RTS perturbations:
  call allocate_3D_structure(is_3D,is_rotate,smax)
  call lay808_mantle_model(is_crust,is_cmb,idisc,&
       laymod,laymod)
  call lay808_mantle_basis(r,mantle_basis)
  dS2_0(:,:,:)  = dS2(:,:,:)
  dB2_0(:,:,:)  = dB2(:,:,:)
  drho_0(:,:,:) = drho(:,:,:)
  dS2(:,:,:) = cmplx(0.,0.)
  dB2(:,:,:) = cmplx(0.,0.)
  drho(:,:,:) = cmplx(0.,0.)
  do iper=1,nper
     print*, iper, dper(iper)
     ! set perturbations:
     dS2  = dS2_0*dper(iper)
     dB2  = dB2_0*dper(iper)
     drho = drho_0*dper(iper)
     call set_3D_structure(is_3d,is_3dQ,is_rotate,&
          is_crust,is_cmb,mantle_basis,idisc)
     call coupling_matrix(0,sbuild1,is_rotate,ndim,nmode,&
          swd0,swd1,nmwd,iwdhs,&
          K_mu_wh,K_kp_wh,K_Tro_wh,K_Vro_wh,&
          K_Td_wh,K_Vd_wh,&
          a0,a1,a2,a3)
     do iw=1,nt
        write(*,*)"Processing frequency ",iw,"/",nt, iper,"/",nper
        if (iw .eq. 1) then
           acl_raw(iper,iw,:) = dcmplx(0.d0,0.d0)
        else
           a(:,:)     = dcmplx(0.,0.)
           wr         = wr1 + (iw-1)*dwr
           wf(iw)     = wr - ii*wi
           for_wr(iw) = wr
           do i=1,ndim
              do j=1,ndim
                 a(i,j) = a0(i,j) + wf(iw)*a1(i,j) &
                      + wf(iw)*wf(iw)*a2(i,j)
              enddo
              a(i,i) = a(i,i) - wf(iw)*wf(iw)
           enddo
           ! call solver
           vsw(:) = vs(:)/(ii*wf(iw))
           call solve(wf(iw),wtb,iw,wk,nmode,ndim,nstation,ll,a,vsw,vr,acl_iw)
           acl_raw(iper,iw,:) = -wf(iw)*wf(iw)*acl_iw
        endif
     enddo

     if (dper(iper) .eq. 0) then
        filename = trim(odir)//'1D/'//trim(hdr(1)%name)//'.dat'
     else
        filename = trim(odir)//'SP12/'//trim(hdr(1)%name)//'.dat'
     end if
     
     print*, filename
     open(211,file=trim(filename),action='write',status='replace')
     do i=1,nt
        write(211,*) for_wr(i)*1000/(2.*pi)*f_norm, &
             dble(acl_raw(iper,i,1)), aimag(acl_raw(iper,i,1)), abs(acl_raw(iper,i,1))
     end do
     close(211)
     
  enddo
  
  call deallocate_3D_structure()

  write(*,*) 'Done!'
  
end program driver_forward_only
