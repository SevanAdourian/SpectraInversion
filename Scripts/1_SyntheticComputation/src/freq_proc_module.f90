MODULE freq_proc_module

  ! ================================================================ !
  USE pchip_module

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: process, inv_fourier, hann

  ! ================================================================ !

CONTAINS

  ! ================================================================ !
  ! TEST TEST TEST PROCESSING
!   SUBROUTINE process_new(points,nt_solved,nt_orig,i1,i2,&
!        f1,f2,df0,df,t1,t2,dt,ep,aclraw,&
!        nfout,aclw,ntout,ts,aclt)
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!     ! filter raw acceleration to filtered spectra and   !
!     ! time series                                       !
!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
 
!     IMPLICIT NONE

!     ! input/output:
!     integer, intent(in) :: nt_solved,nt_orig,i1,i2
!     real*8, intent(in) :: f1,f2,df0,df,t1,t2,dt,ep
!     complex*16, dimension(:), intent(in) :: aclraw
!     complex*16, dimension(:), allocatable, intent(out) :: aclw
!     real*8, dimension(:), allocatable, intent(out) :: ts,aclt
!     integer, intent(out) :: nfout,ntout
!     real*8, dimension(:), intent(in) :: points
!     ! real*8, dimension(:), allocatable :: fs            

!     ! internal:
!     integer(kind=8) :: nt0,nt_,ii1,ii2,ratio
!     real*8, parameter :: facf = 0.1d0, fact = 0.5d0
!     integer :: i,j,nac,ne,i0
!     real*8 :: f12,f21,dff
!     real*8 :: t_,t12,t21
!     real*8 :: freq,fcoef
!     complex*16, dimension(:), allocatable :: acl_xl,atmp

!     !-------------------------------
!     logical :: interp = .false.
!     real*8, dimension(:,:), allocatable :: ri_aclw, deriv, work, imag_aclw
!     real*8, dimension(:), allocatable :: freqq,r_aclwq,i_aclwq
!     integer :: err

!     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
!     ! shift the acl into a larger array
!     allocate(acl_xl(nt_orig))
!     acl_xl(:) = dcmplx(0.d0,0.d0)
!     acl_xl(i1:i2) = aclraw
!     f12 = f1+facf*(f2-f1)
!     ! f21 = f2+facf*(f2-f1)
!     f21 = f2-facf*(f2-f1)

!     ! filter spectrum
!     do i=1,nt_orig/2+1
!       freq = (i-1)*df
!       fcoef = hann(freq,f1,f12,f21,f2)
!       acl_xl(i) = fcoef*acl_xl(i)
!     enddo

!     ! do negative frequencies
!     j = 0
!     do i=nt_orig/2+2,nt_orig
!       j = j+1
!       acl_xl(i) = conjg(acl_xl(nt_orig/2-j+1))
!     enddo

!     call inv_fourier(acl_xl,1)
!     acl_xl = acl_xl / (dt*nt_orig)
!     ! acl_xl = acl_xl / (nt_orig)

!     ! SA
!     ! Don't undo exponential decay on the time series
!     ! and save for external output
!     allocate(ts(nt_orig),aclt(nt_orig)) 
!     ts(:)   = 0.d0
!     aclt(:) = 0.d0
!     do i=1,nt_orig
!       t_ = (i-1)*dt
!       if (t_>t2) exit
!       ! acl_xl(i) = acl_xl(i)*exp(+ep*t_)
!       ts(i)   = t_
!       aclt(i) = acl_xl(i)
!     enddo
!     ntout = i-1

!     ! parameters for time domain filter
!     t12 = t1+fact*(t2-t1)
!     t21 = t2-fact*(t2-t1)

!     do i=1,nt_orig
!       t_ = (i-1)*dt
!       fcoef = hann(t_,t1,t12,t21,t2)
!       acl_xl(i) = fcoef*acl_xl(i)
!     enddo

!     ! pad time series if needed
!     nt0 = floor(1.0d0/(df0*dt)) !HL: changed from floor(1.d0/df0*dt)
!     if (nt0>nt_orig) then
!       ne  = log(real(nt0))/log(dble(2.0)) + 1
!       nt0 = 2**ne
!       allocate(atmp(nt_orig))
!       atmp = acl_xl
!       deallocate(acl_xl)
!       allocate(acl_xl(nt0))
!       acl_xl(1:nt_orig) = atmp(:)
!       acl_xl(nt_orig+1:nt0) = dcmplx(0.d0,0.d0)
!       nt_ = nt0
!       dff = dble(1.d0)/(nt_*dt)
!       write(*,*) 'do i get here?'
!       ii1  = max(floor(f1/dff),1)
!       ii2  = floor(f2/dff)+1
!       deallocate(atmp)
!    else
!       dff  = df ! HL added
!       ii1  = max(floor(f1/dff),1)
!       ii2  = floor(f2/dff)+1
!    endif

!     ! do fourier transform
!     call inv_fourier(acl_xl,-1)

!     ! HL: why is the dt factor here?
!     acl_xl = acl_xl*dt



! !!! --- HL comment out:    start
! !
! !##############################################
! !
! !    ! !---- new version
! !
! !    ! ! output filtered acceleration
! !    nac = ii2 - ii1 + 1
! !    allocate(fs(nac),aclw(nac))
! !    i0 = 1
! !    do i=ii1,ii2
! !      freq = (i-1)*dff*1000.d0
! !      fs(i0) = freq
! !      aclw(i0) = acl_xl(i)
! !      i0 = i0+1
! !    enddo
! !
! !    !nfout = nac
! !
! !    ! ! MD 
! !    ! ! setup and call interpolator
! !
! !    deallocate(acl_xl)
! !
! !    ! setup for real part
! !    allocate(ri_aclw(1,nac))
! !
! !    ri_aclw(1,:) = real(aclw(:)) ! real part
! !    ! tets case : ri_aclw ==> sin(t)
! !
! !    allocate(deriv(1,nac))
! !    call DPCHIM(nac,fs,ri_aclw,deriv,1,err)
! !    allocate(r_aclwq(nt_solved))
! !    call DPCHFE(nac,fs,ri_aclw,deriv,1,interp,nt_solved,points,r_aclwq,err) ! do interp
! !    ri_aclw(1,:) = aimag(aclw(:)) !imag part
! !    deallocate(aclw)
! !
! !    call DPCHIM(nac,fs,ri_aclw,deriv,1,err)
! !    allocate(i_aclwq(nt_solved))
! !    call DPCHFE(nac,fs,ri_aclw,deriv,1,interp,nt_solved,points,i_aclwq,err) ! do interp
! !
! !    deallocate(deriv)
! !    allocate(aclw(nt_solved))
! !    aclw = dcmplx(r_aclwq,i_aclwq)
! !    nfout = nt_solved 
! !
! !    !---- new version
! !    
! !!################################################
! !
! !
! !!! --- HL comment out:    end
    
!     ! !! Old downsampling -- HL: if necessary!?
!     ratio = df/dff
!     ! nfout = ((ii2 - ii1) / ratio)  + 1
!     nfout = ((ii2 - ii1) / ratio) 
!     ! allocate(fs(nfout),aclw(nfout))
!     allocate(aclw(nfout))
!     ! downsample
!     i0=1
!     ! do i = ii1,ii2,ratio
!     do i = ii1,ii2-1,ratio 
!        freq = (i-1)*dff*1000.d0
!        ! fs(i0) = freq
!        aclw(i0) = acl_xl(i)
!        i0 = i0 + 1
!     end do
!     nfout = nt_solved

! ! ####################################


!     !write(*,*) 'nfout = ',nfout

! !##### for now
!     ! output filtered acceleration
!     ! nac = ii2 - ii1 + 1
!     ! allocate(fs(nac),aclw(nac))
!     ! i0 = 1
!     ! do i=ii1,ii2
!     !   freq = (i-1)*dff*1000.d0
!     !   fs(i0) = freq
!     !   aclw(i0) = acl_xl(i)
!     !   i0 = i0+1
!     ! enddo
!     ! nfout = nac

!     return

! END SUBROUTINE process_adjoint_forces

  
  ! ================================================================ !

  SUBROUTINE process(points,nt_solved,nt_orig,i1,i2,&
                     f1,f2,df0,df,t1,t2,dt,ep,aclraw,&
                     nfout,aclw,ntout,ts,aclt)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! filter raw acceleration to filtered spectra and   !
    ! time series                                       !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
 
    IMPLICIT NONE

    ! input/output:
    integer, intent(in) :: nt_solved,nt_orig,i1,i2
    real*8, intent(in) :: f1,f2,df0,df,t1,t2,dt,ep
    complex*16, dimension(:), intent(in) :: aclraw
    complex*16, dimension(:), allocatable, intent(out) :: aclw
    real*8, dimension(:), allocatable, intent(out) :: ts,aclt
    integer, intent(out) :: nfout,ntout
    real*8, dimension(:), intent(in) :: points
    ! real*8, dimension(:), allocatable :: fs            

    ! internal:
    integer(kind=8) :: nt0,nt_,ii1,ii2,ratio
    real*8, parameter :: facf = 0.1d0, fact = 0.5d0
    integer :: i,j,nac,ne,i0
    real*8 :: f12,f21,dff
    real*8 :: t_,t12,t21
    real*8 :: freq,fcoef
    complex*16, dimension(:), allocatable :: acl_xl,atmp

    !-------------------------------
    logical :: interp = .false.
    real*8, dimension(:,:), allocatable :: ri_aclw, deriv, work, imag_aclw
    real*8, dimension(:), allocatable :: freqq,r_aclwq,i_aclwq
    integer :: err

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! shift the acl into a larger array
    allocate(acl_xl(nt_orig))
    acl_xl(:) = dcmplx(0.d0,0.d0)
    acl_xl(i1:i2) = aclraw
    f12 = f1+facf*(f2-f1)
    ! f21 = f2+facf*(f2-f1)
    f21 = f2-facf*(f2-f1)


    
    ! filter spectrum
    do i=1,nt_orig/2+1
      freq = (i-1)*df
      fcoef = hann(freq,f1,f12,f21,f2)
      acl_xl(i) = fcoef*acl_xl(i)
    enddo

    ! do negative frequencies
    j = 0
    ! write(*,*) 'check',i1,i2
    ! write(*,*) 'PRE'
    ! acl_xl(1) = cmplx(1,1)
    ! acl_xl(2) = cmplx(2,2)
    ! acl_xl(3) = cmplx(3,3)
    
    ! write(*,*) acl_xl(1)
    ! write(*,*) acl_xl(2)
    ! write(*,*) acl_xl(3)

    
    do i=nt_orig/2+2,nt_orig
      j = j+1
      ! acl_xl(i) = conjg(acl_xl(nt_orig/2-j+1))
      acl_xl(i) = conjg(acl_xl(nt_orig/2-j+1))
   enddo
   ! write(*,*) 'POST'
   ! write(*,*) acl_xl(nt_orig)
   ! write(*,*) acl_xl(nt_orig-1)
   ! write(*,*) acl_xl(nt_orig-2)

   call inv_fourier(acl_xl,1)
   acl_xl = acl_xl / (dt*nt_orig)
   ! acl_xl = acl_xl / (nt_orig)
   
   ! undo exponential decay on the time series
   ! and save for external output
   allocate(ts(nt_orig),aclt(nt_orig)) 
   ts(:)   = 0.d0
   aclt(:) = 0.d0
   do i=1,nt_orig
      t_ = (i-1)*dt
      if (t_>t2) exit
      ! acl_xl(i) = acl_xl(i)*exp(+ep*t_)
      ts(i)   = t_
      aclt(i) = acl_xl(i)
   enddo
   ntout = i-1

    ! parameters for time domain filter
    t12 = t1+fact*(t2-t1)
    t21 = t2-fact*(t2-t1)

    do i=1,nt_orig
      t_ = (i-1)*dt
      fcoef = hann(t_,t1,t12,t21,t2)
      acl_xl(i) = fcoef*acl_xl(i)
    enddo

    ! pad time series if needed
    nt0 = floor(1.0d0/(df0*dt)) !HL: changed from floor(1.d0/df0*dt)
    if (nt0>nt_orig) then
      ne  = log(real(nt0))/log(dble(2.0)) + 1
      nt0 = 2**ne
      allocate(atmp(nt_orig))
      atmp = acl_xl
      deallocate(acl_xl)
      allocate(acl_xl(nt0))
      acl_xl(1:nt_orig) = atmp(:)
      acl_xl(nt_orig+1:nt0) = dcmplx(0.d0,0.d0)
      nt_ = nt0
      dff = dble(1.d0)/(nt_*dt)
      write(*,*) 'do i get here?'
      ii1  = max(floor(f1/dff),1)
      ii2  = floor(f2/dff)+1
      deallocate(atmp)
   else
      dff  = df ! HL added
      ii1  = max(floor(f1/dff),1)
      ii2  = floor(f2/dff)+1
   endif


    ! write(*,*) 'ne = ',ne
    ! write(*,*)'df, dff , ii1, ii2 ', df, dff, ii1, ii2
    ! write(*,*) 'ratio df / dff = ', df/dff

    ! do fourier transform
    call inv_fourier(acl_xl,-1)

    ! HL: why is the dt factor here?
    acl_xl = acl_xl*dt
    ! acl_xl = acl_xl



!!! --- HL comment out:    start
!
!##############################################
!
!    ! !---- new version
!
!    ! ! output filtered acceleration
!    nac = ii2 - ii1 + 1
!    allocate(fs(nac),aclw(nac))
!    i0 = 1
!    do i=ii1,ii2
!      freq = (i-1)*dff*1000.d0
!      fs(i0) = freq
!      aclw(i0) = acl_xl(i)
!      i0 = i0+1
!    enddo
!
!    !nfout = nac
!
!    ! ! MD 
!    ! ! setup and call interpolator
!
!    deallocate(acl_xl)
!
!    ! setup for real part
!    allocate(ri_aclw(1,nac))
!
!    ri_aclw(1,:) = real(aclw(:)) ! real part
!    ! tets case : ri_aclw ==> sin(t)
!
!    allocate(deriv(1,nac))
!    call DPCHIM(nac,fs,ri_aclw,deriv,1,err)
!    allocate(r_aclwq(nt_solved))
!    call DPCHFE(nac,fs,ri_aclw,deriv,1,interp,nt_solved,points,r_aclwq,err) ! do interp
!    ri_aclw(1,:) = aimag(aclw(:)) !imag part
!    deallocate(aclw)
!
!    call DPCHIM(nac,fs,ri_aclw,deriv,1,err)
!    allocate(i_aclwq(nt_solved))
!    call DPCHFE(nac,fs,ri_aclw,deriv,1,interp,nt_solved,points,i_aclwq,err) ! do interp
!
!    deallocate(deriv)
!    allocate(aclw(nt_solved))
!    aclw = dcmplx(r_aclwq,i_aclwq)
!    nfout = nt_solved 
!
!    !---- new version
!    
!!################################################
!
!
!!! --- HL comment out:    end
    
    ! !! Old downsampling -- HL: if necessary!?
    ratio = df/dff
    ! nfout = ((ii2 - ii1) / ratio)  + 1
    nfout = ((ii2 - ii1) / ratio) + 1
    ! allocate(fs(nfout),aclw(nfout))
    allocate(aclw(nfout))
    ! downsample
    i0=1
    ! do i = ii1,ii2,ratio
    do i = ii1,ii2-1,ratio 
       freq = (i-1)*dff*1000.d0
       ! fs(i0) = freq
       aclw(i0) = acl_xl(i)
       i0 = i0 + 1
    end do
    nfout = nt_solved

! ####################################


    !write(*,*) 'nfout = ',nfout

!##### for now
    ! output filtered acceleration
    ! nac = ii2 - ii1 + 1
    ! allocate(fs(nac),aclw(nac))
    ! i0 = 1
    ! do i=ii1,ii2
    !   freq = (i-1)*dff*1000.d0
    !   fs(i0) = freq
    !   aclw(i0) = acl_xl(i)
    !   i0 = i0+1
    ! enddo
    ! nfout = nac

    return

  END SUBROUTINE process

  ! ================================================================ !

  FUNCTION hann(t,t11,t12,t21,t22)

    USE nrtype
    IMPLICIT NONE

    logical(lgt) :: ltmp
    real*8 :: hann
    real*8, intent(in) :: t,t11,t12,t21,t22
    ltmp = (t11 == dble(0.0) .and. t12 == dble(0.0) & 
            .and. t21 == dble(0.0) .and. t22 == dble(0.0))
    if(.not.ltmp) then
      if(t < t11) then
        hann = dble(0.0)
      else if(t >= t11 .and. t < t12) then
        hann = dble((pi_d))*(t-t11)/(t12-t11)
        hann = dble(0.5)*(dble(1.0)-cos(hann))
      else if(t >= t12 .and. t < t21) then
        hann = dble(1.0)
      else if(t >= t21 .and. t < t22) then
        hann = dble((pi_d))*(t22-t)/(t22-t21)
        hann = dble(0.5)*(dble(1.0)-cos(hann))
      else if(t >= t22) then
        hann = dble(0.0)
      end if       
    end if
    
    return

  END FUNCTION hann

  ! ================================================================ !

  SUBROUTINE inv_fourier(data,isign)

    USE nrtype 
    IMPLICIT NONE
    double complex, dimension(:), intent(inout) :: data
    integer, intent(in) :: isign
    integer :: n,i,istep,j,m,mmax_f,n2
    real*8 :: theta
    double complex :: temp
    double complex :: w,wp
    double complex :: ws

    n=size(data)

    call assert(iand(n,n-1) == 0, &
          'n must b a power of 2 in inv_fourier')

    n2=n/2
    j=n2
    do i=1,n-2
      if(j > i) call swap(data(j+1),data(i+1)) ! MD rewrote swap
      m=n2
      do
        if(m < 2 .or. j < m) exit
        j=j-m
        m=m/2
      end do
      j=j+m
    end do
    mmax_f=1
    do 
      if(n <= mmax_f) exit
      istep=2*mmax_f
      theta=dble(pi_d)/(isign*mmax_f)
      wp=cmplx(dble(-2.0)*sin(dble(0.5)*theta)**2, sin(theta),kind = dpc)
      w=cmplx(dble(1.0), dble(0.0),kind = dpc)
      do m=1,mmax_f
        ws=w
        do i=m,n,istep
          j=i+mmax_f
          temp=ws*data(j)
          data(j)=data(i)-temp
          data(i)=data(i)+temp
        end do
        w=w*wp+w
      end do
      mmax_f=istep
    end do

  END SUBROUTINE inv_fourier

  ! ================================================================ !

  SUBROUTINE assert(n1,string)
    ! Assert and swap rewrote for working with inv_fourier
    IMPLICIT NONE
    character(LEN=*), intent(in) :: string
    logical, intent(in) :: n1
    if(.not. n1) then
      write(*,*) 'Error: assert', string
      stop 'Program ended in assert '
    end if
  END SUBROUTINE assert

  ! ================================================================ !

  SUBROUTINE swap(a,b)
    IMPLICIT NONE
    double complex, intent(inout) :: a,b 
    double complex  :: temp
    temp = a
    a = b
    b = temp 
  END SUBROUTINE swap

  ! ================================================================ !

END MODULE freq_proc_module
