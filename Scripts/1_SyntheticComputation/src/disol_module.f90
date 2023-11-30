MODULE disol_module
  ! This module contains all the subroutines for
  ! the Direct Iterative Solution method

  IMPLICIT NONE
  include 'coupling.h'

  PRIVATE
  PUBLIC :: fcal,&
       fcal_new,&
       solve,&
       solvek_f,&
       solvek_a,&
       check_mode,&
       solvek_f2
  
  ! HL-BM:
  ! Bug was found here.
  ! All source vectors (for *all* "solve" subroutines
  ! have the "iw" removed from vs/iw.
  ! This is for the earthquake acceleration in the
  ! frequency domain. We do this outside now.
  ! see "HL-BM" in this file:
  
CONTAINS

  SUBROUTINE fcal(f1,f2,dt,tout,mex,qex,df,ep,nt,i1,i2)
    
    implicit none
    
    real*8, intent(inout) :: f1,f2
    real*8, intent(in) :: dt,tout
    integer, intent(in) :: mex,qex
    real*8, intent(out) :: df,ep
    integer, intent(out) :: nt,i1,i2
    integer :: ne
    real*8 :: fn

    ! check that the Nyquist frequency for the given
    ! value of dt lies above f2
    fn = sngl(0.5)/dt
    if(fn < f2) stop &
     ' f2 is greater than the Nyquist frequency for the time step'

    ep = mex/tout
    write(*,*) "ep=",ep
    df = ep/(pi2*qex)
    write(*,*) "df=",df
    nt = sngl(1.0)/(df*dt)
    write(*,*) "nt=",nt
    ne = log(real(nt))/log(sngl(2.0))+1
    write(*,*) "ne=",ne
    nt = 2**ne
    write(*,*) "nt2=",nt

    df = sngl(1.0)/(nt*dt)
    write(*,*) "df=",df
    ! Applying Nyquist here
    ! df = sngl(1.0)/(2*nt*dt)
    i1 = max(floor(f1/df),1)
    write(*,*) "i1=",i1
    f1 = (i1-1)*df

    i2 = floor(f2/df)+2
    write(*,*) "i2=",i2
    f2 = (i2-1)*df

    write(*,*) '.... Frequency parameters:'
    write(*,*) '.... ... df (mHz) =',df
    write(*,*) '.... ... nt       =',nt
    write(*,*) '.... ... f1 (mHz) =',f1
    write(*,*) '.... ... f2 (mHz) =',f2


    return

  END SUBROUTINE fcal

  !===========================================================!
  SUBROUTINE fcal_new(f1,f2,dt,tout,mex,qex,df,ep,nt,i1,i2)
    
    implicit none
    
    real*8, intent(inout) :: f1,f2
    real*8, intent(in) :: dt,tout
    integer, intent(in) :: mex,qex
    real*8, intent(out) :: df,ep
    integer, intent(out) :: nt,i1,i2
    integer :: ne
    real*8 :: fn

    ! check that the Nyquist frequency for the given
    ! value of dt lies above f2
    fn = sngl(0.5)/dt
    if(fn < f2) stop &
         ' f2 is greater than the Nyquist frequency for the time step'
    
    ep = mex/tout
    df = ep/(pi2*qex)
    nt = sngl(1.0)/(df*dt)
    ! ne = log(real(nt))/log(sngl(2.0))+1
    ! nt = 2**ne
    write(*,*) "ep=",ep
    write(*,*) "df=",df
    write(*,*) "nt=",nt

    ! df = sngl(1.0)/(nt*dt)
    i1 = floor(f1/df) + 1
    f1 = (i1-1)*df

    i2 = ceiling(f2/df) + 1
    f2 = (i2-1)*df

    write(*,*) "i1=",i1
    write(*,*) "i2=",i2

    write(*,*) '.... Frequency parameters:'
    write(*,*) '.... ... df (mHz) =',df
    write(*,*) '.... ... nt       =',nt
    write(*,*) '.... ... f1 (mHz) =',f1
    write(*,*) '.... ... f2 (mHz) =',f2


    return

  END SUBROUTINE fcal_new
  
  !===========================================================!
  
  SUBROUTINE solve(w,wtb,iw,wk,nmode,ndim,nr,ll,a,vs,vr,acl_iw)

    USE module_util

    implicit none
    include 'coupling.h'
    
    ! inputs:
    integer, intent(in) :: iw,nmode,ndim,nr
    integer, dimension(:), intent(in) :: ll
    real*8, intent(in) :: wtb
    complex*16, intent(in) :: w
    complex*16, dimension(:,:), intent(in) :: a
    complex*16, dimension(:,:), intent(in) :: vr
    complex*16, dimension(:), intent(in) :: wk 
    complex*16, dimension(:), intent(in) :: vs

    ! outputs:
    complex*16, dimension(:), intent(out) :: acl_iw

    ! internal:
    integer, parameter :: maxit = 10000
    real*8, parameter :: tol = 1.d-6
    integer :: i,md,im,ib,jb,tbd,info
    integer :: it,icount,ir
    integer :: ntb,nindex,mindex
    real*8 :: err
    complex*16 :: cfac,alpha,beta,tmp1,tmp2
    complex*16 :: aclr
    complex*16 :: ii = (0.d0,1.0d0)
    logical :: ltmp 
    integer, dimension(:), allocatable :: ipivot
    integer, dimension(:,:), allocatable :: tbl
    complex*16, dimension(:,:), allocatable :: atb
    complex*16, dimension(:,:), allocatable :: vsd
    complex*16, dimension(:), allocatable :: xs,xr,x0
    complex*16, dimension(:), allocatable ::  s0,sb0,p0,pb0,s1, &
                                              sb1,p1,pb1,vt

    ! allocate and initialize
    allocate(ipivot(ndim))
    allocate(tbl(nmode,4))
    allocate(atb(ndim,ndim))
    allocate(vsd(ndim,1))
    allocate(xs(ndim),xr(ndim),x0(ndim))
    allocate(s0(ndim),sb0(ndim),p0(ndim),pb0(ndim),&
             s1(ndim),sb1(ndim),p1(ndim),pb1(ndim),vt(ndim))


    ! build pre conditioning matrix
    ! determine which modes are in the target block:
    ntb    = 0
    nindex = 0
    mindex = 0
    do im = 1,nmode
      if(check_mode(w,wk(im),wtb)) then
        ntb = ntb + 1
        tbl(ntb,1) = im
        tbl(ntb,2) = nindex
        tbl(ntb,3) = mindex
        tbl(ntb,4) = 2*ll(im) + 1
        mindex = mindex + 2*ll(im) + 1
      endif
      nindex = nindex + 2*ll(im) + 1
    enddo

    ! deal with target block contribution
    if (ntb/=0) then
      ! dimension of the target system
      tbd = 0
      do ib=1,ntb
        tbd = tbd + tbl(ib,4)
      enddo

      ! build the coupling matrix for the 
      ! target block
      do ib = 1,ntb
        do jb = 1,ntb
          atb(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4), &
              tbl(jb,3)+1:tbl(jb,3)+tbl(jb,4)) = &
            a(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4), &
              tbl(jb,2)+1:tbl(jb,2)+tbl(jb,4))
        enddo
      enddo

      ! perform an LU decomposition of the target block
      call zgetrf(tbd,tbd,atb,ndim,ipivot,info)

      ! assemble source vector
      do ib = 1,ntb
        vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
          vs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) ! /(ii*w) !HL-BM 
      enddo
      
      ! solve the linear system for target block
      call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
      do ib = 1,ntb
        xs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
          vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
      end do
    endif

    ! deal with the rest of the matrix
    if (ntb /= nmode) then
      nindex = 0
      do im = 1,nmode
        md = 2*ll(im)+1
        if(.not.check_mode(w,wk(im),wtb)) then
          do i = 1,md
            cfac = 1.0d0/(wk(im)**2-w**2)
            xs(nindex+i) = cfac*vs(nindex+i)!/(ii*w) !! HL-BM
            
          end do
        end if
        nindex = nindex+md
      end do
    end if
 
    ! initial guess for the solution
    x0 = xs

    ! set the initial values for the vectors
    s0 = matmul(a,x0)
    if(ntb /= 0) then
      do ib = 1,ntb
        vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
          s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
      enddo
      call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
      do ib = 1,ntb
        s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
          vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
      enddo
    endif
    if(ntb /= nmode) then
      nindex = 0
      do im = 1,nmode
        md = 2*ll(im)+1
        ltmp = check_mode(w,wk(im),wtb)
        if(.not.ltmp) then
          do i = 1,md
            cfac = 1.0d0/(wk(im)**2-w**2)
            s0(nindex+i) = cfac*s0(nindex+i)
          enddo
        endif
        nindex = nindex+md
      enddo
    endif
    s0 = xs-s0
    sb0 = s0
    p0  = s0
    pb0 = sb0
    xr  = x0

    ! check the initial error
    tmp1 = my_dot(sb0,s0)
    err = sqrt(abs(tmp1/my_dot(xs,xs)))
    it  = 1 ! dummy initialization
    if(err > tol) then

      ! start the iteration loop
      icount = 0
      do it=1,maxit

        ! iterate the vectors
        vt = matmul(a,p0)
        if(ntb /= 0) then
          do ib = 1,ntb
            vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
              vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
          enddo
          call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
          do ib = 1,ntb
            vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
          enddo
        end if
        if(ntb /= nmode) then
          nindex = 0
          do im = 1,nmode
            md = 2*ll(im)+1
            if(.not.check_mode(w,wk(im),wtb)) then
              do i = 1,md
                cfac = 1.0d0/(wk(im)**2-w**2)
                vt(nindex+i) = cfac*vt(nindex+i)
              end do
            end if
            nindex = nindex+md
          enddo
        endif
        alpha = tmp1/my_dot(pb0,vt)
        s1 = s0-alpha*vt

        if(ntb /= 0) then
          do ib = 1,ntb
            vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
              pb0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
          enddo
          call zgetrs('T',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
          do ib = 1,ntb
            vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
          enddo
        endif
        if(ntb /= nmode) then
          nindex = 0
          do im = 1,nmode
            md = 2*ll(im)+1
            if(.not.check_mode(w,wk(im),wtb)) then
              do i = 1,md
                cfac = 1.0d0/(wk(im)**2-w**2)
                vt(nindex+i) = cfac*pb0(nindex+i)
              enddo
            endif
            nindex = nindex+md
          enddo
        endif
        vt = matmul(transpose(a),vt)

        sb1  = sb0-alpha*vt
        tmp2 = my_dot(sb1,s1)
        beta = tmp2/tmp1
        p1   = s1+beta*p0
        pb1  = sb1+beta*pb0

        ! update the solution
        xr = xr+alpha*p0

        ! estimate the error
        vt = matmul(a,xr)
        if(ntb /= 0) then
          do ib = 1,ntb
            vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
              vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
          enddo
          call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
          do ib = 1,ntb
            vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
          enddo
        endif
        if(ntb /= nmode) then
          nindex = 0
          do im = 1,nmode
            md = 2*ll(im)+1
            if(.not.check_mode(w,wk(im),wtb)) then
              do i = 1,md
                cfac = 1.0d0/(wk(im)**2-w**2)
                vt(nindex+i) = cfac*vt(nindex+i)
              enddo
            endif
            nindex = nindex+md
          enddo
        endif
        icount = it
        vt = vt-xs

        err = sqrt(abs(my_dot(vt,vt)/my_dot(xs,xs)))
        ! print *, 'iteration = ',it,' err = ',err/tol
        ! update the vectirs:
        s0 = s1
        sb0 = sb1
        p0 = p1
        pb0 = pb1
        tmp1 = tmp2       
 
        if (err<tol) then
          !print*,'Total iteration = ',it,'for iteration',iw
          exit
        endif
      enddo
    else 
       ! print *, 'no iterations',' err = ',err/tol
    end if

    if (it==maxit) then
      print*,'Max iterations at iteration',iw
    endif

    ! converged (unless max iterations hit)
    ! begin loop over the receivers

    do ir=1,nr
      ! form the product with the receiver vector
       aclr = dot_product(vr(:,ir),xr) ! BM-HL
       !aclr = my_dot(vr(:,ir),xr) ! BM - vr already conjugated
      acl_iw(ir) = aclr
    enddo

    return

  END SUBROUTINE solve

  SUBROUTINE solvek_f(w,wtb,iw,wk,nmode,ndim,nr,ll,a,vs,vr,aclk_iw)

    ! Only difference between solve and solvek is the output.
    ! in solvek, we do not sum over the modes.

    USE module_util

    implicit none
    include 'coupling.h'
   
    ! inputs:
    integer, intent(in) :: iw,nmode,ndim,nr
    integer, dimension(:), intent(in) :: ll
    real*8, intent(in) :: wtb
    complex*16, intent(in) :: w
    complex*16, dimension(:,:), intent(in) :: a
    complex*16, dimension(:,:), intent(in) :: vr
    complex*16, dimension(:), intent(in) :: wk
    complex*16, dimension(:), intent(in) :: vs

    ! outputs:
    complex*16, dimension(:,:), intent(out) :: aclk_iw

    ! internal:
    integer, parameter :: maxit = 10000
    real*8, parameter :: tol = 1.d-6
    integer :: i,md,im,ib,jb,tbd,info,imode
    integer :: it,icount,ir
    integer :: ntb,nindex,mindex
    real*8 :: err
    complex*16 :: cfac,alpha,beta,tmp1,tmp2
    complex*16 :: aclr
    complex*16 :: ii = (0.d0,1.0d0)
    logical :: ltmp
    integer, dimension(:), allocatable :: ipivot
    integer, dimension(:,:), allocatable :: tbl
    complex*16, dimension(:,:), allocatable :: atb
    complex*16, dimension(:,:), allocatable :: vsd
    complex*16, dimension(:), allocatable :: xs,xr,x0
    complex*16, dimension(:), allocatable ::  s0,sb0,p0,pb0,s1, &
                                             sb1,p1,pb1,vt

    ! allocate and initialize
    allocate(ipivot(ndim))
    allocate(tbl(nmode,4))
    allocate(atb(ndim,ndim))
    allocate(vsd(ndim,1))
    allocate(xs(ndim),xr(ndim),x0(ndim))
    allocate(s0(ndim),sb0(ndim),p0(ndim),pb0(ndim),&
             s1(ndim),sb1(ndim),p1(ndim),pb1(ndim),vt(ndim))

    ! build pre conditioning matrix
    ! determine which modes are in the target block:
    ntb    = 0
    nindex = 0
    mindex = 0
    do im = 1,nmode
      if(check_mode(w,wk(im),wtb)) then
        ntb = ntb + 1
        tbl(ntb,1) = im
        tbl(ntb,2) = nindex
        tbl(ntb,3) = mindex
        tbl(ntb,4) = 2*ll(im) + 1
        mindex = mindex + 2*ll(im) + 1
      endif
      nindex = nindex + 2*ll(im) + 1
    enddo

    ! deal with target block contribution
    if (ntb/=0) then
      ! dimension of the target system
      tbd = 0
      do ib=1,ntb
        tbd = tbd + tbl(ib,4)
      enddo

      ! build the coupling matrix for the
      ! target block
      do ib = 1,ntb
        do jb = 1,ntb
          atb(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4), &
              tbl(jb,3)+1:tbl(jb,3)+tbl(jb,4)) = &
            a(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4), &
              tbl(jb,2)+1:tbl(jb,2)+tbl(jb,4))
        enddo
      enddo

      ! perform an LU decomposition of the target block
      call zgetrf(tbd,tbd,atb,ndim,ipivot,info)

      ! assemble source vector
      do ib = 1,ntb
        vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
          vs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) ! /(ii*w) - HL-BM
      enddo
     
      ! solve the linear system for target block
      call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
      do ib = 1,ntb
        xs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
          vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
      end do
    endif

    ! deal with the rest of the matrix
    if (ntb /= nmode) then
      nindex = 0
      do im = 1,nmode
        md = 2*ll(im)+1
        if(.not.check_mode(w,wk(im),wtb)) then
          do i = 1,md
            cfac = 1.0d0/(wk(im)**2-w**2)
            xs(nindex+i) = cfac*vs(nindex+i)!/(ii*w) ! HL-BM
            
          end do
        end if
        nindex = nindex+md
      end do
    end if

    ! initial guess for the solution
    x0 = xs

    ! set the initial values for the vectors
    s0 = matmul(a,x0)
    if(ntb /= 0) then
      do ib = 1,ntb
        vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
          s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
      enddo
      call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
      do ib = 1,ntb
        s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
          vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
      enddo
    endif
    if(ntb /= nmode) then
      nindex = 0
      do im = 1,nmode
        md = 2*ll(im)+1
        ltmp = check_mode(w,wk(im),wtb)
        if(.not.ltmp) then
          do i = 1,md
            cfac = 1.0d0/(wk(im)**2-w**2)
            s0(nindex+i) = cfac*s0(nindex+i)
          enddo
        endif
        nindex = nindex+md
      enddo
    endif
    s0 = xs-s0
    sb0 = s0
    p0  = s0
    pb0 = sb0
    xr  = x0

    ! check the initial error
    tmp1 = my_dot(sb0,s0)
    err = sqrt(abs(tmp1/my_dot(xs,xs)))
    it  = 1 ! dummy initialization
    if(err > tol) then
 
      ! start the iteration loop
      icount = 0
      do it=1,maxit
 
        ! iterate the vectors
        vt = matmul(a,p0)
        if(ntb /= 0) then
          do ib = 1,ntb
            vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
              vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
          enddo
          call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
          do ib = 1,ntb
            vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
          enddo
        end if
      if(ntb /= nmode) then
           nindex = 0
           do im = 1,nmode
            md = 2*ll(im)+1
            if(.not.check_mode(w,wk(im),wtb)) then
              do i = 1,md
                cfac = 1.0d0/(wk(im)**2-w**2)
                vt(nindex+i) = cfac*vt(nindex+i)
              end do
            end if
            nindex = nindex+md
          enddo
        endif
        alpha = tmp1/my_dot(pb0,vt)
        s1 = s0-alpha*vt

        if(ntb /= 0) then
          do ib = 1,ntb
            vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
              pb0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
          enddo
          call zgetrs('T',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
          do ib = 1,ntb
            vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
          enddo
        endif
        if(ntb /= nmode) then
          nindex = 0
          do im = 1,nmode
            md = 2*ll(im)+1
            if(.not.check_mode(w,wk(im),wtb)) then
              do i = 1,md
                cfac = 1.0d0/(wk(im)**2-w**2)
                vt(nindex+i) = cfac*pb0(nindex+i)
              enddo
            endif
            nindex = nindex+md
          enddo
        endif
        vt = matmul(transpose(a),vt)

        sb1  = sb0-alpha*vt
        tmp2 = my_dot(sb1,s1)
        beta = tmp2/tmp1
        p1   = s1+beta*p0
        pb1  = sb1+beta*pb0

        ! update the solution
        xr = xr+alpha*p0

        ! estimate the error
        vt = matmul(a,xr)
        if(ntb /= 0) then
          do ib = 1,ntb
            vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
              vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
          enddo
          call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
          do ib = 1,ntb
            vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
          enddo
        endif
        if(ntb /= nmode) then
          nindex = 0
          do im = 1,nmode
            md = 2*ll(im)+1
            if(.not.check_mode(w,wk(im),wtb)) then
              do i = 1,md
                cfac = 1.0d0/(wk(im)**2-w**2)
                vt(nindex+i) = cfac*vt(nindex+i)
              enddo
            endif
            nindex = nindex+md
          enddo
        endif
        icount = it
        vt = vt-xs
 
        err = sqrt(abs(my_dot(vt,vt)/my_dot(xs,xs)))
        ! print *, 'iteration = ',it,' err = ',err/tol
        ! update the vectirs:
        s0 = s1
        sb0 = sb1
        p0 = p1
        pb0 = pb1
        tmp1 = tmp2

        if (err<tol) then
          !print*,'Total iteration = ',it,'for iteration',iw
          exit
        endif
      enddo
    else
       ! print *, 'no iterations',' err = ',err/tol
    end if

    if (it==maxit) then
      print*,'Max iterations at iteration',iw
    endif

    ! converged (unless max iterations hit)
    ! begin loop over the receivers but do not take the
    ! dot product.
    do ir=1,nr
      do imode=1,ndim
        aclr = vr(imode,ir)*xr(imode)
        aclk_iw(ir,imode) = aclr
      enddo
    enddo

    return

  END SUBROUTINE solvek_f

 SUBROUTINE solvek_a(w,wtb,iw,wk,nmode,ndim,nr,ll,a,vs,vr,aclk_iw)

   ! Only difference between solve and solvek is the output.
   ! in solvek, we do not sum over the modes.

   USE module_util

   implicit none
   include 'coupling.h'
  
   ! inputs:
   integer, intent(in) :: iw,nmode,ndim,nr
   integer, dimension(:), intent(in) :: ll
   real*8, intent(in) :: wtb
   complex*16, intent(in) :: w
   complex*16, dimension(:,:), intent(in) :: a
   complex*16, dimension(:,:), intent(in) :: vr
   complex*16, dimension(:), intent(in) :: wk
   complex*16, dimension(:), intent(in) :: vs

   ! outputs:
   complex*16, dimension(:), intent(out) :: aclk_iw

   ! internal:
   integer, parameter :: maxit = 10000
   real*8, parameter :: tol = 1.d-6
   integer :: i,md,im,ib,jb,tbd,info,imode
   integer :: it,icount,ir
   integer :: ntb,nindex,mindex
   real*8 :: err
   complex*16 :: cfac,alpha,beta,tmp1,tmp2
   complex*16 :: aclr
   complex*16 :: ii = (0.d0,1.0d0)
   logical :: ltmp
   integer, dimension(:), allocatable :: ipivot
   integer, dimension(:,:), allocatable :: tbl
   complex*16, dimension(:,:), allocatable :: atb
   complex*16, dimension(:,:), allocatable :: vsd
   complex*16, dimension(:), allocatable :: xs,xr,x0
   complex*16, dimension(:), allocatable ::  s0,sb0,p0,pb0,s1, &
                                            sb1,p1,pb1,vt

   ! allocate and initialize
   allocate(ipivot(ndim))
   allocate(tbl(nmode,4))
   allocate(atb(ndim,ndim))
   allocate(vsd(ndim,1))
   allocate(xs(ndim),xr(ndim),x0(ndim))
   allocate(s0(ndim),sb0(ndim),p0(ndim),pb0(ndim),&
            s1(ndim),sb1(ndim),p1(ndim),pb1(ndim),vt(ndim))

   ! build pre conditioning matrix
   ! determine which modes are in the target block:
   ntb    = 0
   nindex = 0
   mindex = 0
   do im = 1,nmode
     if(check_mode(w,wk(im),wtb)) then
       ntb = ntb + 1
       tbl(ntb,1) = im
       tbl(ntb,2) = nindex
       tbl(ntb,3) = mindex
       tbl(ntb,4) = 2*ll(im) + 1
       mindex = mindex + 2*ll(im) + 1
     endif
     nindex = nindex + 2*ll(im) + 1
   enddo

   ! deal with target block contribution
   if (ntb/=0) then
     ! dimension of the target system
     tbd = 0
     do ib=1,ntb
       tbd = tbd + tbl(ib,4)
     enddo

     ! build the coupling matrix for the
     ! target block
     do ib = 1,ntb
       do jb = 1,ntb
         atb(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4), &
             tbl(jb,3)+1:tbl(jb,3)+tbl(jb,4)) = &
           a(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4), &
             tbl(jb,2)+1:tbl(jb,2)+tbl(jb,4))
       enddo
     enddo

     ! perform an LU decomposition of the target block
     call zgetrf(tbd,tbd,atb,ndim,ipivot,info)

     ! assemble source vector
     do ib = 1,ntb
       vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
         vs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) ! /(ii*w) - HL-BM
     enddo
    
     ! solve the linear system for target block
     call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
     do ib = 1,ntb
       xs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
         vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
     end do
   endif

   ! deal with the rest of the matrix
   if (ntb /= nmode) then
     nindex = 0
     do im = 1,nmode
       md = 2*ll(im)+1
       if(.not.check_mode(w,wk(im),wtb)) then
          do i = 1,md
           cfac = 1.0d0/(wk(im)**2-w**2)
          xs(nindex+i) = cfac*vs(nindex+i)!/(ii*w) ! HL-BM
         end do
       end if
       nindex = nindex+md
     end do
   end if

   ! initial guess for the solution
   x0 = xs

   ! set the initial values for the vectors
   s0 = matmul(a,x0)
   if(ntb /= 0) then
     do ib = 1,ntb
       vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
         s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
     enddo
     call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
     do ib = 1,ntb
       s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
         vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
     enddo
   endif
   if(ntb /= nmode) then
     nindex = 0
     do im = 1,nmode
       md = 2*ll(im)+1
       ltmp = check_mode(w,wk(im),wtb)
       if(.not.ltmp) then
         do i = 1,md
           cfac = 1.0d0/(wk(im)**2-w**2)
           s0(nindex+i) = cfac*s0(nindex+i)
         enddo
       endif
       nindex = nindex+md
     enddo
   endif
   s0 = xs-s0
   sb0 = s0
   p0  = s0
   pb0 = sb0
   xr  = x0

   ! check the initial error
   tmp1 = my_dot(sb0,s0)
   err = sqrt(abs(tmp1/my_dot(xs,xs)))
   it  = 1 ! dummy initialization
   if(err > tol) then

     ! start the iteration loop
     icount = 0
     do it=1,maxit

       ! iterate the vectors
       vt = matmul(a,p0)
       if(ntb /= 0) then
         do ib = 1,ntb
           vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
             vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
         enddo
         call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
         do ib = 1,ntb
           vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
             vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
         enddo
       end if
     if(ntb /= nmode) then
          nindex = 0
          do im = 1,nmode
           md = 2*ll(im)+1
           if(.not.check_mode(w,wk(im),wtb)) then
             do i = 1,md
               cfac = 1.0d0/(wk(im)**2-w**2)
               vt(nindex+i) = cfac*vt(nindex+i)
             end do
           end if
           nindex = nindex+md
         enddo
       endif
       alpha = tmp1/my_dot(pb0,vt)
       s1 = s0-alpha*vt

       if(ntb /= 0) then
         do ib = 1,ntb
           vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
             pb0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
         enddo
         call zgetrs('T',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
         do ib = 1,ntb
           vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
             vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
         enddo
       endif
       if(ntb /= nmode) then
         nindex = 0
         do im = 1,nmode
           md = 2*ll(im)+1
           if(.not.check_mode(w,wk(im),wtb)) then
             do i = 1,md
               cfac = 1.0d0/(wk(im)**2-w**2)
               vt(nindex+i) = cfac*pb0(nindex+i)
             enddo
           endif
           nindex = nindex+md
         enddo
       endif
       vt = matmul(transpose(a),vt)

       sb1  = sb0-alpha*vt
       tmp2 = my_dot(sb1,s1)
       beta = tmp2/tmp1
       p1   = s1+beta*p0
       pb1  = sb1+beta*pb0

       ! update the solution
       xr = xr+alpha*p0

       ! estimate the error
       vt = matmul(a,xr)
       if(ntb /= 0) then
         do ib = 1,ntb
           vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
             vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
         enddo
         call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
         do ib = 1,ntb
           vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
             vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
         enddo
       endif
       if(ntb /= nmode) then
         nindex = 0
         do im = 1,nmode
           md = 2*ll(im)+1
           if(.not.check_mode(w,wk(im),wtb)) then
             do i = 1,md
               cfac = 1.0d0/(wk(im)**2-w**2)
               vt(nindex+i) = cfac*vt(nindex+i)
             enddo
           endif
           nindex = nindex+md
         enddo
       endif
       icount = it
       vt = vt-xs

       err = sqrt(abs(my_dot(vt,vt)/my_dot(xs,xs)))
       ! print *, 'iteration = ',it,' err = ',err/tol
       ! update the vectirs:
       s0 = s1
       sb0 = sb1
       p0 = p1
       pb0 = pb1
       tmp1 = tmp2

       if (err<tol) then
         !print*,'Total iteration = ',it,'for iteration',iw
         exit
       endif
     enddo
   else
      ! print *, 'no iterations',' err = ',err/tol
   end if

   if (it==maxit) then
     print*,'Max iterations at iteration',iw
   endif

   ! converged (unless max iterations hit)
   ! begin loop over the receivers but do not take the
   ! dot product.
   do imode=1,ndim
     aclr = vr(imode,1)*xr(imode)
     aclk_iw(imode) = aclr
   enddo

   return

 END SUBROUTINE solvek_a

  SUBROUTINE solvek_f2(w,wtb,iw,wk,nmode,ndim,nr,ll,a,vs,aclk_iw)

    ! is similar to solvek_f but doesn't do product with receiver function
    ! use w/ mat_op2

    USE module_util

    implicit none
    include 'coupling.h'
   
    ! inputs:
    integer, intent(in) :: iw,nmode,ndim,nr
    integer, dimension(:), intent(in) :: ll
    real*8, intent(in) :: wtb
    complex*16, intent(in) :: w
    complex*16, dimension(:,:), intent(in) :: a
    !complex*16, dimension(:,:), intent(in) :: vr
    complex*16, dimension(:), intent(in) :: wk
    complex*16, dimension(:), intent(in) :: vs

    ! outputs:
    ! complex*16, dimension(:,:), intent(out) :: aclk_iw
    complex*16, dimension(:), intent(out) :: aclk_iw

    ! internal:
    integer, parameter :: maxit = 10000
    real*8, parameter :: tol = 1.d-6
    integer :: i,md,im,ib,jb,tbd,info,imode
    integer :: it,icount,ir
    integer :: ntb,nindex,mindex
    real*8 :: err
    complex*16 :: cfac,alpha,beta,tmp1,tmp2
    complex*16 :: aclr
    complex*16 :: ii = (0.d0,1.0d0)
    logical :: ltmp
    integer, dimension(:), allocatable :: ipivot
    integer, dimension(:,:), allocatable :: tbl
    complex*16, dimension(:,:), allocatable :: atb
    complex*16, dimension(:,:), allocatable :: vsd
    complex*16, dimension(:), allocatable :: xs,xr,x0
    complex*16, dimension(:), allocatable ::  s0,sb0,p0,pb0,s1, &
                                             sb1,p1,pb1,vt

    ! allocate and initialize
    allocate(ipivot(ndim))
    allocate(tbl(nmode,4))
    allocate(atb(ndim,ndim))
    allocate(vsd(ndim,1))
    allocate(xs(ndim),xr(ndim),x0(ndim))
    allocate(s0(ndim),sb0(ndim),p0(ndim),pb0(ndim),&
             s1(ndim),sb1(ndim),p1(ndim),pb1(ndim),vt(ndim))

    ! build pre conditioning matrix
    ! determine which modes are in the target block:
    ntb    = 0
    nindex = 0
    mindex = 0
    do im = 1,nmode
      if(check_mode(w,wk(im),wtb)) then
        ntb = ntb + 1
        tbl(ntb,1) = im
        tbl(ntb,2) = nindex
        tbl(ntb,3) = mindex
        tbl(ntb,4) = 2*ll(im) + 1
        mindex = mindex + 2*ll(im) + 1
      endif
      nindex = nindex + 2*ll(im) + 1
    enddo

    ! deal with target block contribution
    if (ntb/=0) then
      ! dimension of the target system
      tbd = 0
      do ib=1,ntb
        tbd = tbd + tbl(ib,4)
      enddo

      ! build the coupling matrix for the
      ! target block
      do ib = 1,ntb
        do jb = 1,ntb
          atb(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4), &
              tbl(jb,3)+1:tbl(jb,3)+tbl(jb,4)) = &
            a(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4), &
              tbl(jb,2)+1:tbl(jb,2)+tbl(jb,4))
        enddo
      enddo

      ! perform an LU decomposition of the target block
      call zgetrf(tbd,tbd,atb,ndim,ipivot,info)

      ! assemble source vector
      do ib = 1,ntb
         vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
              vs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) ! /(ii*w) - HL-BM
      enddo
     
      ! solve the linear system for target block
      call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
      do ib = 1,ntb
        xs(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
          vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
      end do
    endif

    ! deal with the rest of the matrix
    if (ntb /= nmode) then
      nindex = 0
      do im = 1,nmode
        md = 2*ll(im)+1
        if(.not.check_mode(w,wk(im),wtb)) then
           do i = 1,md
            cfac = 1.0d0/(wk(im)**2-w**2)
            xs(nindex+i) = cfac*vs(nindex+i)! /(ii*w) - HL-BM
            
          end do
        end if
        nindex = nindex+md
      end do
    end if

    ! initial guess for the solution
    x0 = xs

    ! set the initial values for the vectors
    s0 = matmul(a,x0)
    if(ntb /= 0) then
      do ib = 1,ntb
        vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
          s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
      enddo
      call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
      do ib = 1,ntb
        s0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
          vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
      enddo
    endif
    if(ntb /= nmode) then
      nindex = 0
      do im = 1,nmode
        md = 2*ll(im)+1
        ltmp = check_mode(w,wk(im),wtb)
        if(.not.ltmp) then
          do i = 1,md
            cfac = 1.0d0/(wk(im)**2-w**2)
            s0(nindex+i) = cfac*s0(nindex+i)
          enddo
        endif
        nindex = nindex+md
      enddo
    endif
    s0 = xs-s0
    sb0 = s0
    p0  = s0
    pb0 = sb0
    xr  = x0

    ! check the initial error
    tmp1 = my_dot(sb0,s0)
    err = sqrt(abs(tmp1/my_dot(xs,xs)))
    it  = 1 ! dummy initialization
    if(err > tol) then
 
      ! start the iteration loop
      icount = 0
      do it=1,maxit
 
        ! iterate the vectors
        vt = matmul(a,p0)
        if(ntb /= 0) then
          do ib = 1,ntb
            vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
              vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
          enddo
          call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
          do ib = 1,ntb
            vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
          enddo
        end if
      if(ntb /= nmode) then
           nindex = 0
           do im = 1,nmode
            md = 2*ll(im)+1
            if(.not.check_mode(w,wk(im),wtb)) then
              do i = 1,md
                cfac = 1.0d0/(wk(im)**2-w**2)
                vt(nindex+i) = cfac*vt(nindex+i)
              end do
            end if
            nindex = nindex+md
          enddo
        endif
        alpha = tmp1/my_dot(pb0,vt)
        s1 = s0-alpha*vt

        if(ntb /= 0) then
          do ib = 1,ntb
            vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
              pb0(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
          enddo
          call zgetrs('T',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
          do ib = 1,ntb
            vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
          enddo
        endif
        if(ntb /= nmode) then
          nindex = 0
          do im = 1,nmode
            md = 2*ll(im)+1
            if(.not.check_mode(w,wk(im),wtb)) then
              do i = 1,md
                cfac = 1.0d0/(wk(im)**2-w**2)
                vt(nindex+i) = cfac*pb0(nindex+i)
              enddo
            endif
            nindex = nindex+md
          enddo
        endif
        vt = matmul(transpose(a),vt)

        sb1  = sb0-alpha*vt
        tmp2 = my_dot(sb1,s1)
        beta = tmp2/tmp1
        p1   = s1+beta*p0
        pb1  = sb1+beta*pb0

        ! update the solution
        xr = xr+alpha*p0

        ! estimate the error
        vt = matmul(a,xr)
        if(ntb /= 0) then
          do ib = 1,ntb
            vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1) = &
              vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4))
          enddo
          call zgetrs('N',tbd,1,atb,ndim,ipivot,vsd,ndim,info)
          do ib = 1,ntb
            vt(tbl(ib,2)+1:tbl(ib,2)+tbl(ib,4)) = &
              vsd(tbl(ib,3)+1:tbl(ib,3)+tbl(ib,4),1)
          enddo
        endif
        if(ntb /= nmode) then
          nindex = 0
          do im = 1,nmode
            md = 2*ll(im)+1
            if(.not.check_mode(w,wk(im),wtb)) then
              do i = 1,md
                cfac = 1.0d0/(wk(im)**2-w**2)
                vt(nindex+i) = cfac*vt(nindex+i)
              enddo
            endif
            nindex = nindex+md
          enddo
        endif
        icount = it
        vt = vt-xs
 
        err = sqrt(abs(my_dot(vt,vt)/my_dot(xs,xs)))
        ! print *, 'iteration = ',it,' err = ',err/tol
        ! update the vectirs:
        s0 = s1
        sb0 = sb1
        p0 = p1
        pb0 = pb1
        tmp1 = tmp2

        if (err<tol) then
          !print*,'Total iteration = ',it,'for iteration',iw
          exit
        endif
      enddo
    else
       ! print *, 'no iterations',' err = ',err/tol
    end if

    if (it==maxit) then
      print*,'Max iterations at iteration',iw
    endif

    ! ! converged (unless max iterations hit)
    ! ! begin loop over the receivers but do not take the
    ! ! dot product.
    ! do ir=1,nr
    !   do imode=1,ndim
    !     aclr = vr(imode,ir)*xr(imode)
    !     aclk_iw(ir,imode) = aclr
    !   enddo
    ! enddo

      do imode = 1,ndim ! put back solution in vector
        aclk_iw(imode) = xr(imode)
      end do


    return

  END SUBROUTINE solvek_f2

  function check_mode(w,wm,wtb)

    implicit none
    logical(kind(.true.)) :: check_mode
    complex*16 :: w,wm
    real*8:: wtb

!    if(abs(real(/w(im)-w)) <= wtb) then
    if(abs(real(wm-w)) <= wtb) then
       check_mode = .true.
    else
       check_mode = .false.
    end if

    return
  end function check_mode

   



END MODULE disol_module

