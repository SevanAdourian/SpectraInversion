MODULE receiver_source_module

  ! ================================================================ !
  ! Build receiver and source vectors                                !
  ! ================================================================ !

  USE futil_module, ONLY : get_mode_fun,&
                           lgndr

  IMPLICIT NONE

  PRIVATE 
  PUBLIC :: build_r,&
            build_s

CONTAINS

  SUBROUTINE build_r(r,nmodes,nn,types,ll,rlat,rlon,nu)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! build receiver vector                             !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! receiver terms are in complex spherical harmonics !
    ! (note, subroutine receiver only provides 0<=m<=l  !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    IMPLICIT NONE
    include "coupling.h"
    ! input/output:
    integer, intent(in) :: nmodes,nn(:),types(:),ll(:)
    real*8, intent(in)  :: rlat,rlon,nu(3)
    complex*16, intent(out) :: r(:) 
    ! internal:
    integer :: mode,n,type,l,istart,i,m,ierr
    complex*16, dimension(:), allocatable :: rec

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    istart=1
    do mode=1,nmodes
      n=nn(mode)
      type=types(mode)
      l=ll(mode)
      ALLOCATE(rec(l+1),STAT=ierr)
      IF (ierr/=0) STOP 'allocate error 1 in SUB build_r, STOP!'

      call receiver(n,type,l,rlat,rlon,nu,rec)
      do m=-l,l
        i=istart+l+m
          if(m.lt.0) then
            r(i)=((-1.0)**m)*conjg(rec(-m+1))
          else
            r(i)=rec(m+1)
          endif
        enddo
        istart=istart+2*l+1
        DEALLOCATE(rec)
      enddo

      return
  END SUBROUTINE build_r

  ! ================================================================ !

  SUBROUTINE receiver(n,type,l0,lat,lon,nu,r)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! calculate r = nu dot s^*                          !
    !                                                   !
    ! Eq. (D7) invokes D4-D6 and D1 in DT1998           !
    ! Position of receiver is at top of solid Earth     !
    ! i.e., rr = r_surface (w/o ocean) and r_ocean- in  !
    ! model with ocean                                  !   
    !  inputs:                                          !
    !  - nu   = unit polarization vector (vert,n-s,e-w) !
    !  - s    = PREM eigenvector at receiver            !
    !  - n    = overtone number                         !
    !  - type = S is 1; T is 0                          !
    !  - lat  = receiver lat (deg)                      !
    !  - lon  = receiver lon (deg)                      !
    !                                                   !
    ! Note - this is dimensional                        !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    USE modellayer_module, ONLY : NR
    USE eigfun_module, ONLY : eiginfo, eigac
    IMPLICIT NONE
    include "coupling.h"
    ! input/output:
    integer, intent(in) :: n,type,l0
    real*8, intent(in) :: lat,lon,nu(3)
    complex*16, intent(out) :: r(:)
    ! internal:
    character(len=1) :: ctype
    integer :: m,i,l
    real*8 :: om,q,gv,av,ahs,aht,ap
    real*8 :: theta,phi,sint,cost,cosect
    real*8 :: x(LMAX+1),dx(LMAX+1)
    complex*16 :: expphi,expmphi
    real*8, dimension(:), pointer :: radius,&
                          u,du,v,dv,w,dw,p,dp
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    ! note: av (U*), ahs (V*/k), aht (W/k), ap(P)
    ! where gravity effect is taken into account U -> U*
    ! see DT1998 chap. 10.4

    r=CMPLX(0.0,0.0)
    l=l0
    if(type .eq. 0) then
       ctype = 't'
    else
       ctype = 's'
    endif

    
    NULLIFY(radius,u,du,v,dv,w,dw,p,dp)
    call get_mode_fun(n,ctype,l,om,q,av,ahs,aht,ap,&
                      radius,u,du,v,dv,w,dw,p,dp)
    IF (SIZE(radius)/=NR .OR. SIZE(u)/=NR .OR. SIZE(w)/=NR) STOP &
       ' error in dimension of eigenfunctions in SUB receiver, STOP!' 

    theta=PI/2.0-lat*PI/180.0
    phi=lon*PI/180.0
    sint=sin(theta)
    cost=cos(theta)
    cosect=1.0/sint
    call lgndr(l,cost,sint,x,dx) 
    IF (l/=l0) STOP 'please check SUB lgndr, STOP'
    expphi=cexp(cmplx(0.0,phi))
    expmphi=cmplx(1.0,0.0)

    ! normalizations to dimensionalize:
    av  = (av/RA) * om * a_norm
    ahs = (ahs/RA) * om * a_norm
    ! added aht normalization
    aht = (aht/RA) * om * a_norm !Nov9

    !!bm!! BM23 - make receiver vector very simple:
    
!!bm!!    do m=0,l
!!bm!!      i=m+1
!!bm!!      r(i)= ( nu(1)*av*x(i) + &
!!bm!!              nu(2)*(ahs*dx(i)+cmplx(0.0,REAL(m))*aht*cosect*x(i)) + &
!!bm!!              nu(3)*(cmplx(0.0,REAL(m))*ahs*cosect*x(i)-aht*dx(i)) &
!!bm!!            ) * expmphi
!!bm!!      expmphi=expmphi*expphi
!!bm!!      r(i)=conjg(r(i))
!!bm!!    enddo
    
    do m=0,l
       i=m+1
       r(i) = nu(1)*u(NR)*x(i) * expmphi
       expmphi = expmphi*expphi
       !r(i) = conjg(r(i)) !HL-BM
    enddo
    
    DEALLOCATE(radius,u,du,v,dv,w,dw,p,dp)

    return
  
  END SUBROUTINE receiver

  ! ================================================================ !

  SUBROUTINE build_s(s,nmodes,nn,types,ll,&
                     elat,elon,depth,moment_tensor)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! source term, s_lm, in terms of complex spherical  !
    ! harmonics (-l<=m<=l) (13.237)                     !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    IMPLICIT NONE
    include "coupling.h"
    ! input/output:
    integer, intent(in) :: nmodes,nn(:),types(:),ll(:)
    real*8, intent(in) :: elat,elon,depth,moment_tensor(6)
    complex*16, intent(out) :: s(:)
    ! internal:
    integer :: mode,n,type,l,istart,i,m, ierr
    complex*16, dimension(:), allocatable :: sou 
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
       
    istart=1
    do mode=1,nmodes
       n=nn(mode)
       type=types(mode)
       l=ll(mode)
       ALLOCATE(sou(l+1),STAT=ierr)
       IF(ierr/=0) STOP 'allocate error 1 in SUB build_s, STOP'
       call source(n,type,l,elat,elon,depth,moment_tensor,sou)
       do m=-l,l
          i=istart+l+m
          if(m.lt.0) then
             s(i)=((-1.0)**m)*conjg(sou(-m+1))
          else
             s(i)=sou(m+1)
          endif
       enddo
       istart=istart+2*l+1
       DEALLOCATE(sou)
    enddo

    return

  END SUBROUTINE build_s

  ! ================================================================ !

  SUBROUTINE source(n,type,l0,lat,lon,depth,moment_tensor,s)

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! calculate s = M : E^*                             !
    !  - M = moment tensor                              !
    !        M(1) -> M_rr                               !
    !        M(2) -> M_tt                               !
    !        M(3) -> M_pp                               !
    !        M(4) -> M_rt                               !
    !        M(5) -> M_rp                               !
    !        M(6) -> M_tp                               !
    !  - E^* = complex conjugate of strain tensor       !
    !  - n = overtone number                            !
    !  - l0 = spherical harmonic degree                 !
    !  - lat/lon/depth = of event (deg/km)              !
    !                                                   !
    ! s is dimensional                                  !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    USE modellayer_module, ONLY : NR, NOC
    USE eigfun_module
      
    IMPLICIT NONE
    include "coupling.h"
    ! input/output:
    integer, intent(in) :: n,type,l0
    real*8, intent(in) :: lat,lon,depth,moment_tensor(6)
    complex*16, intent(out) :: s(:)
    ! internal:
    integer :: m,i,j,l
    real*8 :: om,q,gv,av,ahs,aht,ap
    real*8, dimension(:), pointer :: r,u,du,v,dv
    real*8, dimension(:), pointer :: w,dw,p,dp
    real*8 :: bl,rs
    real*8 :: us,dus,vs,dvs,ws,dws,xs,zs
    real*8 :: theta,phi,sint,cost,cosect,cott
    real*8, dimension(:), allocatable :: x,dx
    complex*16 :: expphi,expmmphi,e(6)
    character*1 :: ctype

    s = cmplx(0.0,0.0)
    l=l0
    
    ALLOCATE(x(l+1),dx(l+1),STAT=i)
    IF (i/=0) STOP 'allocate error 1 in SUB source, STOP'
    IF (SIZE(s)<l+1) STOP &
      'error 2: enlarge dimension of s in SUB source, STOP'
    NULLIFY(r,u,du,v,dv,w,dw,p,dp)

    if(type .eq. 0) then
      ctype = 't'
    else
      ctype = 's'
    endif

    call get_mode_fun(n,ctype,l,om,q,av,ahs,aht,ap,r,&
                      u,du,v,dv,w,dw,p,dp)
    IF (SIZE(r)/=NR .OR. SIZE(u)/=NR .OR. SIZE(w)/=NR) &
      STOP 'Dimension error in SUB source, STOP'

    bl=l*(l+1.0)
    rs=1.0-depth*1000.0/RA ! non-dimensionalized source radius
    i=NOC

    ! find model r just below source:
    do while(r(i).gt.rs) 
      i=i-1
    enddo
    ! rs coincides with a model radius:
    if(r(i).eq.rs) then 
      us=u(i)
      dus=du(i)
      vs=v(i)
      dvs=dv(i)
      ws=w(i)
      dws=dw(i)
      xs=dvs+(us-vs)/rs
      zs=dws-ws/rs
      write(*,*)'coincide model radius'
    else ! interpolate to find us, vs, etc. at the source radius
      call interpolate(r(i),u(i),du(i),r(i+1),&
                       u(i+1),du(i+1),rs,us,dus)
      call interpolate(r(i),v(i),dv(i),r(i+1),&
                       v(i+1),dv(i+1),rs,vs,dvs)
      call interpolate(r(i),w(i),dw(i),r(i+1),&
                       w(i+1),dw(i+1),rs,ws,dws)
      xs=dvs+(us-vs)/rs ! Auxiliary variable D20
      zs=dws-ws/rs
    endif

    theta=PI/2.0-lat*PI/180.0
    phi=lon*PI/180.0
    sint=sin(theta)
    cost=cos(theta)
    cosect=1.0/sint
    cott=cost/sint
    call lgndr(l,cost,sint,x,dx) 
    IF (l/=l0) STOP 'please check SUB lgndr, STOP'
    expphi=cexp(cmplx(0.0))
    expmmphi=cmplx(1.0,0.0)
    ! compute strain tensor
    do m=0,l  
      i=m+1
      !(D.14)*
      e(1)=cmplx(dus*x(i),0.0)   
      !(D.15)*     
      e(2)=cmplx((us*x(i)-vs*(cott*dx(i)-((m*cosect)**2.0-bl)*x(i)))/rs,0.0) &
           +cmplx(0.0,-m*ws*cosect*(dx(i)-cott*x(i))/rs)  
      !(D.16)* 
      e(3)=cmplx((us*x(i)+vs*(cott*dx(i)-((m*cosect)**2.0)*x(i)))/rs,0.0) &
           +cmplx(0.0,m*ws*cosect*(dx(i)-cott*x(i))/rs) 
      !2x(D.17)*
      e(4)=cmplx(xs*dx(i),0.0)+cmplx(0.0,-m*zs*cosect*x(i)) 
      !2x(D.18)*
      e(5)=cmplx(0.0,-m*xs*cosect*x(i))-cmplx(zs*dx(i),0.0) 
      !2(D.19)*
      e(6)=cmplx(0.0,-2.0*m*vs*cosect*(dx(i)-cott*x(i))/rs) &
           +cmplx(ws*(2.0*cott*dx(i)-(2.0*((m*cosect)**2.0)-bl)*x(i))/rs,0.0) 
      s(i)=cmplx(0.0,0.0)
      do j=1,6
        s(i)=s(i)+cmplx(moment_tensor(j),0.0)*e(j)  !(D.21)
      enddo
      s(i)=s(i)*expmmphi
      expmmphi=expmmphi*expphi
    enddo

    DEALLOCATE(x,dx)
    DEALLOCATE(r,u,du,v,dv,w,dw,p,dp)

    return

  END SUBROUTINE source

  ! ================================================================ !

  SUBROUTINE interpolate(r1,f1,df1,r2,f2,df2,rs,fs,dfs)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !
    ! interpolate finds the value of a function fs and  !
    ! its derivative dfs at the radius rs by            !
    ! interpolating given values of the function and    !
    ! its derivative just below the source at radius r1 !
    ! given by f1 and df1, and values just above the    !
    ! ar radius r2, given by f2 and df2.                !
    ! Note: r2 > rs > r1.                               !
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    IMPLICIT NONE
    ! input/output:
    real*8, intent(in) :: r1,f1,df1,r2,f2,df2,rs
    real*8, intent(out) :: fs,dfs
    ! internal:
    real*8 :: h,dr    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !

    h=rs-r1
    dr=r2-r1
    if((h.lt.0.0).or.(dr.lt.0.0)) then
      write(6,"('wrong input in interpolate')")
      call exit(1)
    endif
    ! 2nd order Taylor expansion:
    fs=f1+h*df1+0.5*h*h*(df2-df1)/dr
    ! 1st order Taylor expansion:  
    dfs=df1+h*(df2-df1)/dr 
    
    return

  END SUBROUTINE interpolate

  ! ================================================================ !

END MODULE receiver_source_module
