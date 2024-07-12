 MODULE futil_module
!
! Last modified by H.Y. Yang, Mar 2014
!
  IMPLICIT NONE

  INTERFACE intgrl
    MODULE PROCEDURE intgrl_sp
    MODULE PROCEDURE intgrl_dp
  END INTERFACE intgrl

  INTERFACE deriv
    MODULE PROCEDURE deriv_sp
    MODULE PROCEDURE deriv_dp
  END INTERFACE deriv

  INTERFACE checknaninf
    MODULE PROCEDURE checknaninf_r4
    MODULE PROCEDURE checknaninf_r4_2d
    MODULE PROCEDURE checknaninf_cpx4
    MODULE PROCEDURE checknaninf_cpx4_2d
  END INTERFACE checknaninf

  INTERFACE RconvertC
    MODULE PROCEDURE RconvertC_1D
    MODULE PROCEDURE RconvertC_2D_savemeo
!    MODULE PROCEDURE RconvertC_2D
  END INTERFACE RconvertC

  ! type for argument of sub RconvertC
  TYPE istruczero_type
    LOGICAL :: real
    LOGICAL :: imag
  END TYPE istruczero_type


  PRIVATE
  PUBLIC :: get_mode_fun, &
            intgrl, &
            deriv, &
            thrj, &
            fact, &
            lgndr, &
            splbasis, &
            choles, &
            checknaninf, &
            RconvertC, &
            create_Usub, &
            istruczero_type, &
            check_BiOrth_eigvec, &
            inquire_matrix_type, &
            norm_SEP_eigvec, &
            check_transpose, &
            dirac_delta

CONTAINS

  SUBROUTINE get_mode_fun(nord,type,l,omega,q,aver,ahors,ahort,ap,r,u,du,v,dv,w,dw,p,dp)
!     
!     read eigenfunctions of mode (nord,l) and corresponding
!     eigenfrequency(omega), quality factor (q) and other information.
!     
!     H.Y. Yang, June 2013
!
!     Note that 
!     (1) format of the binary file of eigenfunctions is "HY's format"
!   
!     (2) these eigenfunctions follow the normalization rules
!          defined in Dahlen & Tromp (1998).
!       output eigfunctions are modified to be consistency in Appendix D:
!        u = U
!        du = dU/dr
!      * v = V/k   
!      * dv = (dV/dr)/k
!      * w = W/k
!      * dw = (dW/dr)/k
!        p = P
!        dp = dP/dr
!     (3) aver, ahors ahort are calculated to consider the effects of gravity on
!    seismograms (see Chap 10.4)
!       aver  = U+Ufree+Upot "at the top of the solid earth" (,rather than
!                             "surface" becuase seismometer cannon't be placed
!                             in the water)
!       ahors = (V+Vtil+Vpot)/k at sthe top of the solid earth 
!             *Dividing k is to make consistent with output of V/k
!       ahort = W/k at the top of the solid earth
!             *Dividing k is to make consistent with output of W/k
!       ap = P at the top of the solid earth
!       Note that if the 1D earth model has sea and receiver is set to be at
!       surface, then ahors ~ 0
!
!     (4) omega in rad/s
!     (5) r is non-dimensional
!   

       USE modellayer_module, ONLY : NR, NOC
       USE prem_module, ONLY : g

       IMPLICIT NONE

       INCLUDE "coupling.h"  ! use GRAV_SURFACE, NRMAX, NHEAD, funTor,
                             ! funSph, funRad, RA,GRAV,RHOAV

       INTEGER,           INTENT(IN   ) :: nord, l
       CHARACTER(LEN=*),  INTENT(IN   ) :: type
       REAL*8,              INTENT(  OUT) :: omega, q
       REAL*8, INTENT(OUT) ::  aver, ahors, ahort,ap
       REAL*8, POINTER :: r(:), u(:), du(:), v(:), dv(:), w(:), dw(:), p(:), dp(:)
       REAL*8  :: gsurf
       INTEGER :: jcom, llmin, llmax, nmin, nmax
       REAL :: wmin1, wmax1, wgrav
       INTEGER :: n, nicfun, nocfun, ifanis
       REAL :: trefl
       REAL :: model_store(NRMAX,9)
       INTEGER :: lbottom, nrfun, mrec, nmodes
       CHARACTER(LEN=8) :: kfmt
       REAL*8, ALLOCATABLE :: eig(:)

       INTEGER :: ntbl, ltbl
       REAL*8    :: gv, err
       REAL*8  :: oms, k, av, ah

       CHARACTER(LEN=128) :: bin_file     
       INTEGER :: ierr, i, j
       LOGICAL :: ishit
       ishit=.FALSE.


       SELECT CASE (type) 
       CASE ('T','t')
          bin_file = TRIM(funTor)
       CASE ('S','s')
          IF (l == 0) THEN
            bin_file = TRIM(funRad)            
          ELSE 
            bin_file = TRIM(funSph)
          END IF
       CASE DEFAULT
         STOP 'mode type is out of range in SUB get_mode_fun, STOP!!'
       END SELECT 

       OPEN(UNIT=1,FILE=bin_file,STATUS='OLD',FORM='UNFORMATTED', &
            ACCESS='SEQUENTIAL',IOSTAT=ierr)
       IF (ierr>0) THEN
         WRITE(*,*) TRIM(bin_file), ' reading error, STOP'
         STOP
       END IF
       READ(UNIT=1) kfmt
       IF ( kfmt == 'VFUN1.02' ) THEN
         READ(UNIT=1) jcom,wmin1,wmax1,llmin,llmax,wgrav
       ELSE IF ( kfmt == 'VFUN1.04') THEN
         READ(UNIT=1) jcom,wmin1,wmax1,nmodes,nmin,nmax,llmin,llmax,wgrav
       ELSE
         STOP 'reading eigfunction error, STOP' 
       END IF

       READ(UNIT=1) n, nicfun, nocfun, ifanis, trefl, ((model_store(i,j),i=1,n),j=1,9)
       model_store(:,1)=model_store(:,1)/RA

       IF (n/=NR) STOP 'error1 in SUB get_model_fun, STOP'
       IF (.NOT.ASSOCIATED(g)) STOP &
          ' g has to be defined before get_model_fun, STOP'

       IF (ASSOCIATED(r)) NULLIFY(r)
       IF (ASSOCIATED(u)) NULLIFY(u)
       IF (ASSOCIATED(du)) NULLIFY(du)
       IF (ASSOCIATED(v)) NULLIFY(v)
       IF (ASSOCIATED(dv)) NULLIFY(dv)
       IF (ASSOCIATED(w)) NULLIFY(w)
       IF (ASSOCIATED(dw)) NULLIFY(dw)
       IF (ASSOCIATED(p)) NULLIFY(p)
       IF (ASSOCIATED(dp)) NULLIFY(dp)

       ALLOCATE(r(n),u(n),du(n),v(n),dv(n),w(n),dw(n),p(n),dp(n),STAT=ierr)
       IF (ierr/=0 .OR. n==0) STOP &
           'allocate error in sub get_mode_fun, STOP!'

       IF (jcom==3 .OR. jcom==1) THEN
         nrfun=n
         mrec=6*nrfun+nhead
         lbottom=0
       ELSE IF (jcom==2) THEN
         nrfun=n-nocfun 
         mrec=2*nrfun+nhead
         lbottom=nocfun
       END IF
       ALLOCATE(eig(mrec),STAT=ierr)
       IF (ierr/=0) STOP &
          'allocate error of eig in SUB get_mode_fun, STOP'

       DO
         READ(UNIT=1,IOSTAT=ierr) (eig(j),j=1,mrec)
         IF (ierr<0) THEN
           EXIT
         ELSE IF (ierr >0) THEN
           STOP 'reading funfile error in sub get_mode_fun, STOP'
         END IF

         ntbl=NINT(eig(1))
         ltbl=NINT(eig(2))
         !=========================================!
         omega=eig(3)
         omega = omega/f_norm ! this is different from MYTCC
         !=========================================!
         q=eig(4)
         gv=eig(5)
         err=eig(6)

         IF ((ntbl == nord).AND.(ltbl == l)) THEN
            ishit = .TRUE.
            k = SQRT((eig(2)+1)*eig(2))
            oms =  eig(3)/f_norm
            eig(1:mrec-nhead)=eig(nhead+1:mrec)*oms  
            u=0.0
            du=0.0
            v=0.0
            dv=0.0
            p=0.0
            dp=0.0
            w=0.0
            dw=0.0
            ahors=0.0
            ahort=0.0
            aver=0.0
            ap=0.0
            ! spheroidal and radial modes
            IF (jcom == 1 .OR. jcom==3) THEN
              gsurf=g(NOC)
              av=eig(NOC)+(2.d0*gsurf*eig(NOC)+&
                          (ltbl+1.d0)*eig(4*nrfun+NOC)) &
                          /(oms**2)/DBLE(model_store(NOC,1))
              aver=av
              ap=eig(4*nrfun+NOC)

              IF (jcom==3) THEN
                ah=eig(2*nrfun+NOC)-(gsurf*eig(NOC)+&
                                     eig(4*nrfun+NOC))*k &
                                    /(oms**2)/DBLE(model_store(NOC,1))
                ah=ah/k       !make consistent with the definition below
                ahors=ah

                DO i =1,nrfun
                  r(i) =model_store( i,1)  
                  u(i) =eig(    i)  
                  du(i)=eig(  nrfun+i)  
                  v(i) =eig(2*nrfun+i)/k
                  dv(i)=eig(3*nrfun+i)/k
                  p(i) =eig(4*nrfun+i)  
                  dp(i)=eig(5*nrfun+i)   
                END DO

              ELSE IF ( jcom == 1) THEN
                DO i =1,nrfun
                  r(i) =model_store( i,1)  
                  u(i) =eig(    i)  
                  du(i)=eig(  nrfun+i)  
                  p(i) =eig(4*nrfun+i)  
                  dp(i)=eig(5*nrfun+i)  
                END DO
              END IF
             
                                  ! toroidal modes
            ELSE IF (jcom == 2) THEN
              ah=eig(NOC-lbottom)
              ah=ah/k   !make consistent with the definition below
              ahort=ah

              DO i=1,nrfun
                r(lbottom+i) =model_store(lbottom+i,1)
                w(lbottom+i) =eig(  i)/k
                dw(lbottom+i)=eig(nrfun+i)/k
              END DO
              DO i=1,lbottom
                r(i)=model_store(i,1)
              END DO

            ELSE
              STOP 'jcom in binary eigenfunction file is out of range, STOP' 
            END IF

            EXIT
         END IF ! end if matched mode 
       END DO
 
       IF (.NOT.ishit) THEN
         WRITE(*,'(A,I5,A,I5,A)')'Cannot find mode ',nord,type,l, ' in binary file, STOP!!'
         STOP
       ENDIF
       DEALLOCATE(eig)
       CLOSE(UNIT=1)

  END SUBROUTINE get_mode_fun

  SUBROUTINE intgrl_sp(sum1,r,nir,ner,f,s1,s2,s3,sum1r)
!
!       Computes the integral of f[i]*r[i]*r[i] from i=nir to i=ner for
!       radii values as in model stored in f
!       sum1: result at r=r(ner)
!       sum1r: result as a function of r(nir:ner), but the dimension
!              of the array is the same as r (i.e. NR), rather than
!              ner-nir+1 
!
        USE modellayer_module, ONLY:NR,kdis
        IMPLICIT NONE
        REAL,   INTENT(IN    ) :: r(:),f(:)
        INTEGER,INTENT(IN    ) :: nir,ner
        REAL,   INTENT(   OUT) :: sum1
        REAL,   INTENT(   OUT) :: s1(:),s2(:),s3(:)
        REAL,   INTENT(   OUT), OPTIONAL :: sum1r(:)
  
        INTEGER :: ndis
                  ! DIMENSION(NR)
        REAL, DIMENSION(:),ALLOCATABLE :: yprime, sum0  
        REAL :: third,fifth,sixth
        REAL :: rji,dsum 
        INTEGER :: n 
 
        INTEGER :: ierr,i,j,nir1

       third = 1.0/3.0
       fifth = 1.0/5.0
       sixth = 1.0/6.0

       n=SIZE(r)
       IF (n/=SIZE(f).OR.n/=NR) STOP 'dimension inconsistent in SUB intgrl, STOP'
       ALLOCATE(yprime(n),sum0(n),STAT=ierr)
       IF (ierr/=0.OR.n==0) STOP 'allocate error 1 in SUB intgrl, STOP'

                        ! use discontinuities r-
       IF (.NOT. ASSOCIATED(kdis)) STOP 'kdis has to be defined previously, STOP' 
       ndis=SIZE(kdis)
   
                        ! calculate cubic splines 
	call deriv(f,yprime,n,r,ndis,kdis,s1,s2,s3)

        sum0=0.0
	nir1 = nir + 1
	DO i=nir1,ner
	  j = i-1
	  rji = r(i) - r(j)
               ! Note that "r(j)**2" at the beginning of RHS represent 
               !  the integral is multiplied by r^2
          dsum=r(j)*r(j)*rji*(f(j)+rji*(.50*s1(j)+rji*(third*s2(j)+rji* &
                  .250*s3(j))))+2.0*r(j)*rji*rji*(.50*f(j)+rji*(third*s1(j)+rji* &
                 (.250*s2(j)+rji*fifth*s3(j))))+rji*rji*rji*(third*f(j)+rji* &
                 (.250*s1(j)+rji*(fifth*s2(j)+rji*sixth*s3(j))))
          sum0(i)=sum0(j)+dsum 
        END DO

        sum1=sum0(ner)
        IF (PRESENT(sum1r)) THEN
          sum1r=sum0
          IF (nir/=1) THEN
             WRITE(*,*) 'Wanring!! beginning entry in r for SUB intgrl is not "1" '
          END IF
        END IF
       
        DEALLOCATE(yprime,sum0)

	return
  END SUBROUTINE intgrl_sp 

  SUBROUTINE deriv_sp(yr4,yprime,n,rr4,ndis,kdis,s1r4,s2r4,s3r4)
!     calculate the coefficients (s1,s2,s3) of cubic spline
!     as a function of r
!
!     yprime=s1
!    
!     Numerical error occurrs in the portion of y=constant,
!     in which s1,s2,s3 have to be zero but the numerical results are not.
!     This problem can be resolved by changing data type of variables
!     from REAL*4 to REAL*8.
!     Even though inaccurate s1,s2,s3 do not affect integration (SUB intgrl).
!
! 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: ndis, kdis(:)      !kdis(ndis)
      REAL,    INTENT(IN) :: yr4(:), rr4(:)     !y(n),r(n)
      REAL,    INTENT(OUT):: yprime(n),s1r4(n),s2r4(n),s3r4(n)

      REAL*8, DIMENSION(n) :: y,r,s1,s2,s3
      REAL*8 :: yy(3), f(3,n)

      INTEGER :: j1,j2,nd,ndp,i,j,k
      REAL*8    :: h,h2,a0,b0,b1,ha,h2a,h3a,hb,h2b,s13,s21,s32,y0
      equivalence (yy(1),y0)
      yy=0.0
      r=DBLE(rr4)
      y=DBLE(yr4)

      ndp=ndis+1
      do 3 nd=1,ndp
      if(nd.eq.1) go to 4
      if(nd.eq.ndp) go to 5
      j1=kdis(nd-1)+1
      j2=kdis(nd)-2
      go to 6
    4 j1=1
      j2=kdis(1)-2
      go to 6
    5 j1=kdis(ndis)+1
      j2=n-2
    6 if((j2+1-j1).gt.0) go to 11
      j2=j2+2
      y0=(y(j2)-y(j1))/(r(j2)-r(j1))
	s1(j1)=yy(1)
	s1(j2)=yy(1)
	s2(j1)=yy(2)
	s2(j2)=yy(2)
	s3(j1)=yy(3)
	s3(j2)=yy(3)
      go to 3
   11 a0=0.0
      if(j1.eq.1) go to 7
      h=r(j1+1)-r(j1)
      h2=r(j1+2)-r(j1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
      b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
      go to 8
    7 b0=0.0
    8 b1=b0
      if(j2 .gt. n)write(0,'("error:deriv:j2= ",i5)')j2
      do 1 i=j1,j2
      h=r(i+1)-r(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.0*a0
      h3a=2.0*h-3.0*a0
      h2b=h2*b0
      s1(i)=h2/ha
      s2(i)=-ha/(h2a*h2)
      s3(i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.0*y0*ha)/(h*h3a)
      a0=s3(i)
    1 b0=f(3,i)
      i=j2+1
      h=r(i+1)-r(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.*h-a0)
      s1(i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=r(j2)-r(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      s3(i)=(y0*h2a+h2b)/(h*h2*(h-2.0*a0))
!     the following statements were expanded to prevent register overflow on unix
      s13=s1(i)*s3(i)
      s2(i)=f(1,i)-s13
      do 2 j=j1,j2
      k=i-1
      s32=s3(k)*s2(i)
      s1(i)=f(3,k)-s32
      s21=s2(k)*s1(i)
      s3(k)=f(2,k)-s21
      s13=s1(k)*s3(k)
      s2(k)=f(1,k)-s13
    2 i=k
      s1(i)=b1
      j2=j2+2
	s1(j2)=yy(1)
	s2(j2)=yy(2)
	s3(j2)=yy(3)
    3 continue

      s1r4=SNGL(s1)
      s2r4=SNGL(s2)
      s3r4=SNGL(s3)
      do 20 i=1,n
   20 yprime(i)=s1r4(i)
      return
  END SUBROUTINE deriv_sp

  FUNCTION thrj(j1,j2,j3,m1,m2,m3)
!
!     evaluate of Wigner 3-j coefficients,
!       but always outputs thrj in the form of single precision
!
!     single/double/quad precision
!   
!    tested by H.Y., Jan 2014
!    j<10: single precision
!    j<20: double precision. otherwise quad precision
!
     IMPLICIT NONE
     INTEGER, PARAMETER ::  r8b = SELECTED_REAL_KIND( p=33, r=4931 )  !quad precision
!                           r8b = SELECTED_REAL_KIND( p=6,  r=37  )  !single
!                           r8b = SELECTED_REAL_KIND( p=15, r=307 )  !double


     REAL :: thrj
     INTEGER, INTENT(IN) :: j1,j2,j3,m1,m2,m3

     REAL(KIND=r8b), ALLOCATABLE :: y(:)  !DIMENSION(m+1)
     INTEGER :: ja,jb,jc,ma,mb,mc,lm1,lm2,m,mx,my,n
     INTEGER :: ierr
     REAL(KIND=r8b) :: ss,alpha,beta,gamma
 
      thrj=0.0
      
      if(j1+j2-j3.lt.0.or.j2+j3-j1.lt.0.or.j3+j1-j2.lt.0) return
      if(j1-iabs(m1).lt.0) return
      if(j2-iabs(m2).lt.0) return
      if(j3-iabs(m3).lt.0) return
      if(m1+m2+m3.ne.0) return
       ! to avoid numerical error, add this criteria (H.Y. 2013)
       ! eqn (C.219) in D&T (1998)
      IF((m1==0).AND.(m2==0).AND.(m3==0).AND.MOD(j1+j2+j3,2)/=0) RETURN
      
!-----
!     use symmetries to make j3 largest of j1,j2,j3

      jc=max(j1,j2,j3)
      
      if(jc.eq.j3) then
      ja=j1
      jb=j2
      jc=j3
      ma=m1
      mb=m2
      mc=m3
      elseif(jc.eq.j2) then
      ja=j3
      jb=j1
      jc=j2
      ma=m3
      mb=m1
      mc=m2
      else
      ja=j2
      jb=j3
      jc=j1
      ma=m2
      mb=m3
      mc=m1
      endif

      lm2=-jb
      if(ja+mc-jb.lt.0) lm2=-mc-ja
      lm1=-mc-lm2
      m=lm1+jb+mc+1
      if(ja-jb-mc.lt.0) m=lm1+ja+1
      ALLOCATE(y(m+1),STAT=ierr)
      IF(ierr/=0 .OR.m<=0) STOP 'allocate error in SUB thrj, STOP'
      y=0.0_r8b

      y(1)=0.0_r8b
      y(2)=1.0_r8b
      ss=1.0_r8b
      if(m.eq.1) goto 20

      do 10 n=2,m
      mx=lm1-n+2
      my=lm2+n-2
      alpha=sqrt(real((ja-mx+1)*(ja+mx)*(jb+my+1)*(jb-my),KIND=r8b))
      beta=real(jc*(jc+1)-ja*(ja+1)-jb*(jb+1)-2*mx*my,KIND=r8b)
      gamma=sqrt(real((ja+mx+1)*(ja-mx)*(jb-my+1)*(jb+my),KIND=r8b))
      y(n+1)=(beta*y(n)-gamma*y(n-1))/alpha
      ss=ss+y(n+1)*y(n+1)
10    continue

20    n=lm1-ma+2

      thrj=y(n)*((-1)**(-ja+jb+mc))/sqrt(ss*real(2*jc+1,KIND=r8b))

      DEALLOCATE(y)
      return
  END FUNCTION thrj

  SUBROUTINE fact(n1,n2,a)
!
!       Finds log[n1(factorial)/n2(factorial)], returned in a.
!
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: n1, n2
        REAL,    INTENT(OUT) :: a

        INTEGER :: m1,m2,m,i
        REAL :: b        

	m1=n1
	m2=n2
	if(n1.eq.0) n1=1
	if(n2.eq.0) n2=1
	if(n2.gt.n1) go to 10
	a=-alog(float(n2))
	m=n1-n2+1
	do 5 i=1,m
5       a=a+alog(float(n2+i-1))
	go to 20
10      b=alog(float(n1))
	m=n2-n1+1
	do 15 i=1,m
15      b=b-alog(float(n1+i-1))
	a=b
20      n1=m1
	n2=m2
	return
  END SUBROUTINE fact

  SUBROUTINE lgndr(l,c,s,x,dx)
!
!    computes legendre function x(l,m,theta)
!    theta=colatitude,c=cos(theta),s=sin(theta),l=angular order,
!    sin(theta) restricted so that sin(theta).ge.1.e-7
!    x(1) contains m=0, x(2) contains m=1, x(k+1) contains m=k
!    m=azimuthal(longitudenal) order 0.le.m.le.l
!    dx=dx/dtheta
!    SUBROUTINE originally came from Physics Dept. Princeton through
!    Peter Davis
!    modified to run stably on the Perkin-Elmer by Jeffrey Park 11/23/85
!
!----------------------------------------------------
!    calculate X_lm(theta)= 
!             (-1)^m * sqrt((2l+1)/4*pi) * sqrt[(l-m)!/(l+m)!] * P_lm(cos_theta)
!                                            (B.58)
!                 where P_lm is the "associated lengendre function"
!
!              dX_lm/d_theta
!----------------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(INOUT) :: l          
      REAL*8,   INTENT(IN   ) :: c
      REAL*8,   INTENT(INOUT) :: s
      REAL*8,   INTENT(  OUT) :: x(:),dx(:)  !x(l+1),dx(l+1)

      REAL*8 :: tol, rfpi, root3, boeing

      REAL*8    :: c1,c2,sos,cot,ct,ss,g3,g2,g1,f3,f2,f1,w,v,y,z,d,t,&
                 stom,fac
      INTEGER :: lsave,lp1,i,mp1,m,lpsafex,lpsafe,mmm,maxsin
      tol=1.e-5
      rfpi=0.282094791773880  !sqrt(1/4(pi))
      root3=1.73205080756890
      boeing=0.707106781186550
      x(:)=0.0 ! set all to zero
      dx(:)=0.0 ! set all to zero

      if(s.ge.1.0-tol) s=1.0-tol
      lsave=l
      if(l.lt.0) THEN
        l=-1-l
        WRITE(*,*) 'are you sure l is negtive in SUB lgndr?? STOP'
        STOP
      END IF
      if(l.gt.0) go to 1
      x(1)=rfpi
      dx(1)=0.0
      l=lsave
      return
    1 if(l.ne.1) go to 2
      c1=root3*rfpi
      c2=boeing*c1
      x(1)=c1*c
      x(2)=-c2*s
      dx(1)=-c1*s
      dx(2)=-c2*c
      l=lsave
      return
    2 sos=s
      if(s.lt.tol) s=tol
      cot=c/s
      ct=2.0*c
      ss=s*s
      lp1=l+1
      g3=0.0
      g2=1.0
      f3=0.0
!  evaluate m=l value, sans (sin(theta))**l
      do 100 i=1,l
  100 g2=g2*(1.0-1.0/(2.0*i))
      g2=rfpi*sqrt((2*l+1)*g2)
      f2=l*cot*g2
      x(lp1)=g2
      dx(lp1)=f2
      w=0.0
      v=1.0
      y=2.0*l
      z=y+1.0
      d=sqrt(v*y)
      t=0.0
      mp1=l
      m=l-1
!  these recursions are similar to ordinary m-recursions, but since we
!  have taken the s**m factor out of the xlm's, the recursion has the powers
!  of sin(theta) instead
    3 g1=-(ct*mp1*g2+ss*t*g3)/d
      f1=(mp1*(2.0*s*g2-ct*f2)-t*ss*(f3+cot*g3))/d-cot*g1
      x(mp1)=g1
      dx(mp1)=f1
      if(m.eq.0) go to 4
      mp1=m
      m=m-1
      v=v+1.0
      y=y-1.0
      t=d
      d=sqrt(v*y)
      g3=g2
      g2=g1
      f3=f2
      f2=f1
       go to 3
    4 maxsin=-72.0/log10(s)
!  maxsin is the max exponent of sin(theta) without underflow
      lpsafe=min0(lp1,maxsin)
      stom=1.0
      fac=sign(1.0,(l/2)*2-l+.50)
!  multiply xlm by sin**m
      do 5 m=1,lpsafe
      x(m)=fac*x(m)*stom
      dx(m)=fac*dx(m)*stom
      stom=stom*s
    5 continue
!  set any remaining xlm to zero
      if(maxsin.le.l) then
      mmm=maxsin+1
      do 200 m=mmm,lp1
      x(m)=0.0
  200 dx(m)=0.0
      endif
      s=sos
      l=lsave
      return
  END SUBROUTINE lgndr 

  FUNCTION splbasis(x,x0,dx,np,i)
      IMPLICIT NONE
      REAL :: splbasis
      INTEGER, INTENT(IN) :: np,i
      REAL,    INTENT(IN) :: x,x0,dx
!
      INTEGER :: nx,interval
      REAL    :: xdiff,xd,h,hsq,hcu,xdsq,xdcu,value    

      nx=np-1
      xdiff=x-x0
      interval=1+int((x-x0)/dx)
      if(abs(xdiff-float(nx)*dx).lt.0.0000001) interval=nx
      xd=x-x0-float(interval-1)*dx
      h=1./dx
      hsq=1./dx**2
      hcu=1./dx**3
      xdsq=xd**2
      xdcu=xd**3
!
!---- return the value of the i-th basis element
!
      value=0.
      if(i.eq.0) then
        if(interval.eq.1) then
          value=0.25*hcu*xdcu-1.5*h*xd+1.5
        else if(interval.eq.2) then
          value=-0.25*hcu*xdcu+0.75*hsq*xdsq-0.75*h*xd+0.25
        else
          value=0.
        endif
      else if(i.eq.1) then
        if(interval.eq.1) then
          value=-0.5*hcu*xdcu+1.5*h*xd
        else if(interval.eq.2) then
          value=0.75*hcu*xdcu-1.5*hsq*xdsq+1.
        else if(interval.eq.3) then
          value=-0.25*hcu*xdcu+0.75*hsq*xdsq-0.75*h*xd+0.25
        else
         value=0.
        endif
      else if(i.gt.1.and.i.lt.nx-1) then
        if(interval.eq.i-1) then
          value=0.25*hcu*xdcu
        else if(interval.eq.i) then
          value=-0.75*hcu*xdcu+0.75*hsq*xdsq+0.75*h*xd+0.25
        else if(interval.eq.i+1) then
          value=0.75*hcu*xdcu-1.5*hsq*xdsq+1.
        else if(interval.eq.i+2) then
          value=-0.25*hcu*xdcu+0.75*hsq*xdsq-0.75*h*xd+0.25
        else
          value=0.
        endif
      else if(i.eq.nx-1) then
        if(interval.eq.nx-2) then
          value=0.25*hcu*xdcu
        else if(interval.eq.nx-1) then
          value=-0.75*hcu*xdcu+0.75*hsq*xdsq+0.75*h*xd+0.25
        else if(interval.eq.nx) then
          value=0.5*hcu*xdcu-1.5*hsq*xdsq+1.
        else
          value=0.
        endif
      else if(i.eq.nx) then
        if(interval.eq.nx-1) then
          value=0.25*hcu*xdcu
        else if(interval.eq.nx) then
          value=-0.25*hcu*xdcu+0.75*hsq*xdsq+0.75*h*xd+0.25
        else
          value=0.
        endif
      endif
      splbasis=value
      return
  END FUNCTION splbasis

  SUBROUTINE choles(a,g,b,y,x,n,nono)
!
!      dimension a(n*(n+1)/2),g(n*(n+1)/2)
!      dimension b(n),y(n),x(n)
!
!        a= row-wise p.d. symm. system  n*(n+1)/2
!        g= cholesky storage
!        b= r.h.s. vector               n
!        y= temp. vector
!        x= answer vector
!        n= system dimension
!        nono .gt. 0 is the level at which p.d. failed
!
!        (a,g) and (b,y,x) may be equivalenced.
!
!----------------------------------------------------------

   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: n
   REAL,    INTENT(IN)  :: a(n*(n+1)/2)
   REAL,    INTENT(OUT) :: g(n*(n+1)/2)  
   REAL,    INTENT(IN)  :: b(n)
   REAL,    INTENT(OUT) :: y(n)
   REAL,    INTENT(OUT) :: x(n)
   INTEGER, INTENT(OUT) :: nono

   INTEGER :: i,j,k,jmin,jmax,kmax,kz,kj 
   REAL    :: sg,gkz
!-----
!     first compute cholesky decomposition

      nono=0
      
      if(a(1).le.0.) then
      nono=1
      return
      endif
      
      g(1)=sqrt(a(1))
      y(1)=b(1)/g(1)

      do 400 i=2,n
      
      kz=(i*(i-1))/2
      g(kz+1)=a(kz+1)/g(1)
      sg=g(kz+1)**2
      y(i)=b(i)-g(kz+1)*y(1)
      
      if(i.gt.2) then
      
      jmax=i-1
      
      do 200 j=2,jmax
      
      gkz=a(kz+j)
      kj=(j*(j-1))/2
      kmax=j-1
      
      do 100 k=1,kmax
      gkz=gkz-g(kz+k)*g(kj+k)
100   continue

      g(kz+j)=gkz/g(kj+j)
      y(i)=y(i)-g(kz+j)*y(j)
      sg=sg+g(kz+j)**2
      
200   continue

      endif
      
      gkz=a(kz+i)-sg
      
      if(gkz.le.0.) then
      nono=i
      return
      endif
      
      g(kz+i)=sqrt(gkz)
      y(i)=y(i)/g(kz+i)
      
400   continue

      kz=(n*(n-1))/2
      x(n)=y(n)/g(kz+n)
      if(n.le.1) return

!-----
!     compute solution for particular rhs
      
      do 600 k=2,n
      
      i=n+1-k
      x(i)=y(i)
      jmin=i+1
      
      do 500 j=jmin,n
      kj=(j*(j-1))/2
      x(i)=x(i)-g(kj+i)*x(j)
500   continue

      kz=(i*(i+1))/2
      x(i)=x(i)/g(kz)
      
600   continue

      return
  END SUBROUTINE choles

  SUBROUTINE checknaninf_r4(x,cx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(:)
    CHARACTER(LEN=*) :: cx
   
    REAL, PARAMETER :: hugeconst=HUGE(1.0)
    LOGICAL :: isfinite
   
    isfinite= .NOT. (ANY(isnan(x)) .OR. ANY(ABS(x)>hugeconst))
    IF (.NOT. isfinite) THEN 
    WRITE(*,*) 'NAN or infinite occurs in ', TRIM(cx), ' STOP'
      OPEN(UNIT=102,FILE='naninf.check')
      WRITE(UNIT=102,FMT='(E14.6)') x
      CLOSE(UNIT=102)
      STOP
    END IF
   END SUBROUTINE checknaninf_r4
  
   SUBROUTINE checknaninf_r4_2d(x,cx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(:,:)
    CHARACTER(LEN=*) :: cx
   
    REAL, PARAMETER :: hugeconst=HUGE(1.0)
    LOGICAL :: isfinite
    INTEGER :: n, i
   
    isfinite=.TRUE.
    n=SIZE(x,DIM=2)
    i=1
    DO 
     IF (i > n .OR. isfinite .eqv. .FALSE.) EXIT
     isfinite= .NOT. (ANY(isnan(x)) .OR. &
                      ANY(ABS(x)>hugeconst))
     i=i+1
    END DO

    IF (.NOT. isfinite) THEN
      WRITE(*,*) 'NAN or infinite occurs in ', TRIM(cx), ' STOP'
      OPEN(UNIT=102,FILE='naninf.check')
      WRITE(UNIT=102,FMT='(E14.6)') x
      CLOSE(UNIT=102)
      STOP
     END IF
   END SUBROUTINE checknaninf_r4_2d


   SUBROUTINE checknaninf_cpx4(x,cx)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: x(:)
    CHARACTER(LEN=*) :: cx
   
    REAL, PARAMETER :: hugeconst=HUGE(1.0)
    LOGICAL :: isfinite
   
    isfinite= .NOT. (ANY(isnan(REAL(x))) .OR. &
                     ANY(isnan(AIMAG(x))) .OR. &
                     ANY(ABS(REAL(x))>hugeconst) .OR. &
                     ANY(ABS(AIMAG(x))>hugeconst))
    IF (.NOT. isfinite) THEN 
    WRITE(*,*) 'NAN or infinite occurs in ', TRIM(cx), ' STOP'
      OPEN(UNIT=102,FILE='naninf.check')
      WRITE(UNIT=102,FMT='(2E14.6)') x
      CLOSE(UNIT=102)
      STOP
    END IF
   END SUBROUTINE checknaninf_cpx4

   SUBROUTINE checknaninf_cpx4_2d(x,cx)
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: x(:,:)
    CHARACTER(LEN=*) :: cx

    REAL, PARAMETER :: hugeconst=HUGE(1.0)
    REAL, DIMENSION(SIZE(x,DIM=1)) :: xreal, ximag
    INTEGER :: i, n 
    LOGICAL :: isfinite

    isfinite=.TRUE.
    n=SIZE(x,DIM=2)
    i=1
    DO 
     IF (i > n .OR. isfinite .eqv. .FALSE.) EXIT
     xreal=REAL(x(:,i))
     ximag=AIMAG(x(:,i))
     isfinite= .NOT. (ANY(isnan(xreal)) .OR. &
                      ANY(isnan(ximag)) .OR. &
                      ANY(ABS(xreal)>hugeconst) .OR. &
                      ANY(ABS(ximag)>hugeconst))
      i=i+1
    END DO
    IF (.NOT. isfinite) THEN 
      WRITE(*,*) 'NAN or infinite occurs in ', TRIM(cx), ' STOP!'
      OPEN(UNIT=102,FILE='naninf.check')
      WRITE(UNIT=102,FMT='(2E14.6)') x
      CLOSE(UNIT=102)
      STOP
    END IF
   END SUBROUTINE checknaninf_cpx4_2d

  SUBROUTINE intgrl_dp(sum1,r,nir,ner,f,s1,s2,s3,sum1r)
!
!       Computes the integral of f[i]*r[i]*r[i] from i=nir to i=ner for
!       radii values as in model stored in f
!       sum1: result at r=r(ner)
!       sum1r: result as a function of r(nir:ner), but the dimension
!              of the array is the same as r (i.e. NR), rather than
!              ner-nir+1 

        USE modellayer_module, ONLY:NR,kdis
        IMPLICIT NONE
        REAL(KIND=8),   INTENT(IN    ) :: r(:),f(:)
        INTEGER,INTENT(IN    ) :: nir,ner
        REAL(KIND=8),   INTENT(   OUT) :: sum1
        REAL(KIND=8),   INTENT(   OUT) :: s1(:),s2(:),s3(:)
        REAL(KIND=8),   INTENT(   OUT), OPTIONAL :: sum1r(:)
  
        INTEGER :: ndis
                  ! DIMENSION(NR)
        REAL(KIND=8), DIMENSION(:),ALLOCATABLE :: yprime, sum0  
        REAL(KIND=8) :: third,fifth,sixth
        REAL(KIND=8) :: rji,dsum 
        INTEGER :: n 
 
        INTEGER :: ierr,i,j,nir1

       third = 1.d0/3.d0
       fifth = 1.d0/5.d0
       sixth = 1.d0/6.d0

       n=SIZE(r)
       IF (n/=SIZE(f).OR.n/=NR) STOP 'dimension inconsistent in SUB intgrl, STOP'
       ALLOCATE(yprime(n),sum0(n),STAT=ierr)
       IF (ierr/=0.OR.n==0) STOP 'allocate error 1 in SUB intgrl, STOP'

                        ! use discontinuities r-
       IF (.NOT. ASSOCIATED(kdis)) STOP 'kdis has to be defined previously, STOP' 
       ndis=SIZE(kdis)
   
                        ! calculate cubic splines 
	call deriv(f,yprime,n,r,ndis,kdis,s1,s2,s3)

        sum0=0.0d0
	nir1 = nir + 1
	DO i=nir1,ner
	  j = i-1
	  rji = r(i) - r(j)
               ! Note that "r(j)**2" at the beginning of RHS represent 
               !  the integral is multiplied by r^2
          dsum=r(j)*r(j)*rji*(f(j)+rji*(0.5d0*s1(j)+rji*(third*s2(j)+rji* &
                  0.25d0*s3(j))))+2.0d0*r(j)*rji*rji*(0.5d0*f(j)+rji*(third*s1(j)+rji* &
                 (0.25d0*s2(j)+rji*fifth*s3(j))))+rji*rji*rji*(third*f(j)+rji* &
                 (0.25d0*s1(j)+rji*(fifth*s2(j)+rji*sixth*s3(j))))
          sum0(i)=sum0(j)+dsum 
        END DO

        sum1=sum0(ner)
        IF (PRESENT(sum1r)) THEN
          sum1r=sum0
          IF (nir/=1) THEN
             WRITE(*,*) 'Wanring!! beginning entry in r for SUB intgrl is not "1" '
          END IF
        END IF
       
        DEALLOCATE(yprime,sum0)

	return
  END SUBROUTINE intgrl_dp 

  SUBROUTINE deriv_dp(y,yprime,n,r,ndis,kdis,s1,s2,s3)
!     calculate the coefficients (s1,s2,s3) of cubic spline
!     as a function of r
!
!     yprime=s1
!    
!     Numerical error occurrs in the portion of y=constant,
!     in which s1,s2,s3 have to be zeor but the numerical results are not.
!     This problem can be resolved by changing data type of variables
!     from REAL*4 to REAL*8.
!     Even though inaccurate s1,s2,s3 do not affect integration (SUB intgrl).
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: ndis, kdis(:)  !kdis(ndis)
      REAL*8,  INTENT(IN) :: y(:), r(:)     !y(n),r(n)
      REAL*8,  INTENT(OUT):: yprime(n),s1(n),s2(n),s3(n)

      REAL*8 :: yy(3), f(3,n)

      INTEGER :: j1,j2,nd,ndp,i,j,k
      REAL*8    :: h,h2,a0,b0,b1,ha,h2a,h3a,hb,h2b,s13,s21,s32,y0
      equivalence (yy(1),y0)
      yy=0.0

      ndp=ndis+1
      do 3 nd=1,ndp
      if(nd.eq.1) go to 4
      if(nd.eq.ndp) go to 5
      j1=kdis(nd-1)+1
      j2=kdis(nd)-2
      go to 6
    4 j1=1
      j2=kdis(1)-2
      go to 6
    5 j1=kdis(ndis)+1
      j2=n-2
    6 if((j2+1-j1).gt.0) go to 11
      j2=j2+2
      y0=(y(j2)-y(j1))/(r(j2)-r(j1))
	s1(j1)=yy(1)
	s1(j2)=yy(1)
	s2(j1)=yy(2)
	s2(j2)=yy(2)
	s3(j1)=yy(3)
	s3(j2)=yy(3)
      go to 3
   11 a0=0.d0
      if(j1.eq.1) go to 7
      h=r(j1+1)-r(j1)
      h2=r(j1+2)-r(j1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
      b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
      go to 8
    7 b0=0.d0
    8 b1=b0
      if(j2 .gt. n)write(0,'("error:deriv:j2= ",i5)')j2
      do 1 i=j1,j2
      h=r(i+1)-r(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.d0*a0
      h3a=2.d0*h-3.d0*a0
      h2b=h2*b0
      s1(i)=h2/ha
      s2(i)=-ha/(h2a*h2)
      s3(i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
      a0=s3(i)
    1 b0=f(3,i)
      i=j2+1
      h=r(i+1)-r(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.d0*h-a0)
      s1(i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=r(j2)-r(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      s3(i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
!     the following statements were expanded to prevent register overflow on unix
      s13=s1(i)*s3(i)
      s2(i)=f(1,i)-s13
      do 2 j=j1,j2
      k=i-1
      s32=s3(k)*s2(i)
      s1(i)=f(3,k)-s32
      s21=s2(k)*s1(i)
      s3(k)=f(2,k)-s21
      s13=s1(k)*s3(k)
      s2(k)=f(1,k)-s13
    2 i=k
      s1(i)=b1
      j2=j2+2
	s1(j2)=yy(1)
	s2(j2)=yy(2)
	s3(j2)=yy(3)
    3 continue

      do 20 i=1,n
   20 yprime(i)=s1(i)
      return
  END SUBROUTINE deriv_dp

  SUBROUTINE RconvertC_1D(iforward,ndim,nmodes,ll,aa)
!
!  convert 1D matrix in terms of real spherical basis to complex spherical basis
!   Supposing that r_real is 1D matrix, U is transformation matrix from real
!   to complex basis, then
!
!    r_{complex}= U * r_{real}                          [iforward=1]
!
!    r_{real}= U^h * r_{complex},  h denotes hermitian  [ifortward=-1]
!
    IMPLICIT NONE
    INTEGER,               INTENT(IN   ) :: iforward
    INTEGER,               INTENT(IN   ) :: ndim
    INTEGER,               INTENT(IN   ) :: nmodes
    INTEGER, DIMENSION(:), INTENT(IN   ) :: ll  
    COMPLEX, DIMENSION(:), INTENT(INOUT) :: aa  !DIMENSION((2*ll_{k}+1),k=1,nmodes) 

    INTEGER              :: ndim_subU
    COMPLEX, ALLOCATABLE :: usub(:,:)
 
    INTEGER :: i, ierr, j1, j2

    j1=1
    DO i = 1,nmodes 
      ndim_subU=2*ll(i)+1
      j2=j1+ndim_subU-1
      ALLOCATE(usub(ndim_subU,ndim_subU),STAT=ierr)
      IF (ierr/=0) STOP 'allocate error in SUB:[RconvertC_1D], STOP'

                 ! create transformation matrix U for each multiplet 
      CALL create_Usub(ll(i),iforward,usub)

                 ! perform transformation on 1D-array aa
      aa(j1:j2)=MATMUL(usub,aa(j1:j2)) 

      j1=j2+1
      DEALLOCATE(usub)
    END DO

    IF (j2/=ndim) STOP &
      'dimension inconstent in sub [RconvertC_1D], STOP'    
    
  END SUBROUTINE RconvertC_1D

  SUBROUTINE RconvertC_2D(iforward,ndim,nmodes,ll,aa,istruczero)
!
!  convert 2D matrix in terms of real spherical basis to complex spherical basis
!   Supposing that R_real is 2D matrix, U is transformation matrix from real
!   to complex basis, then
!
!    R_{complex}= U * R_{real} * U^h                    [iforward=1]
!
!    R_{real}= U^h * R_{complex} * U                  [ifortward=-1]
!   
!    h denotes hermitian
!
    IMPLICIT NONE
    INTEGER               ,INTENT(IN   ) :: iforward                              
    INTEGER               ,INTENT(IN   ) :: ndim
    INTEGER               ,INTENT(IN   ) :: nmodes
    INTEGER,DIMENSION(:)  ,INTENT(IN   ) :: ll
    COMPLEX,DIMENSION(:,:),INTENT(INOUT) :: aa 
                    !truncate real or/and imaginary part to zero
    TYPE(istruczero_type), OPTIONAL,     INTENT(IN)    :: istruczero

    COMPLEX, ALLOCATABLE :: usub(:,:)
    COMPLEX, DIMENSION(ndim,ndim) :: uall
    COMPLEX, DIMENSION(ndim,ndim) :: ctmp, c2tmp
    INTEGER :: ndim_subU
    INTEGER :: i,k,j1,j2,ierr
    LOGICAL :: istest

    IF (iforward/=1 .AND. iforward/=-1) &
    STOP 'iforward is out of range in SUB[RconvertC_2D], STOP!'

    ! for the 2D matrix, we cannot simply generate a sub transformation matrix,
    ! i.e.(D.176). Instead, a combination of sub transformation matrix is
    ! requried, saying uall (D.175)
    uall=CMPLX(0.0,0.0)
    ctmp=CMPLX(0.0,0.0)
    c2tmp=CMPLX(0.0,0.0)

      !j1,j2 keep tracking of entry in huge matrix uall
    j1=1
    DO i = 1,nmodes
      ndim_subU=2*ll(i)+1
      j2=j1+ndim_subU-1
      ALLOCATE(usub(ndim_subU,ndim_subU),STAT=ierr)
      IF (ierr/=0) STOP 'allocate error in SUB:[RconvertC_1D], STOP'

                 ! create transformation matrix U (iforward=1) 
                 ! or U^h (iforward=-1) for each multiplet 
      CALL create_Usub(ll(i),iforward,usub)

                 ! assign sub matrix U/U^h to a huge matrix, uall
      uall(j1:j2,j1:j2)=usub

      j1=j2+1
      DEALLOCATE(usub)
    END DO

     ! if iforward=1:  Uall= U:  U*aa*U^h
     ! if iforward=-1: Uall=U^h: U^h*aa*U
      DO i=1,ndim
        DO k=1,ndim
          ctmp(i,k)=SUM(uall(i,1:ndim)*aa(1:ndim,k))
        END DO
      END DO
      DO i=1,ndim
        DO k=1,ndim
          c2tmp(i,k)=SUM(ctmp(i,1:ndim)*CONJG(uall(k,1:ndim)))
!          aa(i,k)=SUM(ctmp(i,1:ndim)*CONJG(uall(k,1:ndim)))
        END DO
      END DO

      IF (PRESENT(istruczero) .AND. iforward==-1) THEN
        ! set the imaginary part to be zero or not
        IF (istruczero%imag)THEN
          WRITE(*,*) 'ola1'
          WHERE(AIMAG(c2tmp) /= 0.0 )
            c2tmp = CMPLX(REAL(c2tmp),0.0)
          END WHERE
          IF (ANY(AIMAG(c2tmp)/=0.0)) STOP '[RconvertC_2D]: fail 1, STOP'
        END IF
 
        ! set the real part to be zero or not
        IF (istruczero%real)THEN
          WRITE(*,*) 'ola2'
          WHERE(REAL(c2tmp) /= 0.0 )
            c2tmp = CMPLX(0.0,AIMAG(c2tmp))
          END WHERE
          IF (ANY(REAL(c2tmp)/=0.0)) STOP '[RconvertC_2D]: fail 2, STOP'
        END IF
      END IF

      ! 
      aa(1:ndim,1:ndim)=c2tmp

  END SUBROUTINE RconvertC_2D

  SUBROUTINE RconvertC_2D_savemeo(iforward,ndim,nmodes,ll,aa,istruczero)
!>
!>  Like RconvertC_2D except using another way to doing mulitiplication
!>  (save memoery)
!>
!>  convert 2D matrix in terms of real spherical basis to complex spherical basis
!>   Supposing that R_real is 2D matrix, U is transformation matrix from real
!>   to complex basis, then
!>
!>    R_{complex}= U * R_{real} * U^h                    [iforward=1]
!>
!>    R_{real}= U^h * R_{complex} * U                  [ifortward=-1]
!>   
!>    h denotes hermitian
!>
    IMPLICIT NONE
    INTEGER               ,INTENT(IN   ) :: iforward                              
    INTEGER               ,INTENT(IN   ) :: ndim
    INTEGER               ,INTENT(IN   ) :: nmodes
    INTEGER,DIMENSION(:)  ,INTENT(IN   ) :: ll
    COMPLEX,DIMENSION(:,:),INTENT(INOUT) :: aa 
                    !truncate real or/and imaginary part to zero
    TYPE(istruczero_type), OPTIONAL,     INTENT(IN)    :: istruczero

    INTEGER :: ndims(nmodes)
    COMPLEX, ALLOCATABLE :: usub(:,:)
    COMPLEX, ALLOCATABLE :: ctmp(:,:),rightmat(:,:)
    COMPLEX, ALLOCATABLE :: uall(:,:,:)
    INTEGER :: lmax,ndim_max
    INTEGER :: imode,jmode,ndimi,ndimj,ierr,i1,i2,j1,j2

    IF (iforward/=1 .AND. iforward/=-1) &
    STOP 'iforward is out of range in SUB[RconvertC_2D_savemeo], STOP!'

    ! for the 2D matrix, we cannot simply generate a sub transformation matrix,
    ! i.e.(D.176). Instead, a combination of sub transformation matrix is
    ! requried, saying uall (D.175)
    lmax=MAXVAL(ll(1:nmodes))
    ndim_max=2*lmax+1
    ALLOCATE(uall(ndim_max,ndim_max,nmodes),&
             ctmp(ndim_max,ndim_max),&
             rightmat(ndim_max,ndim_max),STAT=ierr)
    IF (ierr/=0) STOP 'allocate error 1 in SUB [RconvertC_2D_savemeo], STOP'
    uall=CMPLX(0.0,0.0)

    DO imode=1,nmodes
      ndimi=2*ll(imode)+1
      ALLOCATE(usub(ndimi,ndimi),STAT=ierr)
      IF (ierr/=0) STOP 'allocate error 2 in SUB:[RconvertC_2D_savemeo], STOP'

              ! create transformation matrix U (iforward=1) 
              ! or U^h (iforward=-1) for each multiplet 
      CALL create_Usub(ll(imode),iforward,usub)

             ! assign sub matrix U/U^h to a huge matrix, uall
      uall(1:ndimi,1:ndimi,imode)=usub

      ndims(imode)=ndimi
      DEALLOCATE(usub)
    END DO
!    WRITE(*,*) 'ndims:', ndims(1:nmodes)

     !j1,j2,i1,i2 keep tracking of entry in global matrix aa
    j1=1 
    DO jmode=1,nmodes
      ndimj=ndims(jmode)
      j2=j1+ndimj-1
      rightmat(1:ndimj,1:ndimj)=TRANSPOSE(CONJG(uall(1:ndimj,1:ndimj,jmode)))

      i1=1
      DO imode=1,nmodes
        ndimi=ndims(imode)
        i2=i1+ndimi-1 
!        WRITE(*,*) 'i,j:', i1,i2,j1,j2 
   
        ctmp(1:ndimi,1:ndimj)= MATMUL(aa(i1:i2,j1:j2),rightmat(1:ndimj,1:ndimj))
        aa(i1:i2,j1:j2)=MATMUL(uall(1:ndimi,1:ndimi,imode),&
                               ctmp(1:ndimi,1:ndimj))
 
        i1=i2+1
        ctmp=CMPLX(0.0,0.0)
      END DO
      j1=j2+1
      rightmat=CMPLX(0.0,0.0)
    END DO

    DEALLOCATE(ctmp,uall,rightmat)

    IF (PRESENT(istruczero) .AND. iforward==-1) THEN
        ! set the imaginary part to be zero or not
        IF (istruczero%imag)THEN
          WHERE(AIMAG(aa) /= 0.0 )
            aa = CMPLX(REAL(aa),0.0)
          END WHERE
          IF (ANY(AIMAG(aa)/=0.0)) STOP '[RconvertC_2D_savemeo]: fail 1, STOP'
        END IF
 
        ! set the real part to be zero or not
        IF (istruczero%real)THEN
          WHERE(REAL(aa) /= 0.0 )
            aa = CMPLX(0.0,AIMAG(aa))
          END WHERE
          IF (ANY(REAL(aa)/=0.0)) STOP '[RconvertC_2D_savemeo]: fail 2, STOP'
        END IF
      END IF
  END SUBROUTINE RconvertC_2D_savemeo

  SUBROUTINE create_Usub(ll,iforward,uu)
!  
!  create transformation matrix from real to complex spherical basis
!  or vice versa 
!
!
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: ll
    INTEGER, INTENT(IN)  :: iforward ! =1:U, =-1:U^h, h denotes hermitian
    COMPLEX, INTENT(OUT) :: uu(:,:)

    REAL :: rtmp
    INTEGER :: m,i1,i2
    REAL :: twosq
  
    i1=2*ll+1 
    IF ((SIZE(uu,DIM=1)/=i1).OR. &
        (SIZE(uu,DIM=2)/=i1)) STOP &
       'size error in sub [create_Usub], STOP'
    twosq=1.0/SQRT(2.0)
    uu=CMPLX(0.0,0.0)

       
    SELECT CASE (iforward)
      CASE (1)
                                   ! create U (D.176) in D.T. 1998
        DO m=-ll,ll
          i1= m+ll+1
          i2=-m+ll+1
          IF (m<0) THEN
            rtmp=twosq*(-1)**(-m)
            uu(i1,i1)=CMPLX(rtmp,0.0)   !left top corner
            uu(i1,i2)=CMPLX(0.0,rtmp)   !right top corner
          ELSE IF (m==0) THEN
            uu(i1,i1)=CMPLX(1.0,0.0)      !center
          ELSE IF (m>0) THEN
            uu(i1,i1)=CMPLX(0.0,-twosq) !right bottom corner
            uu(i1,i2)=CMPLX(twosq,0.0)  !left bottom corner
          END IF
        END DO

      CASE (-1)
                                    ! U^h
        DO m=-ll,ll
          i1= m+ll+1
          i2=-m+ll+1
          IF (m<0) THEN
            rtmp=twosq*(-1)**(-m)
            uu(i1,i1)=CMPLX(rtmp,0.0)   !left top corner
            uu(i2,i1)=CMPLX(0.0,-rtmp)  !left bottom corner
          ELSE IF (m==0) THEN
            uu(i1,i1)=CMPLX(1.0,0.0)      !center
          ELSE IF (m>0) THEN
            uu(i1,i1)=CMPLX(0.0,twosq)  !right bottom corner
            uu(i2,i1)=CMPLX(twosq,0.0)  !right top corner
          END IF
        END DO

      CASE DEFAULT
        STOP 'iforward is out of range in SUB:[create_Usub], STOP'

    END SELECT    

  END SUBROUTINE create_Usub

  SUBROUTINE check_BiOrth_eigvec(QR,QL,dist)
!----------------------------------------------------------
! calculate bi-orthogoanlity of QR and QL, i.e. (QL)^h * QR
! and meausre the distance from identical matrix
!
! {| (QL^h)*QR - I |^2 }/n^2
! For FC case, (QL^h*QL) has not to be idential matrix.
!
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: QR(:,:), QL(:,:)
    REAL, INTENT(OUT):: dist(2)

    COMPLEX, DIMENSION(SIZE(QR,1),SIZE(QR,2)) :: I_pred, II
    INTEGER :: ndim
    INTEGER :: i

    ndim=SIZE(QR,1)
    I_pred=MATMUL(CONJG(TRANSPOSE(QL)),QR)
    II=eyes(ndim)
      !average distance for the entire matrix
    II=(ABS(I_pred-II))**2
    dist(1)=SUM(II)/ndim**2 
      !average distance for the diagonal terms of the entire matrix
    dist(2)=0.0
    DO i=1,ndim
      dist(2)=dist(2)+II(i,i)
    END DO
    dist(2)=dist(2)/ndim 
  CONTAINS
    FUNCTION eyes(ndim)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ndim
      COMPLEX, DIMENSION(ndim,ndim) :: eyes

      INTEGER :: i

      eyes=CMPLX(0.0,0.0)
      DO i=1,ndim
          eyes(i,i)=CMPLX(1.0,0.0)
      END DO
    END FUNCTION eyes
  END SUBROUTINE check_BiOrth_eigvec

  SUBROUTINE norm_SEP_eigvec(y,x,norm)
!>---------------------------------------------------------------------
!> normalize eigenvectors for Standard Eigenvalue Problem, which
!>  has the form,
!>
!>        y^h * H * x = lambda,
!>       
!>     where h denotes hermitian
!>           lambda: eigenvalues
!>           x: right eigenvector  (QR)
!>           y: left eigenvector   (QL)
!>      
!>--------------------------------------------------------------------- 
    IMPLICIT NONE
    COMPLEX, DIMENSION(:), INTENT(INOUT) :: y, x
    COMPLEX,               INTENT(  OUT) :: norm
    norm=DOT_PRODUCT(y,x)
    x=x/norm
  END SUBROUTINE norm_SEP_eigvec

  SUBROUTINE inquire_matrix_type(A,itype)
!>--------------------------------------------------------------------
!> inquire the matrix type
!>  itype =1: real symmetric
!>        =2: hermitian
!>        =3: complex symmetric
!>        =4: general complex
!>---------------------------------------------------------------------
    IMPLICIT NONE
    COMPLEX, INTENT(IN) :: A(:,:)
    INTEGER, INTENT(OUT):: itype

    COMPLEX :: B(SIZE(A,1),SIZE(A,2)), diagA(SIZE(A,1))
    INTEGER :: n
    INTEGER :: i, j
    LOGICAL :: isreal

    isreal=.FALSE.
    n=SIZE(A,1)
    IF( n/=SIZE(A,2) ) STOP &
           'A should be a square matrix in SUB [inquire_matrix]'

    DO i=1,n
      diagA(i)=A(i,i)
    END DO

    isreal = ALL (AIMAG(A)==0.0 )
    B=TRANSPOSE(A)
    IF (isreal) THEN
      IF (ALL(B == A) ) THEN
         !real symmetric
        itype =1
      ELSE 
        B=CONJG(B)
        IF (ALL(B==A)) THEN 
          !hermitian
          itype = 2
        END IF
      END IF
    ELSE
      IF (ALL(B==A)) THEN
        ! complex symmetric
        itype = 3
      ELSE
        ! general complex
        itype = 4
      END IF
    END IF
 
  END SUBROUTINE inquire_matrix_type

  SUBROUTINE check_transpose(A,isT,isW)
    ! real and imaginary part of V_total have to be symmetric, respectively
      IMPLICIT NONE
      COMPLEX, DIMENSION(:,:), INTENT(IN   ) :: A
      LOGICAL,                 INTENT(  OUT) :: isT(2)
      LOGICAL,                 INTENT(IN   ), OPTIONAL :: isW
    
      REAL, DIMENSION(SIZE(A,1),SIZE(A,2)) :: rtmpa
      LOGICAL                  :: isnormal
      INTEGER :: i,j,n

      isT=.TRUE.
      isnormal=.TRUE.
      IF (PRESENT(isW)) THEN
        IF (isW) THEN
          isnormal=.FALSE.
        END IF
      END IF
    
      rtmpa=REAL(A)
      IF (ANY(rtmpa/=TRANSPOSE(rtmpa))) THEN
        isT(1)=.FALSE.
      END IF

      rtmpa=AIMAG(A)
      IF (isnormal) THEN
        IF (ANY(rtmpa/=TRANSPOSE(rtmpa))) THEN
          isT(2)=.FALSE.
        END IF
      ELSE
        IF (ANY(rtmpa/=-TRANSPOSE(rtmpa))) THEN
          isT(2)=.FALSE.
!          WRITE(*,*) 'number:', COUNT(rtmpa/=-TRANSPOSE(rtmpa))
!          n=SIZE(rtmpa,1)
!          DO i=1,n
!            DO j=i,n
!              IF (rtmpa(i,j)/=-rtmpa(j,i)) THEN
!                WRITE(*,*) 'i,j,pos:', i,j,rtmpa(i,j),-rtmpa(j,i)
!              END IF
!            END DO
!          END DO
        END IF
      END IF

    END SUBROUTINE check_transpose

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! <SA> ADDING DIRAC DELTA FUNCTIONAL FOR ADJOINT PURPOSES
    ! subroutine dirac_delta(center_value, increment, size, dirac)
    !   !
    !   implicit none
    !   !
    !   integer, intent(in) :: size
    !   real*8, intent(in) :: center_value, increment
    !   !
    !   complex*16, intent(out) :: dirac(:)
    !   !
    !   real*8 :: lower_limit, upper_limit, integral
    !   real*8 :: dirac_re(size), dirac_im(size)
    !   integer :: i
      
    !   ! Find index at which dirac should be centered
    !   lower_limit = center_value - increment * (size - 1) / 2.0
    !   upper_limit = center_value + increment * (size - 1) / 2.0
     
    !   ! Nullify dirac array
    !   dirac_re(:) = 0.d0
    !   dirac_im(:) = 0.d0
      
    !   ! Populate dirac
    !   print*, center_value, increment, 1.0/increment, lower_limit, upper_limit
      
    !   ! dirac_re = merge(1.0 / increment, dirac_re, &
    !    !     (lower_limit <= (0.0 + epsilon(1.0))) .and. ((0.0 - epsilon(1.0)) <= upper_limit))

    !   do i = 1, size
    !      if (abs((i - 1) * increment - center_value) < increment / 2.0) then
    !         dirac_re(i) = 1.0 / increment
    !      else
    !         dirac_re(i) = 0.0
    !      end if
    !   end do


    !   ! Just checking it integrates to 1
    !   integral = sum(dirac_re) * increment
    !   write(*,*) "Integral of dirac delta is (should be one)", integral 
    !   !
    !   dirac = cmplx(dirac_re, dirac_im)
    ! end subroutine dirac_delta
    ! </SA>

    subroutine dirac_delta(index, increment, ep, size, dirac)
      !
      implicit none
      !
      integer, intent(in) :: index, size
      real*8, intent(in) :: increment, ep
      !
      complex*16, intent(out) :: dirac(:)
      !
      real*8 :: integral
      real*8 :: dirac_re(size), dirac_im(size)
      !integer :: i
      
      ! Find index at which dirac should be centered
      ! lower_limit = center_value - increment * (size - 1) / 2.0
      ! upper_limit = center_value + increment * (size - 1) / 2.0
     
      ! Nullify dirac array
      dirac_re(:) = 0.d0
      dirac_im(:) = 0.d0
      
      ! Populate dirac
      ! print*, center_value, increment, 1.0/increment, lower_limit, upper_limit
      
      ! dirac_re = merge(1.0 / increment, dirac_re, &
       !     (lower_limit <= (0.0 + epsilon(1.0))) .and. ((0.0 - epsilon(1.0)) <= upper_limit))

      ! do i = 1, size
      !    if (abs((i - 1) * increment - center_value) < increment / 2.0) then
      !       dirac_re(i) = 1.0 / increment
      !    else
      !       dirac_re(i) = 0.0
      !    end if

      ! end do
      print*, index, 1.d0/(2.d0*increment)
      dirac_re(index) = 1./(2.d0*increment)
      dirac_im(index) = -1.d0/(2.d0*increment)
      ! dirac_im(index) = 0.d0

      ! Just checking it integrates to 1
      integral = sum(dirac_re) * increment
      write(*,*) "Integral of dirac delta is (should be 0.5)", integral 
      !
      dirac = cmplx(dirac_re, dirac_im)
    end subroutine dirac_delta

END MODULE futil_module
