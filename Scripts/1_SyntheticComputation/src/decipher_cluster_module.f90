MODULE decipher_cluster_module

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: decipher_cluster

CONTAINS

  SUBROUTINE decipher_cluster(cluster,&
       nmodes,ndim,n,type,l,om,gma,Qval)

    ! HL MAR 2018

    ! Simplified version (only 1 cluster-- full coupling!)
    !
    !     INTENT(IN) :: cluster (0T3-0S4-1S2)
    !     INTENT(OUT) :: nmodes  (3)
    !                    ndim: dimension of splitting matrix 
    !                           sum_i { 2*l_i+1 }
    !                    n, type, l  
    !                    om : eigenfrequency in rad/s
    !
    !     Note that "type" in current version can only be "1" (spheroidal
    !                     and radial) or "0" (toroidal).

    USE modellayer_module, ONLY : NR
    USE futil_module, ONLY : get_mode_fun
    USE eigfun_module 
    IMPLICIT NONE
    include "coupling.h"

    CHARACTER(LEN=*), INTENT(IN) :: cluster
    INTEGER, INTENT(OUT) :: nmodes,ndim
    INTEGER, INTENT(OUT) :: n(:),type(:),l(:)  ! DIMENSION(MMAX)
    REAL*8,    INTENT(OUT) :: om(:),gma(:),Qval(:)! DIMENSION(MMAX)

    LOGICAL :: issame
    INTEGER :: lcluster,m1,m2,lmode,ind,i,ierr
    INTEGER :: lnblnk
    REAL*8    :: gv,av,ahs,aht,ap
    REAL*8, DIMENSION(MMAX) :: q
    CHARACTER(LEN=10) :: mode
    CHARACTER(LEN=1)  :: ctype
    ! Dimension: NR
    REAL*8, POINTER :: r(:), u(:), du(:), v(:), dv(:), w(:), dw(:),&
         p(:), dp(:)

    lcluster=lnblnk(cluster)
    nmodes=0
    m1=1
    m2=m1+2
    ndim=0
    
    ! set up n, l, type, nmodes
    do while(m2.le.lcluster)
       nmodes=nmodes+1
       do while((cluster(m2:m2).ne.'-').and.(m2.le.lcluster))
          m2=m2+1
       enddo
       m2=m2-1
       mode=cluster(m1:m2)
       lmode=lnblnk(mode)
       if(index(mode(1:lmode),'S').gt.0) then
          type(nmodes)=1
          ind=index(mode(1:lmode),'S')
       else
          type(nmodes)=0
          ind=index(mode(1:lmode),'T')
       endif
       if(ind.eq.0) then
          write(6,"('Unknown mode type')")
          call exit(1)
       endif
       n(nmodes)=0
       do i=1,ind-1
          n(nmodes)=n(nmodes)+idigit(mode(i:i))*10**(ind-1-i)
       enddo
       l(nmodes)=0
       do i=ind+1,lmode
          l(nmodes)=l(nmodes)+idigit(mode(i:i))*10**(lmode-i)
       enddo
       if(nmodes.gt.MMAX) then
          write(6,"('number of modes in cluster exceeds ',i3)") MMAX
          call exit(1)
       endif
       m1=m2+2
       m2=m1+2
    enddo
    
    ! read eigenfunctions
    CALL allocate_eigfun(NR,nmodes) ! Allocate ONLY

    DO i=1,nmodes
       if(type(i) .eq. 0) then
          ctype = 't'
       else
          ctype = 's'
       endif
       
      NULLIFY(r,u,du,v,dv,w,dw,p,dp)

      call get_mode_fun(n(i),ctype,l(i),om(i),&
          q(i),av,ahs,aht,ap,r,u,du,v,dv,w,dw,p,dp)
       
      IF (SIZE(r)/=SIZE(eigfun_radius) .OR. SIZE(dw)/=SIZE(dw)) &
            STOP 'dimension error in SUB decipher_cluster, STOP'
       eigfun_radius(:) =r
       eigfun_u( :,i)   =u
       eigfun_du(:,i)   =du
       eigfun_v( :,i)   =v    ! v/k
       eigfun_dv(:,i)   =dv   ! dv/k
       eigfun_w( :,i)   =w    ! w/k
       eigfun_dw(:,i)   =dw   ! dw/k
       eigfun_p( :,i)   =p
       eigfun_dp(:,i)   =dp  
       eiginfo(i)%w=om(i) 
       eiginfo(i)%q=q(i)
       eiginfo(i)%n=n(i)
       eiginfo(i)%l=l(i)
       eiginfo(i)%itype=type(i)
       eigac(i)%u=av
       eigac(i)%v=ahs
       eigac(i)%w=aht
       eigac(i)%p=ap
    
!!!! no need to have a maximum bound on modes read in:
!       if(l(nmodes).gt.LMAX) then
!          write(6,"('angular degree of mode exceeds ',i3)") LMAX
!          call exit(1)
!       endif
!!!! ----------------------------------------------------!!!!!!

       ndim=ndim+(2*l(i)+1)
    END DO

    do i=1,nmodes
       gma(i) = 0.5*om(i)/q(i)
       Qval(i) = q(i)
    enddo

    DEALLOCATE(r,u,du,v,dv,w,dw,p,dp)
    return

  CONTAINS
    FUNCTION idigit(ch)
!  
      INTEGER :: idigit
      CHARACTER(LEN=1) :: ch
      !
      if(ch.eq.'0') then
         idigit=0
         return 
      else if(ch.eq.'1') then
         idigit=1
         return 
      else if(ch.eq.'2') then
         idigit=2
         return 
      else if(ch.eq.'3') then
         idigit=3
         return 
      else if(ch.eq.'4') then
         idigit=4
         return 
      else if(ch.eq.'5') then
         idigit=5
         return 
      else if(ch.eq.'6') then
         idigit=6
         return 
      else if(ch.eq.'7') then
         idigit=7
         return 
      else if(ch.eq.'8') then
         idigit=8
         return 
      else if(ch.eq.'9') then
         idigit=9
         return 
      else
        write(6,"('not an iteger')")
        call exit(1)
     endif
     return
   END FUNCTION idigit

 END SUBROUTINE decipher_cluster
 
END MODULE decipher_cluster_module
