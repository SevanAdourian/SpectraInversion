MODULE splitting_seismo_module

  ! ================================================================ !
  ! Build all coupling matrices                                      !
  ! ================================================================ !
  ! MD/HL -- 
  ! This builds the coupled matrices *without* the forcing 
  ! frequency and so must be saved in individual matrices
 
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  coupling_matrix,&
             woodhouse_kernels,&
             woodhouse_kernels_store,&
             adjoint_kernels,&
             adjoint_kernels_int,&
             adjoint_kernels_test,&
             b

CONTAINS

  SUBROUTINE coupling_matrix(sbuild0,sbuild1,is_rotate,ndim,nmode,&
                             swd0,swd1,nmwd,iwdhs,&
                             K_mu_wh,K_kp_wh,K_Tro_wh,K_Vro_wh,&
                             K_Td_wh,K_Vd_wh,&
                             V_mat,W_mat,T_mat,dV_mat)

    USE futil_module, ONLY : intgrl,thrj
    USE prem_module
    USE mantle_model_module_DS
    USE modellayer_module, ONLY : NR,NOC
    USE eigfun_module, ONLY : r=>eigfun_radius,eiginfo
    USE anelastic_module

    IMPLICIT NONE

    include 'coupling.h'
    ! inputs
    logical, intent(in) :: is_rotate
    integer, intent(in) :: sbuild0,sbuild1,ndim,nmode
    integer, intent(in) :: swd0,swd1,nmwd
    ! integer, dimension(0:1,0:120,0:120), intent(in) :: iwdhs
    ! complex*16, dimension(nmwd,nmwd,nr,swd0:swd1), &
    !   intent(in) :: K_mu_wh,K_kp_wh,K_Tro_wh,K_Vro_wh
    ! complex*16, dimension(nmwd,nmwd,nbnd,swd0:swd1), &
    !   intent(in) :: K_Td_wh,K_Vd_wh
    integer, allocatable, dimension(:,:,:), intent(in) :: iwdhs
    complex*16, allocatable, dimension(:,:,:,:), &
      intent(in) :: K_mu_wh,K_kp_wh,K_Tro_wh,K_Vro_wh
    complex*16, allocatable, dimension(:,:,:,:), &
      intent(in) :: K_Td_wh,K_Vd_wh

    ! outputs
    complex*16, dimension(ndim,ndim), intent(out) :: &
               V_mat,W_mat,T_mat,dV_mat

    ! local variables
    character(len=1) :: sym,symp
    integer :: ieig,i,ivec1,ivec2,ivec,l,m,s,t,j,id,i1,i2, &
               jeig,lp,mp,jvec1,jvec2,jvec,ibnd,type,typep,n,np
    integer :: idxi, idxj
    real*8, parameter :: eps = 1.e-7
    real*8 :: om,q,av,ahs,aht,ap,fac,b0,b1,b2,nul,nus,nulp,rho1,rho2,d20
    real*8 :: sumr,sumi,fctr,fcti,fac1,fac2,bom,int1,int2
    real*8, dimension(nr) :: fun,s1,s2,s3,funi,funr
    complex*16 :: ii = cmplx(0.,1.)
    complex*16 :: sumc,jump,WW,VV
    complex*16, dimension(nr,0:sbuild1)   :: K_mu_s
    complex*16, dimension(nr,0:sbuild1)   :: K_kap_s
    complex*16, dimension(nr,0:sbuild1)   :: K_Vrho_s
    complex*16, dimension(nr,0:sbuild1)   :: K_Trho_s
    complex*16, dimension(nbnd,0:sbuild1) :: K_Vd_s
    complex*16, dimension(nbnd,0:sbuild1) :: K_Td_s

    !MD added
    real*8, dimension(nr) :: funrT, funiT
    real*8 :: sumTr, sumTi
    complex*16 :: sumT
    real*8 :: sumdVr, sumdVi, sumdV
    real*8, dimension(nr) :: funrdV, funidV
    complex*16 :: w1 ! Jan31
    !MD added end

    ! initialize the coupling matrices      
    T_mat(:,:)  = dcmplx(0.d0,0.d0)
    W_mat(:,:)  = dcmplx(0.d0,0.d0) 
    V_mat(:,:)  = dcmplx(0.d0,0.d0)
    dV_mat(:,:) = dcmplx(0.d0,0.d0)

    ! set the non-dimensional rotation frequency
    bom = PI2/(24.d0*3600)
    bom = bom/f_norm

    !-----------------------------------------------------!
    !     build up the splitting matrix block by block    !
    !-----------------------------------------------------!

    ! begin loop over the first multiplet index
    ivec2 = 0
    do ieig = 1,nmode
       ! get degree of the multiplet
       l   = eiginfo(ieig)%l
       om  = eiginfo(ieig)%w
       nul = sqrt((2*l+1)/(4.*pi))
       q = eiginfo(ieig)%q    ! Jan31
       !q = q  * 2.17   !!! HL why is this here?
       q = 1/q 
       w1 = dcmplx(om,.5*om*q)
       ! HL - let's keep it non-dimensional 
       ! w1 = w1*f_norm --> Mrina;
       !w1 = dcmplx(om,0.0) ! HLCC elastic
       ! set singlet indices
       ivec1 = ivec2 + 1
       ivec2 = ivec1 + 2*l

       type  = eiginfo(ieig)%itype
       n     = eiginfo(ieig)%n
       if(type == 0) then
          sym = 'T'
       else
          sym = 'S'
       end if

       ! HL: attempt to keep it all non-dimensional
      
       ! add in the spherical reference terms:
       ! for T, this is forcing frequency, so we keep this
       ! for later.
       do m = -l,l
         ivec = ivec1+l+m
         V_mat(ivec,ivec) = V_mat(ivec,ivec) + w1*w1
       enddo

       ! begin loop over the second multiplet index
       jvec2 = 0
       do jeig = 1,nmode
         ! get degree of the multiplet
         lp = eiginfo(jeig)%l
         nulp = sqrt((2*lp+1)/(4.*pi))
         ! set the singlet indices
         jvec1 = jvec2 + 1
         jvec2 = jvec1 + 2*lp

         typep = eiginfo(jeig)%itype
         np    = eiginfo(jeig)%n
         if(typep == 0) then
           symp = 'T'
         else
           symp = 'S'
         end if

         ! print *, '(',n,sym,l,')  x  (',np,symp,lp,')'

         ! ================================================== !
         ! HL (March 2023)
         ! We used to calcualte Woodhouse kernels but now
         ! these are precalculated:
         !
         ! old code:
         ! !get the Woodhouse kernels for the multiplet-pair
           !call woodhouse_kernels(ieig,jeig,sbuild0,sbuild1,&
           !                       K_mu_s,K_kap_s,&
           !                       K_Trho_s,K_Td_s,K_Vrho_s,K_Vd_s) 
         ! new code:
         idxi = iwdhs(type,n,l)
         idxj = iwdhs(typep,np,lp)
         K_mu_s(:,:)   = K_mu_wh(idxi,idxj,1:nr,0:sbuild1)
         K_kap_s(:,:)  = K_kp_wh(idxi,idxj,1:nr,0:sbuild1)
         K_Trho_s(:,:) = K_Tro_wh(idxi,idxj,1:nr,0:sbuild1)
         K_Vrho_s(:,:) = K_Vro_wh(idxi,idxj,1:nr,0:sbuild1)
         K_Td_s(:,:)   = K_Td_wh(idxi,idxj,1:nbnd,0:sbuild1)
         K_Vd_s(:,:)   = K_Vd_wh(idxi,idxj,1:nbnd,0:sbuild1)
         !
         ! ================================================== !


          
        ! begin loop over structural degree
        do s = sbuild0,sbuild1
          nus = sqrt((2*s+1)/(4.*pi))
          fac = 4.*pi*nul*nus*nulp 
          ! begin loop over structural orders
          do t = -s,s
            !--------------------------------------------------------!
            !   perform the radial integrals summing contributions   !
            !--------------------------------------------------------!

             ! For T (D.43)
            funrT(:) = real(rho_st(:,s,t)*K_Trho_s(:,s))
            funiT(:) = aimag(rho_st(:,s,t)*K_Trho_s(:,s))
            call intgrl(sumTr,r,1,nr,funrT,s1,s2,s3)
            call intgrl(sumTi,r,1,nr,funiT,s1,s2,s3)
            sumT = sumTr + ii*sumTi 

            ! for V - rho terms
            funr(:) =  real(rho_st(:,s,t)*K_Vrho_s(:,s))
            funi(:) = aimag(rho_st(:,s,t)*K_Vrho_s(:,s))
            call intgrl(sumr,r,1,nr,funr,s1,s2,s3)
            call intgrl(sumi,r,1,nr,funi,s1,s2,s3)
            sumc = sumr + ii*sumi
            !MD added end

            ! mu terms -- are only for V
            funr(:) =  real(mu_st(:,s,t)*K_mu_s(:,s))
            funi(:) = aimag(mu_st(:,s,t)*K_mu_s(:,s))
            call intgrl(sumr,r,1,nr,funr,s1,s2,s3)
            call intgrl(sumi,r,1,nr,funi,s1,s2,s3)
            sumc = sumc + sumr + ii*sumi

            ! kappa terms -- are only for V
            funr(:) =  real(kap_st(:,s,t)*K_kap_s(:,s))
            funi(:) = aimag(kap_st(:,s,t)*K_kap_s(:,s))
            call intgrl(sumr,r,1,nr,funr,s1,s2,s3)
            call intgrl(sumi,r,1,nr,funi,s1,s2,s3)
            sumc = sumc + sumr + ii*sumi


            !--------------------------------------------------------!
            !   delta perturbation to V -- kappa and mu terms        !
            !--------------------------------------------------------!                

            ! Delta pertubations to V -- kappa terms
            funrdV(:) = real(dkap_st(:,s,t)*K_kap_s(:,s))
            funidV(:) = aimag(dkap_st(:,s,t)*K_kap_s(:,s))
            call intgrl(sumdVr,r,1,nr,funrdV,s1,s2,s3)
            call intgrl(sumdVi,r,1,nr,funidV,s1,s2,s3)
            sumdV = sumdVr + ii*sumdVi

            ! Delta pertubations to V -- mu terms
            funrdV(:) = real(dmu_st(:,s,t)*K_mu_s(:,s))
            funidV(:) = aimag(dmu_st(:,s,t)*K_mu_s(:,s))
            call intgrl(sumdVr,r,1,nr,funrdV,s1,s2,s3)
            call intgrl(sumdVi,r,1,nr,funidV,s1,s2,s3)
            sumdV = sumdV + sumdVr + ii*sumdVi

            !--------------------------------------------------!
            !             deal with boundary terms             !
            !--------------------------------------------------!

            sumT = sumT + sum(K_Td_s(:,s)*d_st(:,s,t))
            sumc = sumc + sum(K_Vd_s(:,s)*d_st(:,s,t))

            ! begin loop over singlets for first multiplet
            do m = -l,l
              ivec = ivec1+l+m
              fac1 = (-1)**m*fac

              ! begin loop over singlets for the second multiplet
              do mp = -lp,lp
                jvec = jvec1+lp+mp
                ! get the geometric terms
                fac2 = fac1*thrj(l,s,lp,-m,t,mp)

                T_mat(ivec,jvec) = T_mat(ivec,jvec) + fac2*sumT
                V_mat(ivec,jvec) = V_mat(ivec,jvec) + fac2*sumc

                ! for delta matrix
                dV_mat(ivec,jvec) = dV_mat(ivec,jvec) + fac2*sumdV
              end do
              ! end loop over singlets for the second multiplet
            end do
            ! end loop over singlets for first multiplet                
          end do
          ! end loop over structural orders
        end do
        ! end loop over structural degrees
 
        !------------------------------------------------------!
        !                 add in rotational terms              !
        !------------------------------------------------------!
        if(is_rotate) then

        ! perform the necessary integrals
          call rotation_integrals(ieig,jeig,int1,int2)

          ! begin loop over singlets for first multiplet
          do m = -l,l
            ivec = ivec1+l+m

            ! begin loop over singlets for second multiplet
            do mp = -lp,lp
              jvec = jvec1+lp+mp

               if(type == typep) then
                 ! like coupling
                 ! D.68
                 WW = m*bom*delta(l,lp)*delta(m,mp)*int1
                 ! D.72 
                 VV = -(2./3.)*l*(l+1)*bom**2*delta(l,lp)*delta(m,mp)*int1&
                      +(-1)**m*sqrt((2*l+1.)*(2*lp+1.))*bom**2&
                      *thrj(l,2,lp,-m,0,mp)*int2
                 if(n == np) then
                   VV = VV + (2./3.)*bom**2*delta(l,lp)*delta(m,mp)
                 end if
               else
                 ! unlike coupling
                 ! D.68
                 if(l == lp+1) then
                   WW = -ii*bom*slm(l,m)*delta(m,mp)*int1
                 else if(l == lp-1) then
                   WW = -ii*bom*slm(lp,m)*delta(m,mp)*int1
                 else
                   WW = 0.
                 end if
                 ! D.72
                   VV = (-1)**m*sqrt((2*l+1.)*(2*lp+1.))*bom**2 &
                         *thrj(l,2,lp,-m,0,mp)*ii*int2
               end if
               V_mat(ivec,jvec) = V_mat(ivec,jvec) + VV 
               W_mat(ivec,jvec) = W_mat(ivec,jvec) + WW              
             end do
             ! end loop over singlets for second multiplet
           end do
           ! end loop over singlets for first multiplet
         end if
       end do
       ! end loop over the second multiplet index
     end do
    ! add loop over the first multiplet index

    return

  END SUBROUTINE coupling_matrix
 
  ! ================================================================ !
  ! Decomissioned 
  ! SUBROUTINE woodhouse_kernels(ieig,jeig,sbuild0,sbuild1,&
  !                              K_mu_s,K_kap_s,&
  !                              K_Trho_s,K_Td_s,K_Vrho_s,K_Vd_s)
 
  !   ! this routine returns the Woodhouse kernels for the 
  !   ! (ieig,jeig)th multiplet pair.
  !   ! The results can be used either in computing matrix 
  !   ! elements or in calculating the associated sensitivity kernels. 
  !   !
  !   ! inputs:
  !   !
  !   ! ieig  -- multiplet index for the first mode
  !   ! jeig  -- multiplet index for the second mode
  !   ! sbuild-- maximum structural degree to consider
  !   ! idisc -- indices for the major discontinuitities
  !   !
  !   ! outputs:
  !   ! frequency independent: HL
  !   ! K_mu_s    -- kernel for shear modulus at each radius
  !   ! K_kap_s   -- kernel for bulk modulus at each radius
  !   ! frequency dependent: HL
  !   ! K_Trho_s  -- kernel for kinetic density at each radius 
  !   ! K_Vrho_s  -- kernel for potential density at each radius 
  !   ! K_Td_s    -- kernel for boundary topography at each 
  !   !              major discontinuity
  !   !
  !   !
  !   ! Note that the kernels are defined for ABSOLUTE and not 
  !   ! relative perturbations

  !   ! MD added
  !   ! Seperating the kernels for T_rho, V_rho, T_d, V_d
  !   ! refer to D43, D44 to see what is meant
  !   ! MD added end
    
  
  !   USE modellayer_module, ONLY:NR,NCMB,NMOHO,kdis
  !   USE prem_module
  !   USE mantle_model_module_DS
  !   USE futil_module, ONLY : intgrl
  !   USE eigfun_module, ONLY : r=>eigfun_radius, &
  !                                eigfun_u, eigfun_du, &
  !                                eigfun_v, eigfun_dv, &
  !                                eigfun_w, eigfun_dw, &
  !                                eigfun_p, eigfun_dp, &
  !                                eiginfo
  
  !   implicit none
  !   include "coupling.h"

  !   ! inputs/outputs
  !   integer, intent(in) :: ieig
  !   integer, intent(in) :: jeig
  !   integer, intent(in) :: sbuild0,sbuild1
  !   complex*16, dimension(nr,0:sbuild1), intent(out) :: K_mu_s
  !   complex*16, dimension(nr,0:sbuild1), intent(out) :: K_kap_s
  !   complex*16, dimension(nr,0:sbuild1), intent(out) :: K_Trho_s ! T_rho
  !   complex*16, dimension(nr,0:sbuild1), intent(out) :: K_Vrho_s ! V_rho
  !   complex*16, dimension(nbnd,0:sbuild1), intent(out) :: K_Td_s ! T_d
  !   complex*16, dimension(nbnd,0:sbuild1), intent(out) :: K_Vd_s ! V_d
    
  !   ! local variables  
  !   integer :: np,typep,lp, n,type,l
  !   integer :: i,k,s,d,ds,ibnd
  !   real*8    :: omp,qp,gvp,avp,ahsp,ahtp
  !   real*8    :: om,q,gv,av,ahs,aht
  !   real*8    :: blp, bl
  !   real*8    :: fp,xp,zp,f,x,z
  !   real*8    :: b1,b2,b3,b4,b5,b6,b7,b8,b9
  !   real*8    :: gs1_r,gs2_r,i_1_r,i_2_r
  !   real*8    :: gs1_i,gs2_i,i_1_i,i_2_i
  !   real*8    :: ks_r,ms_r,kt1s_r,rs1_r,rs2_r
  !   real*8    :: kt2s_r,kt3s_r,kt4s_r,kt5s_r
  !   real*8    :: ms_i,rs1_i,rs2_i,kt1s_i,kt2s_i,kt4s_i,kt5s_i
  !   real*8    :: ts_r,ts_i
  !   real*8    :: kappa0,mu0
  !   real*8    :: c11,c12,c13,c14,c15,kappa1,mu1
  !   real*8    :: b2s_r0,s2s_r0,s2s_i0,rs_r0,rs_i0
  !   real*8    :: b2s_r1,s2s_r1,s2s_i1,rs_r1,rs_i1
  !   complex*16 :: ii = cmplx(0.,1.)  

  !   ! These are of dimension NR:
  !   real*8, pointer :: up(:),dup(:),vp(:),dvp(:),wp(:),dwp(:),pp(:),dpp(:)
  !   real*8, pointer :: u(:), du(:), v(:), dv(:), w(:), dw(:), p(:), dp(:)
  !   real*8, dimension(:), allocatable :: k1_r,k2_r,k1_i,k2_i
  !   complex*16, dimension(:), allocatable :: jump,jumpt
  !   real*8, dimension(:), allocatable :: s1,s2,s3
  !   ! These are of dimsnion 0:NS x NR
  !   real*8, dimension(:),allocatable :: substitution_r,substitution_i
  !   real*8 qinvks_A,qinvms_A
  !   integer :: itmp,ndis
      
  !   IF (.NOT.(ALLOCATED(r).AND.ALLOCATED(eigfun_u)&
  !        .AND.ALLOCATED(eiginfo))) STOP &
  !        'Variables in eigfun_module has to be defined previously'
  
  !   NULLIFY(u, du, v, dv, w, dw, p, dp)
  !   NULLIFY(up,dup,vp,dvp,wp,dwp,pp,dpp)
    
  !   ! define eigenfunctions for mode 1 (p)
  !   om=eiginfo(ieig)%w
  !   q=eiginfo(ieig)%q 
  !   n=eiginfo(ieig)%n
  !   l=eiginfo(ieig)%l
  !   type=eiginfo(ieig)%itype
  !   u=>eigfun_u(:,ieig)
  !   v=>eigfun_v(:,ieig)
  !   w=>eigfun_w(:,ieig)
  !   p=>eigfun_p(:,ieig)
  !   du=>eigfun_du(:,ieig)
  !   dv=>eigfun_dv(:,ieig)
  !   dw=>eigfun_dw(:,ieig)
  !   dp=>eigfun_dp(:,ieig)
  !   bl=l*(l+1.0)
    
  !   ! define eigenfunctions for mode 2
  !   omp=eiginfo(jeig)%w
  !   qp=eiginfo(jeig)%q 
  !   np=eiginfo(jeig)%n
  !   lp=eiginfo(jeig)%l
  !   typep=eiginfo(jeig)%itype
  !   up=>eigfun_u(:,jeig)
  !   vp=>eigfun_v(:,jeig)
  !   wp=>eigfun_w(:,jeig)
  !   pp=>eigfun_p(:,jeig)
  !   dup=>eigfun_du(:,jeig)
  !   dvp=>eigfun_dv(:,jeig)
  !   dwp=>eigfun_dw(:,jeig)
  !   dpp=>eigfun_dp(:,jeig)
  !   blp=lp*(lp+1.0)
    
  !   IF (SIZE(r)/=NR .OR. SIZE(u)/=NR .OR. SIZE(w)/=NR .OR. &
  !        SIZE(up)/=NR.OR. SIZE(wp)/=NR) &
  !        STOP 'dimension of eigenfunctions has to be equal, in SUB &
  !        initialize_mantle STOP'
  !   ALLOCATE(k1_r(NR),k2_r(NR),k1_i(NR),k2_i(NR),&
  !            jump(NR),jumpt(NR),s1(NR),s2(NR),s3(NR),STAT=itmp)
  !   IF (itmp/=0) STOP 'allocate error 0 in mantle.f, STOP '
  !   k1_r=0.0
  !   k2_r=0.0
  !   k1_i=0.0
  !   k2_i=0.0      
    
  !   ALLOCATE(substitution_r(NR),substitution_i(NR), &
  !        STAT=itmp)
  !   IF (itmp/=0) STOP 'allocate error 2 in mantle.f, STOP'
  !   substitution_r=0.0
  !   substitution_i=0.0    

  !   ! initialise the kernels
  !   K_mu_s        = 0.
  !   K_kap_s       = 0.
  !   K_Trho_s(:,:) = 0.
  !   K_Vrho_s(:,:) = 0. 
  !   K_Td_s(:,:)   = 0.
  !   K_Vd_s(:,:)   = 0.

  !   if(type.eq.typep) then ! like-type coupling
  !     do s=sbuild0,sbuild1
  !       ! check the selection rules
  !       if(l < abs(lp-s) .or. l > lp+s) cycle
  !       if(mod(l+s+lp,2) /= 0) cycle
  !       b1=b(l,s,lp,0,1)
  !       b2=b(l,s,lp,2,1)
  !       b3=b(l,s,lp,1,1)
  !       b6=b(lp,l,s,1,1)
  !       b7=b(l,lp,s,1,1)
  !       do i=2,NR
  !         fp=(2.0*up(i)-blp*vp(i))/r(i)
  !         f=(2.0*u(i)-bl*v(i))/r(i)
  !         !\ REAL part of V{phi} in (D.51) due to like-type coupling
  !         gs1_r=0.5*rho(i)*(u(i)*dvp(i)+u(i)*vp(i)/r(i) &
  !               -du(i)*vp(i)-2.0*f*vp(i))*b6/r(i) &
  !               +0.5*rho(i)*(up(i)*dv(i)+up(i)*v(i)/r(i) &
  !               -dup(i)*v(i)-2.0*fp*v(i))*b7/r(i) &
  !               +rho(i)*u(i)*up(i)*s*(s+1.0)*b1/(r(i)*r(i))
  !         !\ REAL part of V{dphi/dr} in (D.52) due to like-type
  !         !  coupling
  !         gs2_r=0.5*rho(i)*u(i)*vp(i)*b6/r(i)+0.5*rho(i)*up(i)*v(i)*b7/r(i) &
  !               -rho(i)*(fp*u(i)+up(i)*f)*b1
  !         !\ the 1st and second integrand/(r^2) of V{rho}_substitution in (D.55)
  !         k1_r(i)=(r(i)**(-s))*((s+1.0)*gs2_r-r(i)*gs1_r)
  !         k2_r(i)=(r(i)**(s+1.0))*(s*gs2_r+r(i)*gs1_r)
  !       enddo
  !       k1_r = k1_r/(r*r)
  !       k2_r = k2_r/(r*r)
  !       k1_r(1)=k1_r(2)
  !       k2_r(1)=k2_r(2)
  !       do i=1,NR
  !         call intgrl(i_1_r,r,i,NR,k1_r,s1,s2,s3)
  !         call intgrl(i_2_r,r,1,i,k2_r,s1,s2,s3)
  !         !\ the 2nd term of V{rho}_substitution in (D.55)
  !         substitution_r(i)= (4.0/(2.0*s+1.0))*((r(i)**s)*i_1_r &
  !                            -(r(i)**(-s-1.0))*i_2_r)
  !       enddo
         
  !       do i=2,NR
  !         fp=(2.0*up(i)-blp*vp(i))/r(i)
  !         xp=dvp(i)+(up(i)-vp(i))/r(i)
  !         zp=dwp(i)-wp(i)/r(i)
  !         f=(2.0*u(i)-bl*v(i))/r(i)
  !         x=dv(i)+(u(i)-v(i))/r(i)
  !         z=dw(i)-w(i)/r(i)
  !         !/ V{kappa} in (D.48)
  !         ks_r=(dup(i)+fp)*(du(i)+f)*b1
  !         !/ REAL part of V{mu} in (D.49)
  !         ms_r=(vp(i)*v(i)+wp(i)*w(i))*b2/(r(i)*r(i))+(x*xp+z*zp)*b3 &
  !              +(1.0/3.0)*(2.0*dup(i)-fp)*(2.0*du(i)-f)*b1
  !         !/ rs1_r is the REAL (like-type coupling) part of V{rho} in (D.50)
  !         rs1_r=((pp(i)*v(i)+p(i)*vp(i))/r(i) &
  !               +0.5*g(i)*(up(i)*v(i)+vp(i)*u(i))/r(i))*b3 &
  !               +(8.0*rho(i)*u(i)*up(i)+dpp(i)*u(i)+dp(i)*up(i) &
  !               -0.5*g(i)*(4.0*u(i)*up(i)/r(i)+up(i)*f+u(i)*fp))*b1
  !         !/ like-type coupling part of V{rho}_substitute in (D.55)
  !         rs2_r=rs1_r+substitution_r(i)
  !         !/ REAL part of T{rho} in *D.46)
  !         ts_r=u(i)*up(i)*b1+(vp(i)*v(i)+wp(i)*w(i))*b3
  !         !/ REAL part of V{d} - refer to Henson, GJI, 1989
  !         !      boundary perturbations of TI earth model:
  !         !      C, 2N, A-N, -F, -L
  !         kt1s_r=-dup(i)*du(i)*b1+du(i)*vp(i)*b6/r(i)+dup(i)*v(i)*b7/r(i)
  !         kt2s_r=0.5*(vp(i)*v(i)+wp(i)*w(i))*b2/(r(i)*r(i))
  !         kt3s_r=fp*f*b1
  !         kt4s_r=-f*vp(i)*b6/r(i)-fp*v(i)*b7/r(i)
  !         kt5s_r=(-xp*x-zp*z+dvp(i)*x+dv(i)*xp+dwp(i)*z+dw(i)*zp)*b3
            
  !         ! d^2*V{d}
  !         jump(i)  = -( c1(i)*kt1s_r+c2(i)*kt2s_r+c3(i)*kt3s_r &
  !                    +c4(i)*kt4s_r+c5(i)*kt5s_r+rho(i)*rs1_r)*r(i)**2
  !         !/ d^2 * T{d} in (D.47)
  !         jumpt(i) = -rho(i)*ts_r*r(i)*r(i)

  !         !-------------------------------------------------!
  !         !    store the combined kernels at each radius    !
  !         !-------------------------------------------------!

  !         ! kappa kernel
  !         K_kap_s(i,s) = ks_r
  !         ! mu kernel
  !         K_mu_s(i,s) = ms_r
  !         ! rho kernels with frequency extracted:
  !         K_Trho_s(i,s) = -ts_r !F1
  !         K_Vrho_s(i,s) = rs2_r
  !         ! boundary kernels
  !         do ibnd=1,nbnd-1
  !           if(i.eq.bnd(ibnd)+1) then
  !             ! boundary kernels with frequency extracted:
  !             K_Td_s(ibnd,s) = -(jumpt(i)-jumpt(i-1)) !F1
  !             K_Vd_s(ibnd,s) = (jump(i)-jump(i-1))
  !           endif
  !         enddo
  !         if(i == NR) then
  !           ! boundary kernels with frequency extracted:
  !           K_Td_s(nbnd,s) = jumpt(i) 
  !           K_Vd_s(nbnd,s) = - jump(i)
  !         end if
  !       enddo
          
  !     enddo

  !   else ! unlike-type coupling
  !      do s=sbuild0,sbuild1
  !        ! check the selection rules
  !        if(l < abs(lp-s) .or. l > lp+s) cycle
  !        if(mod(l+s+lp,2) == 0) cycle
          
  !        b4=b(l,s,lp,2,-1)
  !        b5=b(l,s,lp,1,-1)
  !        b8=b(l,lp,s,1,-1)
  !        b9=b(lp,l,s,1,-1)
          
  !        do i=2,NR
  !          fp=(2.0*up(i)-blp*vp(i))/r(i)
  !          f=(2.0*u(i)-bl*v(i))/r(i)
  !          gs1_i=0.5*rho(i)*(u(i)*dwp(i)+u(i)*wp(i)/r(i) &
  !               -du(i)*wp(i)-2.0*f*wp(i))*b9/r(i) &
  !               -0.5*rho(i)*(up(i)*dw(i)+up(i)*w(i)/r(i) &
  !               -dup(i)*w(i)-2.0*fp*w(i))*b8/r(i)
  !          gs2_i=0.5*rho(i)*u(i)*wp(i)*b9/r(i)-0.5*rho(i)*up(i)*w(i)*b8/r(i)
  !          k1_i(i)=(r(i)**(-s-2.0))*((s+1.0)*gs2_i-r(i)*gs1_i)
  !          k2_i(i)=(r(i)**(s-1.0))*(s*gs2_i+r(i)*gs1_i)
  !        enddo
  !        k1_i(1)=k1_i(2)
  !        k2_i(1)=k2_i(2)
  !        do i=1,NR
  !          call intgrl(i_1_i,r,i,NR,k1_i,s1,s2,s3)
  !          call intgrl(i_2_i,r,1,i,k2_i,s1,s2,s3)
  !          substitution_i(i)= (4.0/(2.0*s+1.0))*((r(i)**s)*i_1_i &
  !                             -(r(i)**(-s-1.0))*i_2_i)
  !        enddo
         
  !        do i=2,NR
  !          fp=(2.0*up(i)-blp*vp(i))/r(i)
  !          xp=dvp(i)+(up(i)-vp(i))/r(i)
  !          zp=dwp(i)-wp(i)/r(i)
  !          f=(2.0*u(i)-bl*v(i))/r(i)
  !          x=dv(i)+(u(i)-v(i))/r(i)
  !          z=dw(i)-w(i)/r(i)
  !          ! imaginary part of DT.49
  !          ms_i=(vp(i)*w(i)-wp(i)*v(i))*b4/(r(i)*r(i))+(xp*z-x*zp)*b5
  !          ! imaginary part of DT.50
  !          rs1_i=((pp(i)*w(i)-p(i)*wp(i))/r(i) &
  !                +0.5*g(i)*(up(i)*w(i)-wp(i)*u(i))/r(i))*b5
  !          ! imaginary part of DT.46
  !          ts_i=(vp(i)*w(i)-wp(i)*v(i))*b5
  !          ! imaginary part of DT.55
  !          rs2_i=rs1_i+substitution_i(i)
  !          ! imaginary parts of boundary kernels
  !          kt1s_i=du(i)*wp(i)*b8/r(i)-dup(i)*w(i)*b9/r(i)
  !          kt2s_i=0.5*(vp(i)*w(i)-wp(i)*v(i))*b4/(r(i)*r(i))
  !          kt4s_i=-f*wp(i)*b8/r(i)+fp*w(i)*b9/r(i)
  !          kt5s_i=(-xp*z+zp*x+dvp(i)*z-dv(i)*zp+dw(i)*xp-dwp(i)*x)*b5
            
  !          ! d^2*V{d} unrelated to physical dispersion
  !          jump(i) = -( c1(i)*kt1s_i+c2(i)*kt2s_i &
  !                    +c3(i)*kt4s_i+c5(i)*kt5s_i+rho(i)*rs1_i)*r(i)**2
  !          ! r^2*T{d} in (D.47)
  !          jumpt(i) = -r(i)*r(i)*rho(i)*ts_i

  !          !-------------------------------------------------!
  !          !    store the combined kernels at each radius    !
  !          !-------------------------------------------------!
             
  !          ! mu kernel
  !          K_mu_s(i,s) = ii*ms_i

  !          ! rho kernels with frequency extracted:
  !          K_Trho_s(i,s) = -ts_i !F1
  !          K_Trho_s(i,s) = ii*K_Trho_s(i,s)
  !          K_Vrho_s(i,s) = rs2_i
  !          K_Vrho_s(i,s) = ii*K_Vrho_s(i,s)

  !          ! boundary kernels
  !          do ibnd=1,nbnd-1
  !            if(i.eq.bnd(ibnd)+1) then
  !             ! boundary kernels with frequency extracted:
  !              K_Td_s(ibnd,s) = -(jumpt(i)-jumpt(i-1)) !F1
  !              K_Td_s(ibnd,s) = ii*K_Td_s(ibnd,s)
  !              K_Vd_s(ibnd,s) = (jump(i)-jump(i-1))
  !              K_Vd_s(ibnd,s) = ii*K_Vd_s(ibnd,s)
  !            endif
  !          enddo
  !          if(i == NR) then
  !           ! boundary kernels with frequency extracted:
  !            K_Td_s(nbnd,s) = jumpt(i)
  !            K_Td_s(ibnd,s) = ii*K_Td_s(ibnd,s)
  !            K_Vd_s(nbnd,s) = -jump(i)
  !            K_Vd_s(ibnd,s) = ii*K_Vd_s(ibnd,s)
  !          end if

  !        enddo
          
  !      enddo
       
  !    endif
 
  !   DEALLOCATE(substitution_r,substitution_i)
  !   DEALLOCATE(k1_r,k2_r,k1_i,k2_i,jump,jumpt,s1,s2,s3)
    
  !   NULLIFY(u, du, v, dv, w, dw, p, dp)
  !   NULLIFY(up,dup,vp,dvp,wp,dwp,pp,dpp)
  !   return

  ! END SUBROUTINE woodhouse_kernels
 
  ! ================================================================ !
  ! new version
  SUBROUTINE woodhouse_kernels(ieig,jeig,sbuild0,sbuild1,&
                               K_mu_s,K_kap_s,&
                               K_Trho_s,K_Td_s,K_Vrho_s,K_Vd_s)

    ! this routine returns the Woodhouse kernels for the
    ! (ieig,jeig)th multiplet pair.
    ! The results can be used either in computing matrix
    ! elements or in calculating the associated sensitivity
    ! kernels.
    !
    ! inputs:
    !
    ! ieig  -- multiplet index for the first mode
    ! jeig  -- multiplet index for the second mode
    ! sbuild-- maximum structural degree to consider
    ! idisc -- indices for the major discontinuitities
    !
    ! outputs:
    ! frequency independent: HL
    ! K_mu_s    -- kernel for shear modulus at each radius
    ! K_kap_s   -- kernel for bulk modulus at each radius
    ! frequency dependent: HL
    ! K_Trho_s  -- kernel for kinetic density at each radius
    ! K_Vrho_s  -- kernel for potential density at each radius
    ! K_Td_s    -- kernel for boundary topography at each
    !              major discontinuity
    !
    !
    ! Note that the kernels are defined for ABSOLUTE and not
    ! relative perturbations

    ! MD added
    ! Seperating the kernels for T_rho, V_rho, T_d, V_d
    ! refer to D43, D44 to see what is meant
    ! MD added end

    USE modellayer_module, ONLY:NR,NCMB,NMOHO,kdis
    USE prem_module
    USE mantle_model_module_DS
    USE futil_module, ONLY : intgrl
    USE eigfun_module, ONLY : r=>eigfun_radius, &
                                 eigfun_u, eigfun_du, &
                                 eigfun_v, eigfun_dv, &
                                 eigfun_w, eigfun_dw, &
                                 eigfun_p, eigfun_dp, &
                                 eiginfo

    implicit none
    include "coupling.h"

    ! inputs/outputs
    integer, intent(in) :: ieig
    integer, intent(in) :: jeig
    integer, intent(in) :: sbuild0,sbuild1
    complex*16, dimension(nr,0:sbuild1), intent(out) :: K_mu_s
    complex*16, dimension(nr,0:sbuild1), intent(out) :: K_kap_s
    complex*16, dimension(nr,0:sbuild1), intent(out) :: K_Trho_s ! T_rho
    complex*16, dimension(nr,0:sbuild1), intent(out) :: K_Vrho_s ! V_rho
    complex*16, dimension(nbnd,0:sbuild1), intent(out) :: K_Td_s ! T_d
    complex*16, dimension(nbnd,0:sbuild1), intent(out) :: K_Vd_s ! V_d

    ! local variables
    integer :: np,typep,lp, n,type,l
    integer :: i,k,s,d,ds,ibnd
    real*8    :: omp,qp,gvp,avp,ahsp,ahtp
    real*8    :: om,q,gv,av,ahs,aht
    real*8    :: blp, bl
    real*8    :: fp,xp,zp,f,x,z
    real*8    :: b1,b2,b3,b4,b5,b6,b7,b8,b9
    real*8    :: gs1_r,gs2_r,i_1_r,i_2_r
    real*8    :: gs1_i,gs2_i,i_1_i,i_2_i
    real*8    :: ks_r,ms_r,kt1s_r,rs1_r,rs2_r
    real*8    :: kt2s_r,kt3s_r,kt4s_r,kt5s_r
    real*8    :: ms_i,rs1_i,rs2_i,kt1s_i,kt2s_i,kt4s_i,kt5s_i
    real*8    :: ts_r,ts_i
    real*8    :: kappa0,mu0
    real*8    :: c11,c12,c13,c14,c15,kappa1,mu1
    real*8    :: b2s_r0,s2s_r0,s2s_i0,rs_r0,rs_i0
    real*8    :: b2s_r1,s2s_r1,s2s_i1,rs_r1,rs_i1
    complex*16 :: ii = cmplx(0.,1.)

    ! These are of dimension NR:
    real*8, pointer :: up(:),dup(:),vp(:),dvp(:),wp(:),dwp(:),pp(:),dpp(:)
    real*8, pointer :: u(:), du(:), v(:), dv(:), w(:), dw(:), p(:), dp(:)
    real*8, dimension(:), allocatable :: k1_r,k2_r,k1_i,k2_i
    complex*16, dimension(:), allocatable :: jump,jumpt
    real*8, dimension(:), allocatable :: s1,s2,s3
    ! These are of dimsnion 0:NS x NR
    real*8, dimension(:),allocatable :: sub_r,sub_i
    real*8 qinvks_A,qinvms_A
    integer :: itmp,ndis
      
    IF (.NOT.(ALLOCATED(r).AND.ALLOCATED(eigfun_u)&
         .AND.ALLOCATED(eiginfo))) STOP &
         'Variables in eigfun_module has to be defined previously'

    NULLIFY(u, du, v, dv, w, dw, p, dp)
    NULLIFY(up,dup,vp,dvp,wp,dwp,pp,dpp)

    ! define eigenfunctions for mode 1 (p)
    om=eiginfo(ieig)%w
    q=eiginfo(ieig)%q
    n=eiginfo(ieig)%n
    l=eiginfo(ieig)%l
    type=eiginfo(ieig)%itype
    u=>eigfun_u(:,ieig)
    v=>eigfun_v(:,ieig)
    w=>eigfun_w(:,ieig)
    p=>eigfun_p(:,ieig)
    du=>eigfun_du(:,ieig)
    dv=>eigfun_dv(:,ieig)
    dw=>eigfun_dw(:,ieig)
    dp=>eigfun_dp(:,ieig)
    bl=l*(l+1.0)

    ! define eigenfunctions for mode 2
    omp=eiginfo(jeig)%w
    qp=eiginfo(jeig)%q
    np=eiginfo(jeig)%n
    lp=eiginfo(jeig)%l
    typep=eiginfo(jeig)%itype
    up=>eigfun_u(:,jeig)
    vp=>eigfun_v(:,jeig)
    wp=>eigfun_w(:,jeig)
    pp=>eigfun_p(:,jeig)
    dup=>eigfun_du(:,jeig)
    dvp=>eigfun_dv(:,jeig)
    dwp=>eigfun_dw(:,jeig)
    dpp=>eigfun_dp(:,jeig)
    blp=lp*(lp+1.0)

    IF (SIZE(r)/=NR .OR. SIZE(u)/=NR .OR. SIZE(w)/=NR .OR. &
         SIZE(up)/=NR.OR. SIZE(wp)/=NR) &
         STOP 'dimension of eigenfunctions has to be equal, in SUB &
         initialize_mantle STOP'
    ALLOCATE(k1_r(NR),k2_r(NR),k1_i(NR),k2_i(NR),&
             jump(NR),jumpt(NR),s1(NR),s2(NR),s3(NR),STAT=itmp)
    IF (itmp/=0) STOP 'allocate error 0 in mantle.f, STOP '
    k1_r=0.0
    k2_r=0.0
    k1_i=0.0
    k2_i=0.0

    ALLOCATE(sub_r(NR),sub_i(NR),STAT=itmp)
    IF (itmp/=0) STOP 'allocate error 2 in mantle.f, STOP'
    sub_r=0.0
    sub_i=0.0

    ! initialise the kernels
    K_mu_s        = 0.
    K_kap_s       = 0.
    K_Trho_s(:,:) = 0.
    K_Vrho_s(:,:) = 0.
    K_Td_s(:,:)   = 0.
    K_Vd_s(:,:)   = 0.

    ! ==================================================== !
    ! BEGIN CALCULCULATIONS
    ! ==================================================== !
    
    ! LIKE-TYPE COUPLING
    if (type.eq.typep) then
      do s=sbuild0,sbuild1
        ! selection rules
        if (l<abs(lp-s).or.l>lp+s) cycle
        if (mod(l+s+lp,2).ne.0) cycle
        b1=b(l,s,lp,0,1) ! B_{l,s,l'}^{0+}
        b2=b(l,s,lp,2,1) ! B_{l,s,l'}^{2+}
        b3=b(l,s,lp,1,1) ! B_{l,s,l'}^{1+}
        b6=b(lp,l,s,1,1) ! B_{l',l,s}^{1+}
        b7=b(l,lp,s,1,1) ! B_{l,l',s}^{1+}
        ! deal with integrals in potential term:
        do i=2,NR
          f  = (2.*u(i)-bl*v(i))/r(i)
          fp = (2.*up(i)-blp*vp(i))/r(i)
          ! Real part of D.51 (like coupling):
          gs1_r = dble(s*(s+1))*rho(i)*u(i)*up(i)*b1/(r(i)*r(i)) &
                + 0.5*rho(i)*(u(i)*dvp(i)-du(i)*vp(i) + &
                  u(i)*vp(i)/r(i)-2.*f*vp(i))*b6/r(i) &
                + 0.5*rho(i)*(dv(i)*up(i)-v(i)*dup(i) + &
                  v(i)*up(i)/r(i)-2.*v(i)*fp)*b7/r(i)
          ! Real part of D.52 (like coupling):
          gs2_r = -rho(i)*(f*up(i)+u(i)*fp)*b1 &
                + 0.5*rho(i)*u(i)*vp(i)*b6/r(i) &
                + 0.5*rho(i)*v(i)*up(i)*b7/r(i)
          ! First and Second integrand/r^2 of D.55:
          k1_r(i) = r(i)**(-s-2)*(dble(s+1)*gs2_r-r(i)*gs1_r)
          k2_r(i) = r(i)**(s-1)*(dble(s)*gs2_r+r(i)*gs1_r)
        enddo
        k1_r(1) = k1_r(2)
        k2_r(1) = k2_r(2)
        ! calculate the substituted phi terms in D.55
        do i=1,NR
          call intgrl(i_1_r,r,i,NR,k1_r,s1,s2,s3)
          call intgrl(i_2_r,r,1,i,k2_r,s1,s2,s3)
          sub_r(i) = i_1_r*r(i)**s - i_2_r*r(i)**(-s-1)
          sub_r(i) = sub_r(i)*4.d0/dble(2*s+1)
        enddo

        ! now compile the kernels
        do i=2,NR
          f  = (2.*u(i)-bl*v(i))/r(i)
          fp = (2.*up(i)-blp*vp(i))/r(i)
          x  = dv(i)-v(i)/r(i)+u(i)/r(i)
          xp = dvp(i)-vp(i)/r(i)+up(i)/r(i)
          z  = dw(i) - w(i)/r(i)
          zp = dwp(i) - wp(i)/r(i)

          ! D.48:
          ks_r  = (du(i)+f)*(dup(i)+fp)*b1
          ! Real part of D.49:
          ms_r  = (2.*du(i)-f)*(2.*dup(i)-fp)*b1/3.d0 &
                + (x*xp+z+zp)*b3 &
                + (v(i)*vp(i)+w(i)*wp(i))*b2/(r(i)*r(i))
          ! Real part of D.50:
          rs1_r = (u(i)*dpp(i)+dp(i)*up(i) - &
                   0.5*g(i)*(4.*u(i)*up(i)/r(i)+f*up(i)+u(i)*fp) + &
                   8.0*rho(i)*u(i)*up(i))*b1 &
                + ((p(i)*vp(i)+v(i)*pp(i)) + &
                   0.5*g(i)*(u(i)*vp(i)+v(i)*up(i)))*b3/r(i)
          ! D.55:
          rs2_r = rs1_r + sub_r(i)
          ! Real part of D.46:
          ts_r  = u(i)*up(i)*b1 + (v(i)*vp(i)+w(i)*wp(i))*b3
          ! Real part of Vd (note this is copied from original
          ! HY code to follow Henson, GJI, 1989 for a TI Earth):
          kt1s_r= -dup(i)*du(i)*b1+du(i)*vp(i)*b6/r(i)+ &
                  dup(i)*v(i)*b7/r(i)
          kt2s_r=0.5*(vp(i)*v(i)+wp(i)*w(i))*b2/(r(i)*r(i))
          kt3s_r=fp*f*b1
          kt4s_r=-f*vp(i)*b6/r(i)-fp*v(i)*b7/r(i)
          kt5s_r=(-xp*x-zp*z+dvp(i)*x+dv(i)*xp+dwp(i)*z+dw(i)*zp)*b3
          ! d^2*V{d}
          jump(i)= -( c1(i)*kt1s_r+c2(i)*kt2s_r+c3(i)*kt3s_r &
                 + c4(i)*kt4s_r+c5(i)*kt5s_r+rho(i)*rs1_r)*r(i)**2
          ! d^2 * T{d} in (D.47)
          jumpt(i) = -rho(i)*ts_r*r(i)*r(i)

          ! ============================================== !
          ! STORE KERNELS
          ! ============================================== !

          ! kap/mu/rho
          K_kap_s(i,s)  = ks_r
          K_mu_s(i,s)   = ms_r
          ! rho kernels with frequency extracted:
          K_Trho_s(i,s) = -ts_r
          K_Vrho_s(i,s) = rs2_r
          ! boundary kernels
          do ibnd=1,nbnd-1
            if (i.eq.bnd(ibnd)+1) then
              ! boundary kernels with frequency extracted:
              K_Td_s(ibnd,s) = -(jumpt(i)-jumpt(i-1))
              K_Vd_s(ibnd,s) = (jump(i)-jump(i-1))
            endif
          enddo
          if (i==NR) then
            ! boundary kernels with frequency extracted:
            K_Td_s(nbnd,s) = jumpt(i)
            K_Vd_s(nbnd,s) = jump(i)
          endif
        enddo
      enddo

      ! UNLIKE-TYPE COUPLING
    else
      do s=sbuild0,sbuild1
        ! check the selection rules
        if (l<abs(lp-s).or.l>(lp+s)) cycle
        if (mod(l+s+lp,2).eq.0) cycle
        b4 = b(l,s,lp,2,-1) ! B_{lsl'}^{2-}
        b5 = b(l,s,lp,1,-1) ! B_{lsl'}^{1-}
        b8 = b(l,lp,s,1,-1) ! B_{ll's}^{1-}
        b9 = b(lp,l,s,1,-1) ! B_{l'ls}^{1-}

        ! deal with integrals in potential term:
        do i=2,NR
          f  = (2.*u(i)-bl*v(i))/r(i)
          fp = (2.*up(i)-blp*vp(i))/r(i)
          ! imaginary part of D.51:
          gs1_i = 0.5*rho(i)*(u(i)*dwp(i)-du(i)*wp(i) + &
                              u(i)*wp(i)/r(i) - &
                              2.*f*wp(i))*b9/r(i) &
                - 0.5*rho(i)*(dw(i)*up(i)-w(i)*dup(i) + &
                              w(i)*up(i)/r(i) - &
                              2.*w(i)*fp)*b8/r(i)
          ! imaginary part of D.52:
          gs2_i = 0.5*rho(i)*u(i)*wp(i)*b9/r(i) &
                - 0.5*rho(i)*w(i)*up(i)*b8/r(i)
          ! integrands of D.55 /(r*r)
          k1_i(i) = (dble(s+1)*gs2_i-r(i)*gs1_i)*r(i)**(-s-2)
          k2_i(i) = (dble(s)*gs2_i+r(i)*gs1_i)*r(i)**(s-1)
        enddo
        k1_i(1) = k1_i(2)
        k2_i(1) = k2_i(2)
        ! calculate the substituted phi terms in D.55
        do i=1,NR
          call intgrl(i_1_i,r,i,NR,k1_i,s1,s2,s3)
          call intgrl(i_2_i,r,1,i,k2_i,s1,s2,s3)
          sub_i(i) = i_1_i*r(i)**s - i_2_i*r(i)**(-s-1)
          sub_i(i) = sub_i(i)*4.0/dble(2*s+1)
        enddo

        ! now compile the kernels
        do i=2,NR
          f  = (2.*u(i)-bl*v(i))/r(i)
          fp = (2.*up(i)-blp*vp(i))/r(i)
          x  = dv(i)-v(i)/r(i)+u(i)/r(i)
          xp = dvp(i)-vp(i)/r(i)+up(i)/r(i)
          z  = dw(i) - w(i)/r(i)
          zp = dwp(i) - wp(i)/r(i)
          ! imaginary part of D.49:
          ms_i  = -(x*zp-z*xp)*b5 &
                  -(v(i)*wp(i)-w(i)*vp(i))*b4/(r(i)*r(i))
          ! imaginary part of D.50:
          rs1_i = -((p(i)*wp(i)-w(i)*pp(i)) + &
                    0.5*g(i)*(u(i)*wp(i)-w(i)*up(i)))*b5/r(i)
          ! imaginary part of D.55:
          rs2_i = rs1_i + sub_i(i)
          ! imaginary part of D.46
          ts_i  = -(v(i)*wp(i)-w(i)*vp(i))*b5
          ! imaginary parts of boundary kernels
          ! (directly copied from HY's code with TI
          kt1s_i=du(i)*wp(i)*b8/r(i)-dup(i)*w(i)*b9/r(i)
          kt2s_i=0.5*(vp(i)*w(i)-wp(i)*v(i))*b4/(r(i)*r(i))
          kt4s_i=-f*wp(i)*b8/r(i)+fp*w(i)*b9/r(i)
          kt5s_i=(-xp*z+zp*x+dvp(i)*z-dv(i)*zp+dw(i)*xp-dwp(i)*x)*b5

          ! d^2*V{d} unrelated to physical dispersion
          jump(i) = -( c1(i)*kt1s_i+c2(i)*kt2s_i &
                     +c3(i)*kt4s_i+c5(i)*kt5s_i+ &
                     rho(i)*rs1_i)*r(i)**2
          ! r^2*T{d} in (D.47)
          jumpt(i) = -r(i)*r(i)*rho(i)*ts_i

          ! ============================================== !
          ! STORE KERNELS
          ! ============================================== !

          ! mu:
          K_mu_s(i,s)  = ii*ms_i
          ! rho kernels with frequency extracted:
          K_Trho_s(i,s) = -ii*ts_i
          K_Vrho_s(i,s) = ii*rs2_i

          ! boundary kernels
          do ibnd=1,nbnd-1
            if (i.eq.bnd(ibnd)+1) then
              ! boundary kernels with frequency extracted:
              K_Td_s(ibnd,s) = -ii*(jumpt(i)-jumpt(i-1))
              K_Vd_s(ibnd,s) = ii*(jump(i)-jump(i-1))
            endif
          enddo
          if (i.eq.NR) then
            ! boundary kernels with frequency extracted:
            K_Td_s(nbnd,s) = ii*jumpt(i)
            K_Vd_s(nbnd,s) = -ii*jump(i)
          endif
        enddo
      enddo
      endif

    ! ==================================================== !
    ! END CALCULATIONS AND DEALLOCATE
    ! ==================================================== !

    DEALLOCATE(sub_r,sub_i)
    DEALLOCATE(k1_r,k2_r,k1_i,k2_i,jump,jumpt,s1,s2,s3)
    NULLIFY(u, du, v, dv, w, dw, p, dp)
    NULLIFY(up,dup,vp,dvp,wp,dwp,pp,dpp)

    return

  END SUBROUTINE woodhouse_kernels

  ! ================================================================ !

  SUBROUTINE woodhouse_kernels_store(sbuild0,sbuild1,nmode,&
                                K_mu_wh,K_kp_wh,&
                                K_Tro_wh,K_Vro_wh,&
                                K_Td_wh,K_Vd_wh)


    ! routine will call subroutine woodhouse_kernels and store them
    ! for adjoint calculations.

    USE modellayer_module, ONLY:NR,NCMB,NMOHO,kdis
    USE prem_module
    USE mantle_model_module_DS
    USE futil_module, ONLY : intgrl

    IMPLICIT NONE
    include 'coupling.h'
    ! input/output:
    integer, intent(in) :: sbuild0,sbuild1,nmode
    complex*16, allocatable, dimension(:,:,:,:), intent(out) :: &
         K_mu_wh,K_kp_wh,&
         K_Tro_wh,K_Vro_wh
    complex*16, allocatable, dimension(:,:,:,:), intent(out) :: &
         K_Td_wh,K_Vd_wh
    ! complex*16, dimension(nmode,nmode,nr,0:sbuild1), intent(out) :: &
    !                                   K_mu_wh,K_kp_wh,&
    !                                   K_Tro_wh,K_Vro_wh
    ! complex*16, dimension(nmode,nmode,nbnd,0:sbuild1), intent(out) :: &
    !                                   K_Td_wh,K_Vd_wh
    ! internal:
    integer :: ieig,jeig
    complex*16, dimension(nr,0:sbuild1)   :: kmu_ij_s,&
                                             kkp_ij_s,&
                                             kTro_ij_s,&
                                             kVro_ij_s
    complex*16, dimension(nbnd,0:sbuild1) :: kTd_ij_s,&
                                             kVd_ij_s                                          

    ! initialize kernels:
    K_mu_wh(:,:,:,:)  = dcmplx(0.d0,0.d0)
    K_kp_wh(:,:,:,:)  = dcmplx(0.d0,0.d0)
    K_Tro_wh(:,:,:,:) = dcmplx(0.d0,0.d0)
    K_Vro_wh(:,:,:,:) = dcmplx(0.d0,0.d0)
    K_Td_wh(:,:,:,:)  = dcmplx(0.d0,0.d0)
    K_Vd_wh(:,:,:,:)  = dcmplx(0.d0,0.d0)

    ! begin build:
    do ieig=1,nmode
      do jeig=1,nmode
        write(*,*) '.... ... modes',ieig,'+',jeig,'of',nmode,'modes'
        call woodhouse_kernels(ieig,jeig,sbuild0,sbuild1,&
                               kmu_ij_s,kkp_ij_s,&
                               kTro_ij_s,kTd_ij_s,&
                               kVro_ij_s,kVd_ij_s)
        K_mu_wh(ieig,jeig,:,:)  = kmu_ij_s
        K_kp_wh(ieig,jeig,:,:)  = kkp_ij_s
        K_Tro_wh(ieig,jeig,:,:) = kTro_ij_s
        K_Vro_wh(ieig,jeig,:,:) = kVro_ij_s
        K_Td_wh(ieig,jeig,:,:)  = kTd_ij_s
        K_Vd_wh(ieig,jeig,:,:)  = kVd_ij_s
      enddo
    enddo

    return

  END SUBROUTINE woodhouse_kernels_store

  ! ================================================================ !
 
  SUBROUTINE rotation_integrals(ieig,jeig,int1,int2)

    ! routine returns the matrix elements for rotational coupling of two mutliplets
    ! result returned as a Coriolis and centrifugal contribuition
    ! note that as the matrix elements are non-zero only along the diagonal m = mp
    ! only these values are returned as a vector with dimensions (1:l+1)

    USE modellayer_module, ONLY:NR,NCMB,NMOHO
    USE prem_module
    USE futil_module, ONLY : intgrl,thrj  
    USE eigfun_module, ONLY : r=>eigfun_radius, &
                                 eigfun_u, eigfun_du, &
                                 eigfun_v, eigfun_dv, &
                                 eigfun_w, eigfun_dw, &
                                 eigfun_p, eigfun_dp, &
                                 eiginfo
  
    implicit none
    include "coupling.h"
  
    ! inputs
    integer, intent(in) :: ieig
    integer, intent(in) :: jeig

    ! outputs
    real*8, intent(out) :: int1,int2


    ! local variables
    integer :: np,typep,lp, n,type,l,m,mp,im1,im2
    integer :: smax,i,k,s,d,ds,ibnd,itmp
    real*8    :: omp,qp,gvp,avp,ahsp,ahtp
    real*8    :: om,q,gv,av,ahs,aht
    real*8    :: blp, bl,fac1,fac2
    real*8 :: fp,xp,zp,f,x,z,gs1,gs2
    real*8 :: b1,b2,b3,b4,b5,b6,b7,b8,b9
    complex*16 :: ii = cmplx(0.,1.)
  
    ! dimension(NR)
    real*8, POINTER :: up(:),dup(:),vp(:),dvp(:),wp(:),dwp(:),pp(:),dpp(:)
    real*8, POINTER :: u(:), du(:), v(:), dv(:), w(:), dw(:), p(:), dp(:)
    real*8, DIMENSION(:), ALLOCATABLE :: s1,s2,s3,intg1,intg2
  
    IF (.NOT.(ALLOCATED(r).AND.ALLOCATED(eigfun_u)&
         .AND.ALLOCATED(eiginfo))) STOP &
         'Variables in eigfun_module has to be defined previously'
  
    NULLIFY(u, du, v, dv, w, dw, p, dp)
    NULLIFY(up,dup,vp,dvp,wp,dwp,pp,dpp)
    
    ! define eigenfunctions for mode 1 (p)
    om=eiginfo(ieig)%w
    q=eiginfo(ieig)%q 
    n=eiginfo(ieig)%n
    l=eiginfo(ieig)%l
    type=eiginfo(ieig)%itype
    u=>eigfun_u(:,ieig)
    v=>eigfun_v(:,ieig)
    w=>eigfun_w(:,ieig)
    p=>eigfun_p(:,ieig)
    du=>eigfun_du(:,ieig)
    dv=>eigfun_dv(:,ieig)
    dw=>eigfun_dw(:,ieig)
    dp=>eigfun_dp(:,ieig)
    bl=l*(l+1.0)
    
    ! define eigenfunctions for mode 2
    omp=eiginfo(jeig)%w
    qp=eiginfo(jeig)%q 
    np=eiginfo(jeig)%n
    lp=eiginfo(jeig)%l
    typep=eiginfo(jeig)%itype
    up=>eigfun_u(:,jeig)
    vp=>eigfun_v(:,jeig)
    wp=>eigfun_w(:,jeig)
    pp=>eigfun_p(:,jeig)
    dup=>eigfun_du(:,jeig)
    dvp=>eigfun_dv(:,jeig)
    dwp=>eigfun_dw(:,jeig)
    dpp=>eigfun_dp(:,jeig)
    blp=lp*(lp+1.0)

    ! initialise values
    int1 = 0.0
    int2 = 0.0

    ! check the selection rules
    if(abs(lp-2) > l .or. l > lp+2) return
    
    IF (SIZE(r)/=NR .OR. SIZE(u)/=NR .OR. SIZE(w)/=NR .OR. &
         SIZE(up)/=NR.OR. SIZE(wp)/=NR) &
         STOP 'dimension of eigenfunctions has to be equal, in SUB &
         initialize_mantle STOP'
    ALLOCATE(s1(NR),s2(NR),s3(NR),intg1(NR),intg2(NR),STAT=itmp)
    IF (itmp/=0) STOP 'allocate error 0 in mantle.f, STOP '
    intg1 = 0.
    intg2 = 0.
           
    if(type.eq.typep) then ! like-type coupling

       ! set the Coriolis integrand
       intg1 = rho*(v*vp + u*vp + v*up + w*wp)
              
       ! set the centrifugal integrand
       b1=b(l,2,lp,0,1)
       b2=b(l,2,lp,2,1)
       b3=b(l,2,lp,1,1)
       b6=b(lp,l,2,1,1)
       b7=b(l,lp,2,1,1)
       do i=2,NR
          fp=(2.0*up(i)-blp*vp(i))/r(i)
          f=(2.0*u(i)-bl*v(i))/r(i)
          ! REAL part of V{phi} in (D.51) due to like-type coupling
          gs1=0.5*rho(i)*(u(i)*dvp(i)+u(i)*vp(i)/r(i) &
              -du(i)*vp(i)-2.0*f*vp(i))*b6/r(i) &
              +0.5*rho(i)*(up(i)*dv(i)+up(i)*v(i)/r(i) &
              -dup(i)*v(i)-2.0*fp*v(i))*b7/r(i) &
              +rho(i)*u(i)*up(i)*2*(2+1.0)*b1/(r(i)*r(i))
          ! REAL part of V{dphi/dr} in (D.52) due to like-type
          !  coupling
          gs2=0.5*rho(i)*u(i)*vp(i)*b6/r(i)+0.5*rho(i)*up(i)*v(i)*b7/r(i) &
              -rho(i)*(fp*u(i)+up(i)*f)*b1          
          ! store the integrand
          intg2(i) =  (1./3.)*r(i)**2*gs1 + (2./3.)*r(i)*gs2
       enddo
       ! perform the integrals
       call intgrl(int1,r,1,nr,intg1,s1,s2,s3)
       call intgrl(int2,r,1,nr,intg2,s1,s2,s3)              
      
    else ! unlike coupling
       ! set the Coriolis integrand
       intg1 =  0.5*(bl-blp-2.)*u*wp + 0.5*(bl-blp+2.)*w*up &
               +0.5*(bl+blp-2.)*(v*wp - w*vp)
       intg1 = rho*intg1
       ! set the centrifugal integrand
       b4=b(l,2,lp,2,-1)
       b5=b(l,2,lp,1,-1)
       b8=b(l,lp,2,1,-1)
       b9=b(lp,l,2,1,-1)
       do i=2,NR
          fp=(2.0*up(i)-blp*vp(i))/r(i)
          f=(2.0*u(i)-bl*v(i))/r(i)
          gs1=0.5*rho(i)*(u(i)*dwp(i)+u(i)*wp(i)/r(i) &
              -du(i)*wp(i)-2.0*f*wp(i))*b9/r(i) &
              -0.5*rho(i)*(up(i)*dw(i)+up(i)*w(i)/r(i) &
              -dup(i)*w(i)-2.0*fp*w(i))*b8/r(i)
          gs2=0.5*rho(i)*u(i)*wp(i)*b9/r(i)-0.5*rho(i)*up(i)*w(i)*b8/r(i)
          intg2(i) =  (1./3.)*r(i)**2*gs1 + (2./3.)*r(i)*gs2
       end do

       ! perform the integrals
       call intgrl(int1,r,1,nr,intg1,s1,s2,s3)
       call intgrl(int2,r,1,nr,intg2,s1,s2,s3)       
       
    end if    
    return

  END SUBROUTINE rotation_integrals

  ! ================================================================ !
  ! Changing adjoint_kernels routine such that they can be computed
  ! between starting and ending indeces. No integration is performed.
  ! For a routine that computes integration in the frequency domain,
  ! see adjoint_kernels_int.
  
  SUBROUTINE adjoint_kernels(nmode,ndim,istation,&
       ukf,uka,wf,dw,nf,&
       sbuild0,sbuild1,&
       swd0,swd1,nmwd,iwdhs,&
       K_mu_wh,K_kp_wh,&
       K_Tro_wh,K_Vro_wh,&
       K_Td_wh,K_Vd_wh,&
       kmu_out,kkp_out,&
       kro_out,kds_out,&
       ind_start, ind_end)

    ! Calculate adjoint kernels.
    USE futil_module, ONLY : thrj
    USE eigfun_module, ONLY : r=>eigfun_radius,eiginfo
    USE modellayer_module, ONLY : NR
    USE mantle_model_module_DS


    IMPLICIT NONE
    include 'coupling.h'
    ! input/output:
    integer, intent(in) :: nmode,ndim,istation,nf
    integer, intent(in) :: sbuild0,sbuild1
    integer, intent(in) :: swd0,swd1,nmwd
    integer, optional, intent(in) :: ind_start,ind_end
    integer, dimension(0:1,0:120,0:120) :: iwdhs
    real*8, intent(in) :: dw
    complex*16, dimension(nf), intent(in) :: wf
    complex*16, dimension(nf,ndim), intent(in) :: ukf,uka
    complex*16, dimension(nmwd,nmwd,nr,swd0:swd1), intent(in) :: &
                                                  K_mu_wh,K_kp_wh,&
                                                  K_Tro_wh,K_Vro_wh
    complex*16, dimension(nmwd,nmwd,nbnd,swd0:swd1), intent(in) :: &
                                                  K_Td_wh,K_Vd_wh
    complex*16, dimension(nr,0:sbuild1,-sbuild1:sbuild1,nf), intent(out) :: &
                                                 kmu_out,kkp_out,&
                                                 kro_out
    complex*16, dimension(nbnd,0:sbuild1,-sbuild1:sbuild1,nf), intent(out) :: &
                                                 kds_out
    ! internal:
    character(len=1) :: sym,symp
    integer :: ieig,jeig,idxi,idxj
    integer :: ivec1,ivec2,ivec
    integer :: jvec1,jvec2,jvec
    integer :: l,lp,m,mp,s,t,ist,i1,i2
    integer :: id,i,ibnd,type,typep,n,np,ifreq
    integer :: start_loop, end_loop
    complex*16, dimension(NR) :: kmu,kkp,kro,krot,krov
    complex*16, dimension(nbnd) :: kd,kdt,kdv
    real, parameter :: eps = 1.e-7
    real :: nus,nul,nulp,nult,rho1,rho2,fac,fac1,fac2
    complex :: jump

    complex :: wki,wkj !!bm!!

    ! Dealing with optional arguments to get the frequency bounds of computation
    start_loop = 1
    end_loop = nf
    if (present(ind_start)) start_loop = ind_start
    if (present(ind_end)) end_loop = ind_end

    write(*,*)"Computing kernels between indeces", start_loop, end_loop
    ! initialize the kernels
    kmu_out(:,:,:,:) = dcmplx(0.d0,0.d0)
    kkp_out(:,:,:,:) = dcmplx(0.d0,0.d0)
    kro_out(:,:,:,:) = dcmplx(0.d0,0.d0)
    kds_out(:,:,:,:) = dcmplx(0.d0,0.d0)

    ivec2 = 0
    do ieig = 1,nmode
      ! degree of multiplet
       l = eiginfo(ieig)%l
       wki = eiginfo(ieig)%w !!bm!!
       !nul = sqrt((2.*l+1.))/4.*pi !!! BM -- is this right ( 4*pi -> (4*pi) )
      nul = sqrt((2.*l+1.)/(4.*pi)) !!! BM -- is this right ( 4*pi -> (4*pi) ) !!! HL UPDATE!!!
      ivec1 = ivec2 + 1
      ivec2 = ivec1 + 2*l
      type = eiginfo(ieig)%itype
      n    = eiginfo(ieig)%n
      if (type == 0) then
        sym = 'T'
      else
        sym = 'S'
      endif

      ! loop over the forward multiplets
      jvec2 = 0
      do jeig = 1,nmode
        ! degree of multiplet
         lp = eiginfo(jeig)%l
         wkj = eiginfo(jeig)%w !!bm!!
        nulp = sqrt((2*lp+1.)/(4.*pi))
        jvec1 = jvec2 + 1
        jvec2 = jvec1 + 2*lp
        typep = eiginfo(jeig)%itype
        np    = eiginfo(jeig)%n
        if(typep == 0) then
          symp = 'T'
        else
          symp = 'S'
        end if
        ! print *, '(',n,sym,l,')  x  (',np,symp,lp,')'

        ! loop over structural degrees
        do s = sbuild0,sbuild1 
           nus = sqrt((2*s+1)/(4.*pi))
           fac = 4.*pi*nul*nus*nulp
           idxi    = iwdhs(type,n,l)
           idxj    = iwdhs(typep,np,lp)
           kmu(:)  = K_mu_wh(idxi,idxj,:,s)
           kkp(:)  = K_kp_wh(idxi,idxj,:,s)
           krot(:) = K_Tro_wh(idxi,idxj,:,s)
           krov(:) = K_Vro_wh(idxi,idxj,:,s)
           kdt(:)  = K_Td_wh(idxi,idxj,:,s)
           kdv(:)  = K_Vd_wh(idxi,idxj,:,s)
           ! begin loop over structural orders
           do t = -s,s
              ! begin loop over singlets for first multiplet
              do m = -l,l
                 ivec = ivec1+l+m
                 fac1 = (-1)**m*fac   
                 ! begin loop over singlets for the second multiplet
                 do mp = -lp,lp
                    jvec = jvec1+lp+mp
                    ! get the geometric terms
                    fac2 = fac1*thrj(l,s,lp,-m,t,mp)
                    ! add in contributions to the kernels and begin
                    ! frequency integral
                    !!bm!! BM23: do not integrate
                    ! do ifreq=1,nf
                    ! Adding the possibility to choose bounds of computation
                    do ifreq=start_loop,end_loop
                       ! HL - krot already has minus in it.
                       ! write(*,*) 'adjoint f',wf(ifreq)
                       kro(:) = wf(ifreq)*wf(ifreq)*krot(:) + krov(:)
                       kd(:)  = wf(ifreq)*wf(ifreq)*kdt(:) + kdv(:)
                       kmu_out(:,s,t,ifreq) = kmu_out(:,s,t,ifreq) + &
                            fac2*conjg(kmu)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec)) 
                       kkp_out(:,s,t,ifreq) = kkp_out(:,s,t,ifreq) + &
                            fac2*conjg(kkp)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))
                       kro_out(:,s,t,ifreq) = kro_out(:,s,t,ifreq) + &
                            fac2*conjg(kro)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))
                       kds_out(:,s,t,ifreq) = kds_out(:,s,t,ifreq) + &
                            fac2*conjg(kd)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))
                    enddo
                    ! end integration over frequency
                 enddo
                 ! end loop over singlets for the second multiplet
              enddo
              ! end loop over singlets for first multiplet
           enddo
           ! end loop over structural orders
        enddo
        ! end loop over structural degrees
     enddo
     ! end loop over forward multiplets
  enddo
  ! end loop over adjoint multiplets

return

  END SUBROUTINE adjoint_kernels

  ! ================================================================ !
  
!   SUBROUTINE adjoint_kernels(nmode,ndim,istation,&
!                              ukf,uka,wf,dw,nf,&
!                              sbuild0,sbuild1,&
!                              swd0,swd1,nmwd,iwdhs,&
!                              K_mu_wh,K_kp_wh,&
!                              K_Tro_wh,K_Vro_wh,&
!                              K_Td_wh,K_Vd_wh,&
!                              kmu_out,kkp_out,&
!                              kro_out,kds_out)

!     ! Calculate adjoint kernels.
!     USE futil_module, ONLY : thrj
!     USE eigfun_module, ONLY : r=>eigfun_radius,eiginfo
!     USE modellayer_module, ONLY : NR
!     USE mantle_model_module_DS


!     IMPLICIT NONE
!     include 'coupling.h'
!     ! input/output:
!     integer, intent(in) :: nmode,ndim,istation,nf
!     integer, intent(in) :: sbuild0,sbuild1
!     integer, intent(in) :: swd0,swd1,nmwd
!     integer, dimension(0:1,0:120,0:120) :: iwdhs
!     real*8, intent(in) :: dw
!     complex*16, dimension(nf), intent(in) :: wf
!     complex*16, dimension(nf,ndim), intent(in) :: ukf,uka
!     complex*16, dimension(nmwd,nmwd,nr,swd0:swd1), intent(in) :: &
!                                                   K_mu_wh,K_kp_wh,&
!                                                   K_Tro_wh,K_Vro_wh
!     complex*16, dimension(nmwd,nmwd,nbnd,swd0:swd1), intent(in) :: &
!                                                   K_Td_wh,K_Vd_wh
!     complex*16, dimension(nr,0:sbuild1,-sbuild1:sbuild1,nf), intent(out) :: &
!                                                  kmu_out,kkp_out,&
!                                                  kro_out
!     complex*16, dimension(nbnd,0:sbuild1,-sbuild1:sbuild1,nf), intent(out) :: &
!                                                  kds_out
!     ! internal:
!     character(len=1) :: sym,symp
!     integer :: ieig,jeig,idxi,idxj
!     integer :: ivec1,ivec2,ivec
!     integer :: jvec1,jvec2,jvec
!     integer :: l,lp,m,mp,s,t,ist,i1,i2
!     integer :: id,i,ibnd,type,typep,n,np,ifreq
!     complex*16, dimension(NR) :: kmu,kkp,kro,krot,krov
!     complex*16, dimension(nbnd) :: kd,kdt,kdv
!     real, parameter :: eps = 1.e-7
!     real :: nus,nul,nulp,nult,rho1,rho2,fac,fac1,fac2
!     complex :: jump

!     complex :: wki,wkj !!bm!!

!     ! initialize the kernels
!     kmu_out(:,:,:,:) = dcmplx(0.d0,0.d0)
!     kkp_out(:,:,:,:) = dcmplx(0.d0,0.d0)
!     kro_out(:,:,:,:) = dcmplx(0.d0,0.d0)
!     kds_out(:,:,:,:) = dcmplx(0.d0,0.d0)

!     ivec2 = 0
!     do ieig = 1,nmode
!       ! degree of multiplet
!        l = eiginfo(ieig)%l
!        wki = eiginfo(ieig)%w !!bm!!
!        !nul = sqrt((2.*l+1.))/4.*pi !!! BM -- is this right ( 4*pi -> (4*pi) )
!       nul = sqrt((2.*l+1.)/(4.*pi)) !!! BM -- is this right ( 4*pi -> (4*pi) ) !!! HL UPDATE!!!
!       ivec1 = ivec2 + 1
!       ivec2 = ivec1 + 2*l
!       type = eiginfo(ieig)%itype
!       n    = eiginfo(ieig)%n
!       if (type == 0) then
!         sym = 'T'
!       else
!         sym = 'S'
!       endif

!       ! loop over the forward multiplets
!       jvec2 = 0
!       do jeig = 1,nmode
!         ! degree of multiplet
!          lp = eiginfo(jeig)%l
!          wkj = eiginfo(jeig)%w !!bm!!
!         nulp = sqrt((2*lp+1.)/(4.*pi))
!         jvec1 = jvec2 + 1
!         jvec2 = jvec1 + 2*lp
!         typep = eiginfo(jeig)%itype
!         np    = eiginfo(jeig)%n
!         if(typep == 0) then
!           symp = 'T'
!         else
!           symp = 'S'
!         end if
!         ! print *, '(',n,sym,l,')  x  (',np,symp,lp,')'

!         ! loop over structural degrees
!         do s = sbuild0,sbuild1 
!           nus = sqrt((2*s+1)/(4.*pi))
!           fac = 4.*pi*nul*nus*nulp
!           idxi    = iwdhs(type,n,l)
!           idxj    = iwdhs(typep,np,lp)
!           kmu(:)  = K_mu_wh(idxi,idxj,:,s)
!           kkp(:)  = K_kp_wh(idxi,idxj,:,s)
!           krot(:) = K_Tro_wh(idxi,idxj,:,s)
!           krov(:) = K_Vro_wh(idxi,idxj,:,s)
!           kdt(:)  = K_Td_wh(idxi,idxj,:,s)
!           kdv(:)  = K_Vd_wh(idxi,idxj,:,s)
!           ! begin loop over structural orders
!           do t = -s,s
!             ! begin loop over singlets for first multiplet
!             do m = -l,l
!               ivec = ivec1+l+m
!               fac1 = (-1)**m*fac   
!               ! begin loop over singlets for the second multiplet
!               do mp = -lp,lp
!                 jvec = jvec1+lp+mp
!                 ! get the geometric terms
!                 fac2 = fac1*thrj(l,s,lp,-m,t,mp)
!                 ! add in contributions to the kernels and begin
!                 ! frequency integral
!                 !!bm!! BM23: do not integrate
!                 do ifreq=1,nf
!                    ! HL - krot already has minus in it.
! !                   write(*,*) 'adjoint f',wf(ifreq)
!                   kro(:) = wf(ifreq)*wf(ifreq)*krot(:) + krov(:)
!                   kd(:)  = wf(ifreq)*wf(ifreq)*kdt(:) + kdv(:)
!                   kmu_out(:,s,t,ifreq) = kmu_out(:,s,t,ifreq) + &
!                     fac2*conjg(kmu)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec)) 
!                   kkp_out(:,s,t,ifreq) = kkp_out(:,s,t,ifreq) + &
!                     fac2*conjg(kkp)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))
!                   kro_out(:,s,t,ifreq) = kro_out(:,s,t,ifreq) + &
!                     fac2*conjg(kro)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))
!                   kds_out(:,s,t,ifreq) = kds_out(:,s,t,ifreq) + &
!                        fac2*conjg(kd)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))
!                 enddo
!                 ! end integration over frequency
!               enddo
!               ! end loop over singlets for the second multiplet
!             enddo
!             ! end loop over singlets for first multiplet
!           enddo
!           ! end loop over structural orders
!         enddo
!         ! end loop over structural degrees
!       enddo
!       ! end loop over forward multiplets
!     enddo
!     ! end loop over adjoint multiplets

!     return

!   END SUBROUTINE adjoint_kernels

! ================================================================ !

  SUBROUTINE adjoint_kernels_int(nmode,ndim,istation,&
       ukf,uka,wf,dw,nf,&
       sbuild0,sbuild1,&
       swd0,swd1,nmwd,iwdhs,&
       K_mu_wh,K_kp_wh,&
       K_Tro_wh,K_Vro_wh,&
       K_Td_wh,K_Vd_wh,&
       kmu_out,kkp_out,&
       kro_out,kds_out,&
       ind_start,ind_end)

    ! <SA> Changing the list of arguments to get the start and end of
    ! frequency integration. Useful when one wants to isolate a mode
    ! specifically. Indeces, to be computed and found prior the call.
    ! </SA>
    
    ! Calculate adjoint kernels.
    USE futil_module, ONLY : thrj
    USE eigfun_module, ONLY : r=>eigfun_radius,eiginfo
    USE modellayer_module, ONLY : NR
    USE mantle_model_module_DS


    IMPLICIT NONE
    include 'coupling.h'
    ! input/output:
    integer, intent(in) :: nmode,ndim,istation,nf
    integer, intent(in) :: sbuild0,sbuild1
    integer, intent(in) :: swd0,swd1,nmwd
    integer, optional, intent(in) :: ind_start,ind_end
    integer, dimension(0:1,0:120,0:120) :: iwdhs
    real*8, intent(in) :: dw
    complex*16, dimension(nf), intent(in) :: wf
    complex*16, dimension(nf,ndim), intent(in) :: ukf,uka
    complex*16, dimension(nmwd,nmwd,nr,swd0:swd1), intent(in) :: &
                                                  K_mu_wh,K_kp_wh,&
                                                  K_Tro_wh,K_Vro_wh
    complex*16, dimension(nmwd,nmwd,nbnd,swd0:swd1), intent(in) :: &
                                                  K_Td_wh,K_Vd_wh
    complex*16, dimension(nr,0:sbuild1,-sbuild1:sbuild1), intent(out) :: &
                                                 kmu_out,kkp_out,&
                                                 kro_out
    complex*16, dimension(nbnd,0:sbuild1,-sbuild1:sbuild1), intent(out) :: &
                                                 kds_out
    ! internal:
    character(len=1) :: sym,symp
    integer :: ieig,jeig,idxi,idxj
    integer :: ivec1,ivec2,ivec
    integer :: jvec1,jvec2,jvec
    integer :: l,lp,m,mp,s,t,ist,i1,i2
    integer :: id,i,ibnd,type,typep,n,np,ifreq
    integer :: start_loop, end_loop
    complex*16, dimension(NR) :: kmu,kkp,kro,krot,krov
    complex*16, dimension(nbnd) :: kd,kdt,kdv
    real, parameter :: eps = 1.e-7
    real :: nus,nul,nulp,nult,rho1,rho2,fac,fac1,fac2
    complex :: jump

    complex :: wki,wkj !!bm!!

    ! Dealing with optional arguments to get the frequency bounds of computation
    start_loop = 1
    end_loop = nf
    if (present(ind_start)) start_loop = ind_start
    if (present(ind_end)) end_loop = ind_end

    ! initialize the kernels
    kmu_out(:,:,:) = dcmplx(0.d0,0.d0)
    kkp_out(:,:,:) = dcmplx(0.d0,0.d0)
    kro_out(:,:,:) = dcmplx(0.d0,0.d0)
    kds_out(:,:,:) = dcmplx(0.d0,0.d0)

    ivec2 = 0
    do ieig = 1,nmode
      ! degree of multiplet
       l = eiginfo(ieig)%l
       wki = eiginfo(ieig)%w !!bm!!
       !nul = sqrt((2.*l+1.))/4.*pi !!! BM -- is this right ( 4*pi -> (4*pi) )
      nul = sqrt((2.*l+1.)/(4.*pi)) !!! BM -- is this right ( 4*pi -> (4*pi) ) !!! HL UPDATE!!!
      ivec1 = ivec2 + 1
      ivec2 = ivec1 + 2*l
      type = eiginfo(ieig)%itype
      n    = eiginfo(ieig)%n
      if (type == 0) then
        sym = 'T'
      else
        sym = 'S'
      endif

      ! loop over the forward multiplets
      jvec2 = 0
      do jeig = 1,nmode
        ! degree of multiplet
         lp = eiginfo(jeig)%l
         wkj = eiginfo(jeig)%w !!bm!!
        nulp = sqrt((2*lp+1.)/(4.*pi))
        jvec1 = jvec2 + 1
        jvec2 = jvec1 + 2*lp
        typep = eiginfo(jeig)%itype
        np    = eiginfo(jeig)%n
        if(typep == 0) then
          symp = 'T'
        else
          symp = 'S'
        end if
        print *, '(',n,sym,l,')  x  (',np,symp,lp,')'

        ! loop over structural degrees
        do s = sbuild0,sbuild1 
          nus = sqrt((2*s+1)/(4.*pi))
          fac = 4.*pi*nul*nus*nulp
          idxi    = iwdhs(type,n,l)
          idxj    = iwdhs(typep,np,lp)
          kmu(:)  = K_mu_wh(idxi,idxj,:,s)
          kkp(:)  = K_kp_wh(idxi,idxj,:,s)
          krot(:) = K_Tro_wh(idxi,idxj,:,s)
          krov(:) = K_Vro_wh(idxi,idxj,:,s)
          kdt(:)  = K_Td_wh(idxi,idxj,:,s)
          kdv(:)  = K_Vd_wh(idxi,idxj,:,s)
          ! begin loop over structural orders
          do t = -s,s
             ! print*, "Structural degree ", t
            ! begin loop over singlets for first multiplet
            do m = -l,l
              ivec = ivec1+l+m
              fac1 = (-1)**m*fac   
              ! begin loop over singlets for the second multiplet
              do mp = -lp,lp
                jvec = jvec1+lp+mp
                ! get the geometric terms
                fac2 = fac1*thrj(l,s,lp,-m,t,mp)
                ! add in contributions to the kernels and begin
                ! frequency integral
                
                ! integrate
                do ifreq=start_loop,end_loop
                  kro(:) = wf(ifreq)*wf(ifreq)*krot(:) + krov(:)
                  kd(:)  = wf(ifreq)*wf(ifreq)*kdt(:) + kdv(:)
                  kmu_out(:,s,t) = kmu_out(:,s,t) + &
                    fac2*conjg(kmu)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))*dw
                  kkp_out(:,s,t) = kkp_out(:,s,t) + &
                    fac2*conjg(kkp)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))*dw
                  kro_out(:,s,t) = kro_out(:,s,t) + &
                    fac2*conjg(kro)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))*dw
                  kds_out(:,s,t) = kds_out(:,s,t) + &
                       fac2*conjg(kd)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))*dw
                enddo
                ! end integration over frequency
              enddo
              ! end loop over singlets for the second multiplet
            enddo
            ! end loop over singlets for first multiplet
          enddo
          ! end loop over structural orders
        enddo
        ! end loop over structural degrees
      enddo
      ! end loop over forward multiplets
    enddo
    ! end loop over adjoint multiplets

    return

  END SUBROUTINE adjoint_kernels_int


  ! ================================================================ !
    SUBROUTINE adjoint_kernels_test(iftest,nftest,&
                             nmode,ndim,istation,&
                             ukf,uka,wf,dw,nf,&
                             sbuild0,sbuild1,&
                             swd0,swd1,nmwd,iwdhs,&
                             K_mu_wh,K_kp_wh,&
                             K_Tro_wh,K_Vro_wh,&
                             K_Td_wh,K_Vd_wh,&
                             kmu_out,kkp_out,&
                             kro_out,kds_out)

    ! Calculate adjoint kernels.
    USE futil_module, ONLY : thrj
    USE eigfun_module, ONLY : r=>eigfun_radius,eiginfo
    USE modellayer_module, ONLY : NR
    USE mantle_model_module_DS


    IMPLICIT NONE
    include 'coupling.h'
    ! input/output:
    integer, intent(in) :: nftest
    integer, dimension(nftest), intent(in) :: iftest
    integer, intent(in) :: nmode,ndim,istation,nf
    integer, intent(in) :: sbuild0,sbuild1
    integer, intent(in) :: swd0,swd1,nmwd
    integer, dimension(0:1,0:120,0:120) :: iwdhs
    real*8, intent(in) :: dw
    complex*16, dimension(nf), intent(in) :: wf
    complex*16, dimension(nf,ndim), intent(in) :: ukf,uka
    complex*16, dimension(nmwd,nmwd,nr,swd0:swd1), intent(in) :: &
                                                  K_mu_wh,K_kp_wh,&
                                                  K_Tro_wh,K_Vro_wh
    complex*16, dimension(nmwd,nmwd,nbnd,swd0:swd1), intent(in) :: &
                                                  K_Td_wh,K_Vd_wh
    complex*16, dimension(nr,0:sbuild1,-sbuild1:sbuild1,nftest), intent(out) :: &
                                                 kmu_out,kkp_out,&
                                                 kro_out
    complex*16, dimension(nbnd,0:sbuild1,-sbuild1:sbuild1,nftest), intent(out) :: &
                                                 kds_out
    ! internal:
    character(len=1) :: sym,symp
    integer :: ieig,jeig,idxi,idxj,itest
    integer :: ivec1,ivec2,ivec
    integer :: jvec1,jvec2,jvec
    integer :: l,lp,m,mp,s,t,ist,i1,i2
    integer :: id,i,ibnd,type,typep,n,np,ifreq
    complex*16, dimension(NR) :: kmu,kkp,kro,krot,krov
    complex*16, dimension(nbnd) :: kd,kdt,kdv
    real, parameter :: eps = 1.e-7
    real :: nus,nul,nulp,nult,rho1,rho2,fac,fac1,fac2
    complex :: jump

    complex :: wki,wkj !!bm!!

    ! initialize the kernels
    kmu_out(:,:,:,:) = dcmplx(0.d0,0.d0)
    kkp_out(:,:,:,:) = dcmplx(0.d0,0.d0)
    kro_out(:,:,:,:) = dcmplx(0.d0,0.d0)
    kds_out(:,:,:,:) = dcmplx(0.d0,0.d0)

    ivec2 = 0
    do ieig = 1,nmode
      ! degree of multiplet
       l = eiginfo(ieig)%l
       wki = eiginfo(ieig)%w !!bm!!
       !nul = sqrt((2.*l+1.))/4.*pi !!! BM -- is this right ( 4*pi -> (4*pi) )
      nul = sqrt((2.*l+1.)/(4.*pi)) !!! BM -- is this right ( 4*pi -> (4*pi) ) !!! HL UPDATE!!!
      ivec1 = ivec2 + 1
      ivec2 = ivec1 + 2*l
      type = eiginfo(ieig)%itype
      n    = eiginfo(ieig)%n
      if (type == 0) then
        sym = 'T'
      else
        sym = 'S'
      endif

      ! loop over the forward multiplets
      jvec2 = 0
      do jeig = 1,nmode
        ! degree of multiplet
         lp = eiginfo(jeig)%l
         wkj = eiginfo(jeig)%w !!bm!!
        nulp = sqrt((2*lp+1.)/(4.*pi))
        jvec1 = jvec2 + 1
        jvec2 = jvec1 + 2*lp
        typep = eiginfo(jeig)%itype
        np    = eiginfo(jeig)%n
        if(typep == 0) then
          symp = 'T'
        else
          symp = 'S'
        end if
        print *, '(',n,sym,l,')  x  (',np,symp,lp,')'

        ! loop over structural degrees
        do s = sbuild0,sbuild1 
          nus = sqrt((2*s+1)/(4.*pi))
          fac = 4.*pi*nul*nus*nulp
          idxi    = iwdhs(type,n,l)
          idxj    = iwdhs(typep,np,lp)
          kmu(:)  = K_mu_wh(idxi,idxj,:,s)
          kkp(:)  = K_kp_wh(idxi,idxj,:,s)
          krot(:) = K_Tro_wh(idxi,idxj,:,s)
          krov(:) = K_Vro_wh(idxi,idxj,:,s)
          kdt(:)  = K_Td_wh(idxi,idxj,:,s)
          kdv(:)  = K_Vd_wh(idxi,idxj,:,s)
          ! begin loop over structural orders
          do t = -s,s
            ! begin loop over singlets for first multiplet
            do m = -l,l
              ivec = ivec1+l+m
              fac1 = (-1)**m*fac   
              ! begin loop over singlets for the second multiplet
              do mp = -lp,lp
                jvec = jvec1+lp+mp
                ! get the geometric terms
                fac2 = fac1*thrj(l,s,lp,-m,t,mp)
                ! add in contributions to the kernels and begin
                ! frequency integral
                !!bm!! BM23: do not integrate
                do ifreq=1,nftest

                   itest = iftest(ifreq)
                   ! HL - krot already has minus in it.
!                   write(*,*) 'adjoint f',wf(ifreq)
                   kro(:) = wf(itest)*wf(itest)*krot(:) + krov(:)
                   kd(:)  = wf(itest)*wf(itest)*kdt(:) + kdv(:)
                   ! kmu_out(:,s,t,ifreq) = kmu_out(:,s,t,ifreq) + &
                   !      fac2*conjg(kmu)*uka(itest,ivec)*conjg(ukf(itest,jvec)) 
                   ! kkp_out(:,s,t,ifreq) = kkp_out(:,s,t,ifreq) + &
                   !      fac2*conjg(kkp)*uka(itest,ivec)*conjg(ukf(itest,jvec))
                   ! kro_out(:,s,t,ifreq) = kro_out(:,s,t,ifreq) + &
                   !      fac2*conjg(kro)*uka(itest,ivec)*conjg(ukf(itest,jvec))
                   ! kds_out(:,s,t,ifreq) = kds_out(:,s,t,ifreq) + &
                   !      fac2*conjg(kd)*uka(itest,ivec)*conjg(ukf(itest,jvec))
                   kmu_out(:,s,t,ifreq) = kmu_out(:,s,t,ifreq) + &
                        fac2*conjg(kmu)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec)) 
                   kkp_out(:,s,t,ifreq) = kkp_out(:,s,t,ifreq) + &
                        fac2*conjg(kkp)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))
                   kro_out(:,s,t,ifreq) = kro_out(:,s,t,ifreq) + &
                        fac2*conjg(kro)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))
                   kds_out(:,s,t,ifreq) = kds_out(:,s,t,ifreq) + &
                        fac2*conjg(kd)*uka(ifreq,ivec)*conjg(ukf(ifreq,jvec))
                enddo
                ! end integration over frequency
              enddo
              ! end loop over singlets for the second multiplet
            enddo
            ! end loop over singlets for first multiplet
          enddo
          ! end loop over structural orders
        enddo
        ! end loop over structural degrees
      enddo
      ! end loop over forward multiplets
    enddo
    ! end loop over adjoint multiplets

    return

  END SUBROUTINE adjoint_kernels_test

! ================================================================ !

  FUNCTION b(lp,s,l,n,sign)

    !   Calculating B(l',s,l,N,sign) defined in (D.42) 
    !   For like-type coupling (i.e. toroidal-toridal,
    !      spheriodal-spheroidal) => B(l,s,l',1 or 2, +) are used
    !   For cross-type coupling (i.e. toridal-spheroidal),
    !      => B(l,s,l',1 or 2, -) are used

    USE futil_module, ONLY : thrj
    IMPLICIT NONE
    REAL*8 :: b
    INTEGER, INTENT(IN) :: lp,s,l,n,sign

    INTEGER :: i,fac1,fac2
    !      REAL :: thrj
    !
    IF ((sign/= -1) .AND. (sign/=1)) THEN
       STOP 'sign in funcion b is out of range, STOP'
    END IF

    fac1=1
    i=lp-n+1
    do while(i.le.(lp+n))
       fac1=fac1*i
       i=i+1
    enddo
    fac2=1
    i=l-n+1
    do while(i.le.(l+n))
       fac2=fac2*i
       i=i+1
    enddo

    b=0.5*(1.0+sign*((-1.0)**(lp+s+l)))*sqrt(float(fac1)*float(fac2)) &
       *((-1.0)**n)*thrj(lp,s,l,-n,0,n)

    return

  END FUNCTION b

! ================================================================ !

  FUNCTION delta(i,j)

    integer, intent(in) :: i,j
    real*8 :: delta
    delta = 0.
    if(i == j) delta = 1.
    return

  end FUNCTION delta

! ================================================================ !

  FUNCTION slm(l,m)
    !   S_lm (D.69)
    !
    IMPLICIT NONE
    REAL*8 :: slm
    INTEGER, INTENT(IN) :: l,m
    !
    slm=sqrt(float((l+m)*(l-m))/((2.0*l+1.0)*(2.0*l-1.0)))
    return
  END FUNCTION slm

! ================================================================ !

END MODULE splitting_seismo_module

