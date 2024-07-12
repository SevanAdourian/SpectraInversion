MODULE anelastic_module

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                                          !!
  !! HL 03/2018                                               !!
  !!                                                          !!
  !! This small subroutine treats the frequency dependence    !!
  !! of physical dispersion and attenuation.                  !!
  !! As we solve the system using the Direct Solution method, !!         
  !! we are solving the system at a fixed frequency omega     !!
  !! (or 'wf').                                               !!
  !!                                                          !!
  !! Some notes to clarify details:                           !!
  !!                                                          !!
  !! (1) The reference model is read in from the output in    !!
  !!     the MINEOS format.                                   !!
  !!     Thus elastic moduli imported are evaluated at wref   !!
  !!     (the reference frequency of the reference model).    !!
  !!     This includes Qkappa and Qmu also.                   !!
  !!                                                          !!
  !! (2) The degenerate eigenfrequencies (om(k) and q(k))     !!
  !!     include 1D physical dispersion.  This is calculated  !!
  !!     within MINEOS and so MINEOS's subroutine must be     !!
  !!     consistent with the frequency dependence adopted     !!
  !!     here.                                                !!
  !!                                                          !!
  !! (3) 3D dispersive effects are taken care of by adjusting !!
  !!     the elastic parameters when forming the kernels      !!
  !!     in rot_ell_module_FC and mantle_module_FC            !!
  !!                                                          !!
  !! (4) 3D attenuation effects are taken care of by forming  !!
  !!     the kernel in mantle_module_FC.  The 3D Q model is   !!
  !!     of the same format as the 3D vs/vb/rho model and     !!
  !!     imported in a similar way                            !!
  !!                                                          !!
  !! (5) The system to invert is thus:                        !!
  !!     S = -wf*wf * T + wf * W + V(wf) + i*A                !!
  !!                                                          !!
  !!                                                          !!
  !! note that 1D elastic basis model can be transversely     !!
  !! isotropic but the anelasticity is just isotropic         !!
  !! kappa = (1/9)*(C+4A-4N+4F)     [DT 8.191]                !!
  !! mu = (1/15)*(C+A+6L+5N-2F)     [DT 8.192]                !!
  !! but anelasticity is applied to isotopic:                 !!
  !! C=A=kappa + 4/3*mu                                       !!
  !! L=N=mu                                                   !!
  !! F=kappa-2/3*mu                                           !!
  !!                                                          !!
  !!                                                          !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  include "coupling.h"
  PRIVATE
  PUBLIC :: freq_dispersion, &
            freq_attenuation

CONTAINS

  SUBROUTINE freq_dispersion(omega,fct)

    ! Use this subroutine to calculate dispersive 
    ! correction.  It will multiply the (dimensionless)
    ! elastic moduli to map it from wref to omega 
    ! (or 'wf')
    !
    ! wref must be provided in coupling.h and match 
    ! the reference frequency in MINEOS
    !
    ! The frequency dependent form must match that 
    ! provided in MINEOS (see the preamble in 
    ! mineos/minos_bran_HY_HL.f)
    !
    ! e.g., mu(w) = mu(wref) + (mu(wref)/Q(wref))*fct
    ! the second term fct is what is returned here

    IMPLICIT NONE
    include "coupling.h"
    
    REAL, INTENT(IN) :: omega
    REAL :: wm,alpha
    REAL, INTENT(OUT) :: fct

    real :: om
    
    wm=1.9415e-3
    alpha=0.3

    ! DA modification to allow for negative frequencies
    ! check inequalities!
    om = abs(omega)*f_norm
    
    if (om.ge.wm) then
       fct = 2.*(log10(om/wref))/PI
    else
       fct = 2.*(log10(wm/wref)) + ((1.-(wm/om))**alpha)/alpha
    endif
    
    return
    
  END SUBROUTINE freq_dispersion

  SUBROUTINE freq_attenuation(omega,fct)

    ! Use this subroutine to calculate the imaginary
    ! part of the modulus due to attenuation. It will
    ! be scaled from the Q (either 1D or 3D, evaluated
    ! at wref) to Q at omega (or "wf").
    ! 
    ! For 1D Q this will already be accounted for by
    ! MINEOS in the imaginary part of the degenerate 
    ! eigenvalues.  For the 3D part it will form the
    ! kernels for Qinvk and Qinvm.  Make sure the
    ! two mappings are consistent.
    !
    ! The reference frequency must be provided in 
    ! coupling.h and called wref
    !
    ! e.g., mu(w) = Re[mu(w)] + i* fct*mu(wref)/Q(wref)
    ! the second term (fct) is what is returned here

    IMPLICIT NONE
    include "coupling.h"

    REAL, INTENT(IN) :: omega
    REAL :: wm,alpha
    REAL, INTENT(OUT) :: fct
    
    real :: om
    
   
    wm=1.9415e-3
    alpha=0.3

    ! DA modification to allow for negative frequencies
    om = abs(omega)*f_norm
        
    if (om.ge.wm) then
       fct = 1.
    else
       fct = (wm/om)**alpha
    endif

    if(omega < 0.) fct = -fct
    
    return
    
  END SUBROUTINE freq_attenuation
  
END MODULE anelastic_module
