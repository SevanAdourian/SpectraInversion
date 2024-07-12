MODULE nrtype

  ! ================================================================ !
  ! Predeclare and define values                                     !
  ! ================================================================ !

  ! symbolic names for kind types of 4-, 2-, and 1-byte integers:
  integer, parameter :: i8b = selected_int_kind(15)  
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: i2b = selected_int_kind(4)
  integer, parameter :: ib1 = selected_int_kind(2)

  ! symbolic names for kinds of single- and double-precision reals:
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)
  
  ! symbolic names for kinds of single- and double-precision complex:
  integer, parameter :: spc = kind((1.0,1.0))
  integer, parameter :: dpc=kind((1.0d00,1.0d0))

  ! symbolic name for kind type of default logical:
  integer, parameter :: lgt=kind(.true.)
  
  ! frequently used mathematical constants (with precision to spare):
  real(sp), parameter :: &
            pi02=1.57079632679489661923132169163975144209858_sp
  real(sp), parameter :: &
            twopi=6.283185307179586476925286766559005768394_sp
  real(sp), parameter :: &
            sqrt2=1.41421356237309504880168872420969807856967_sp
  real(sp), parameter :: &
            euler=0.5772156649015328606065120900824024310422_sp
  real(dp), parameter :: &
            pi_d=3.141592653589793238462643383279502884197_sp
  real(dp), parameter :: &
            pio2_d=1.57079632679489661923132169163975144209858_sp
  real(dp), parameter :: &
            twopi_d=6.283185307179586476925286766559005768394_sp
  real(dp), parameter :: &
            twoopi_d=0.63661977236758134307553505349006_dp
  real(dp), parameter :: &
            fourpi_d=12.56637061435917295385057353311801153679_dp
  real(dp), parameter :: bigg=6.6723e-11_dp  
  real(dp), parameter :: deg2rad = pi_d/180.0_dp
  real(dp), parameter :: rad2deg = 180.0_dp/pi_d
  complex(dpc), parameter :: czero=(0.0_dp,0.0_dp)
  complex(dpc), parameter :: cone=(1.0_dp,0.0_dp)
  complex(dpc), parameter :: ii=(0.0_dp,1.0_dp)

  ! derived data types for sparse matrices, single and double precision
  type sprs2_sp
    integer(i4b) :: n,len
    real(sp), dimension(:), pointer :: val
    integer(i4b), dimension(:), pointer :: irow
    integer(i4b), dimension(:), pointer :: jcol
  end type sprs2_sp

  type sprs2_dp
    integer(i4b) n,len
    real(dp), dimension(:), pointer :: val
    integer(i4b), dimension(:), pointer :: irow
    integer(i4b), dimension(:), pointer :: jcol
  end type sprs2_dp

  ! ================================================================ !

END MODULE nrtype
