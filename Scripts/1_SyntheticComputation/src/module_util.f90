MODULE module_util
  
  INTERFACE my_dot
    MODULE PROCEDURE my_dot_c,my_dot_z
  END INTERFACE 

CONTAINS

  ! ================================================================ !

  FUNCTION my_dot_c(x,y)
    implicit none
    complex :: my_dot_c
    complex, dimension(:), intent(in) :: x,y
    integer :: i,n
      
    my_dot_c = 0.0
    n = size(x,1)
    do i = 1,n
      my_dot_c = my_dot_c + x(i)*y(i)
    end do
      
    return
    
  END FUNCTION my_dot_c

  ! ================================================================ !

  FUNCTION my_dot_z(x,y)
    implicit none
    complex*16 :: my_dot_z
    complex*16, dimension(:), intent(in) :: x,y
    integer :: i,n
      
    my_dot_z = dble(0.0)
    n = size(x,1)
    do i = 1,n
      my_dot_z = my_dot_z + x(i)*y(i)
    end do
 
    return

  END FUNCTION my_dot_z      

  ! ================================================================ !

END MODULE module_util

