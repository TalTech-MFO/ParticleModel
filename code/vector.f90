#include "cppdefs.h"
module mod_vector
  !----------------------------------------------------------------
  ! Vector class. Simply defines a 3D vector.
  ! This could be replaced with an external library if more functionality is needed.
  !----------------------------------------------------------------
  use mod_common
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_vector
  !---------------------------------------------
  type t_vector
    real(rk) :: x, y, z
  end type t_vector
  !===================================================
contains
  !===========================================
  ! [subroutine/function definition]

end module mod_vector
