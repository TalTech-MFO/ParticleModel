#include "cppdefs.h"
module mod_utils
  !----------------------------------------------------------------
  ! Useful functions
  !----------------------------------------------------------------
  use mod_precdefs
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: utils
  !---------------------------------------------
  type :: t_constants
    real(rk) :: pi = 4.0_rk * atan(1.0_rk)
    real(rk) :: gravity = 9.81_rk
    real(rk) :: boltzmann = 1.38064852e-23_rk
  end type t_constants
  !---------------------------------------------
  type :: t_utils
    type(t_constants) :: constants
  contains
    procedure :: normal_random
  end type t_utils
  !---------------------------------------------
  type(t_utils) :: utils !< global instance of t_utils (singleton)
  !===================================================
contains
  !===========================================
  real(rk) function normal_random(this) result(r)
    !---------------------------------------------
    ! Get normally distributed random number from
    ! uniform distribution given by 'call random_number'
    !---------------------------------------------
    class(t_utils), intent(in) :: this
    real(rk) :: r1, r2

    call random_number(r1); r1 = 1 - r1
    call random_number(r2); r2 = 1 - r2

    r = sqrt(-2 * log(r1)) * cos(2 * this%constants%pi * r2)

    return
  end function normal_random

end module mod_utils
