#include "cppdefs.h"
module mod_random
  !----------------------------------------------------------------
  ! Random number generator
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_params, only: pi
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: init_rng, normal_random
  !---------------------------------------------
  ! [variable/type definition]
  !===================================================
contains
  !===========================================
  subroutine init_rng()
    integer, allocatable :: seed(:)
    integer :: n

    call random_seed(size=n)
    allocate (seed(n))
    seed = 9872
    call random_seed(put=seed)
    deallocate (seed)

  end subroutine init_rng
  !===========================================
  real(rk) function normal_random() result(r)
    !---------------------------------------------
    ! Get normally distributed random number from
    ! uniform distribution given by 'call random_number'
    !---------------------------------------------
    real(rk) :: r1, r2

    call random_number(r1); r1 = 1 - r1
    call random_number(r2); r2 = 1 - r2

    r = sqrt(-2 * log(r1)) * cos(2 * pi * r2)

    return
  end function normal_random

end module mod_random
