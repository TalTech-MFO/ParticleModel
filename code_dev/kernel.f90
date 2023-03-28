#include "cppdefs.h"
module mod_kernel
  !----------------------------------------------------------------
  ! [module description]
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use mod_statevector
  use mod_fieldset
  use mod_process_list
  ! use mod_advection
  ! use mod_diffusion ! And so on
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_kernel
  !---------------------------------------------
  type :: t_kernel
    private
    type(t_process_list) :: processes
  contains
    procedure :: init
    procedure :: run
  end type t_kernel
  !---------------------------------------------
  !===================================================
contains
  !===========================================
  subroutine init(this)
    class(t_kernel), intent(inout) :: this

    call this%advection%init()
    call this%diffusion%init()

    return
  end subroutine init
  !===========================================
  function run(this, sv, fieldset, time, dt) result(res)
    class(t_kernel), intent(in) :: this
    real(rk), intent(in) :: sv(:)
    class(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk), intent(in) :: dt
    real(rk) :: res(size(sv))

    res = ZERO

    ! TODO: Processes could have an interface with their name, so that we don't have to call run() on each of them
    res = this%advection%run(sv, fieldset, time, dt) + this%diffusion%run(sv, fieldset, time, dt)

    return
  end function run

end module mod_kernel
