#include "cppdefs.h"
#include "advection.h"
module mod_advection
  use mod_common
  use mod_process
  use mod_fieldset, only: t_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_advection
  !---------------------------------------------
  type, extends(t_process) :: t_advection
    private
    logical :: is_3d
  contains
    procedure :: init => init_advection
    procedure :: run => run_advection
  end type t_advection
  !---------------------------------------------
  !===================================================
contains
  !===========================================
  subroutine init_advection(this)
    class(t_advection), intent(inout) :: this
  end subroutine init_advection
  !===========================================
  function run_advection(this, sv, fieldset, time, dt) result(res)
    class(t_advection), intent(in) :: this
    type(t_fieldset), intent(in)   :: fieldset
    real(rk), intent(in)           :: time, dt
    real(rk), intent(in)           :: sv(:)
    real(rk) :: res(size(sv))
    real(rk) :: u, v, ug, vg

    res = ZERO

    u = fieldset%get("U", time, sv(1), sv(2), sv(3))
    v = fieldset%get("V", time, sv(1), sv(2), sv(3))

    call fieldset%domain%xy2lonlat(u, v, sv(2), ug, vg)

    res(1) = ug
    res(4) = u
    res(2) = vg
    res(5) = v

    if (this%is_3d) then
      res(3) = fieldset%get("W", time, sv(1), sv(2), sv(3))
      res(6) = res(3)
    end if

    return
  end function run_advection
end module mod_advection
