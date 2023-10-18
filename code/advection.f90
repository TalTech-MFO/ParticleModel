#include "cppdefs.h"
#include "particle.h"
module mod_advection
  use mod_common
  use mod_process
  use mod_fieldset, only: t_fieldset
  use mod_statevector, only: t_statearray
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
  function run_advection(this, sa, fieldset, time, dt) result(res)
    class(t_advection), intent(in) :: this
    type(t_fieldset), intent(in)   :: fieldset
    real(rk), intent(in)           :: time, dt
    type(t_statearray), intent(in) :: sa
    real(rk) :: res(sa%nvars)
    real(rk) :: u, v !< u and v velocities in m/s
    real(rk) :: ug, vg !< u and v velocities in degrees/s

    res = ZERO

    ! TODO: Preprocessor definition ST_SUSPENDED to (global) simulation parameter (e.g., sim_global%particle_state%suspended)
    if (sa%state /= ST_SUSPENDED) return

    u = fieldset%get("U", time, sa%current(1), sa%current(2), sa%current(3))
    v = fieldset%get("V", time, sa%current(1), sa%current(2), sa%current(3))

    call fieldset%domain%xy2lonlat(u, v, sa%current(2), ug, vg)

    res(1) = ug
    res(4) = u
    res(2) = vg
    res(5) = v

    if (this%is_3d) then
      res(3) = fieldset%get("W", time, sa%current(1), sa%current(2), sa%current(3))
      res(6) = res(3)
    end if

    return
  end function run_advection
end module mod_advection
