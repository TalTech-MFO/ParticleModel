#include "cppdefs.h"
#include "particle.h"
module mod_kernel
  !----------------------------------------------------------------
  ! Kernel module with basic processes
  !----------------------------------------------------------------
  use mod_common
  use mod_statevector, only: t_statearray
  use mod_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_kernel
  !---------------------------------------------
  type :: t_kernel
    private
    logical :: is_3d
    integer :: vertical_diffusion_method
  contains
    procedure :: init
    procedure :: run
    procedure :: advection
    procedure :: diffusion
  end type t_kernel
  !---------------------------------------------
  !===================================================
contains
  !===========================================
  subroutine init(this)
    class(t_kernel), intent(inout) :: this

    this%is_3d = .false. !! Set this

    return
  end subroutine init
  !===========================================
  function run(this, sa, fieldset, time, dt) result(res)
    class(t_kernel), intent(in) :: this
    type(t_statearray), intent(in) :: sa
    class(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk), intent(in) :: dt
    real(rk) :: res(size(sa%current))

    res = ZERO

    ! TODO: Processes could have an interface with their name, so that we don't have to call run() on each of them
    res = this%advection(sa, fieldset, time, dt) + this%diffusion(sa, fieldset, time, dt)

    return
  end function run
  !===========================================
  function advection(this, sa, fieldset, time, dt)
    class(t_kernel), intent(in) :: this
    type(t_statearray), intent(in) :: sa
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk), intent(in) :: dt
    real(rk) :: advection(sa%nvars)
    real(rk) :: u, v !< u and v velocities in m/s
    real(rk) :: ug, vg !< u and v velocities in degrees/s

    advection = ZERO

    ! TODO: Preprocessor definition ST_SUSPENDED to (global) simulation parameter (e.g., sim_global%particle_state%suspended)
    if (sa%state /= ST_SUSPENDED) return

    u = fieldset%get("U", time, sa%current(1), sa%current(2), sa%current(3))
    v = fieldset%get("V", time, sa%current(1), sa%current(2), sa%current(3))

    call fieldset%domain%xy2lonlat(u, v, sa%current(2), ug, vg)

    advection(1) = ug
    advection(4) = u
    advection(2) = vg
    advection(5) = v

    ! TODO: this%is_3d to global setting (sim_global%run_defs%is_3d or something...)
    if (this%is_3d) then
      advection(3) = fieldset%get("W", time, sa%current(1), sa%current(2), sa%current(3))
      advection(6) = advection(3)
    end if

    return
  end function advection
  !===========================================
  function diffusion(this, sa, fieldset, time, dt)
    class(t_kernel), intent(in) :: this
    type(t_statearray), intent(in) :: sa
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk), intent(in) :: dt
    real(rk) :: diffusion(sa%nvars)
    real(rk) :: u, v !< u and v velocities in m/s
    real(rk) :: ug, vg !< u and v velocities in degrees/s
    real(rk) :: u_rand, v_rand, w_rand !< Random numbers for diffusion
    real(rk) :: ah, kv !< Horizontal and vertical eddy diffusivities

    diffusion = ZERO

    if (sa%state /= ST_SUSPENDED) return

    u_rand = utils%normal_random()
    v_rand = utils%normal_random()
    ah = this%eddy_viscosity(sa, fieldset, time)

    u = u_rand * sqrt(2.0_rk * ah / dt)
    v = v_rand * sqrt(2.0_rk * ah / dt)

    call fieldset%domain%xy2lonlat(u, v, sa%current(2), ug, vg)

    diffusion(1) = ug
    diffusion(2) = vg

    if (this%is_3d) then
      ! TODO: this%vertical_diffusion_method to global setting (sim_global%run_defs%vertical_diffusion_method or something...)
      select case (this%vertical_diffusion_method)
      case (DIFF_VARIABLE)
        kv = fieldset%get("KV", time, sa%current(1), sa%current(2), sa%current(3))
      case (DIFF_DEFAULT)
        ! TODO: diffusion_vert_const to global setting (sim_global%sim_const%diffusion_vert_const or something...)
        kv = diffusion_vert_const
      end select
      w_rand = utils%normal_random()
      diffusion(3) = w_rand * sqrt(2.0_rk * kv / dt)
    end if

    return
  end function diffusion
  !===========================================
  function eddy_viscosity(this, sa, fieldset, time)
    class(t_kernel), intent(in) :: this
    type(t_statearray), intent(in) :: sa
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk) :: eddy_viscosity
    integer :: i, j, k
    real(rk) :: dlon, dlat

    call fieldset%search_indices(t=time, x=sa%current(1), y=sa%current(2), z=sa%current(3), i=i, j=j, k=k)
    dlon = fieldset%domain%dlon%get(i, j)
    dlat = fieldset%domain%dlat%get(i, j)


    return
  end function eddy_viscosity

end module mod_kernel
