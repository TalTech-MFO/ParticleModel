#include "cppdefs.h"
#include "particle.h"
#include "field.h"
module mod_diffusion
  use mod_precdefs
  use mod_errors
  use mod_params
  use time_vars, only: dt
  use field_vars, only: vertical_diffusion_method
  use mod_particle, only: t_particle
  use mod_fieldset, only: t_fieldset
  use mod_physics, only: Ah_Smagorinsky
  use mod_random, only: normal_random
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: diffuse
  !---------------------------------------------
  !===================================================
contains
!===========================================
  subroutine diffuse_2D(p, fieldset, time)

    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: Ah
    real(rk)                        :: i, j
    real(rk)                        :: x0, x1, &
                                       y0, y1

    i = p%ir0
    j = p%jr0

    call fieldset%domain%lonlat2xy(p%lon1, p%lat1, x0, y0)

#if defined(SMAGORINSKY_INTERP_UV)
    Ah = max(Ah_Smagorinsky(fieldset, time, i, j), diffusion_hor_const)
#elif defined(SMAGORINSKY_FULL_FIELD)
    Ah = get_Ah_Smagorinsky_full_field(fieldset, i, j)
#endif

    x1 = x0 + normal_random() * sqrt(2 * Ah * dt)
    y1 = y0 + normal_random() * sqrt(2 * Ah * dt)

    call fieldset%domain%xy2lonlat(x1, y1, p%lon1, p%lat1)
    call fieldset%search_indices(x=p%lon1, y=p%lat1, i=p%i1, j=p%j1, ir=p%ir1, jr=p%jr1)

    return
  end subroutine diffuse_2D
  !===========================================
  subroutine diffuse_3D(p, fieldset, time)

    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    real(rk)                        :: Ah, kv
    real(rk)                        :: i, j, k
    real(rk)                        :: x0, x1, &
                                       y0, y1, &
                                       z0, z1

    dbghead(diffusion :: diffuse_3D)
    i = p%ir0
    j = p%jr0
    k = p%kr0

    debug(i); debug(j); debug(k); debug(time)

    call fieldset%domain%lonlat2xy(p%lon1, p%lat1, x0, y0)
    z0 = p%depth1

    debug(p%lon1); debug(p%lat1); debug(p%depth1)
    debug(x0); debug(y0); debug(z0)

#ifdef DEBUG
    if (k > fieldset%nz) then
      ERROR, "diffusion :: diffuse_3D ", "k >= fieldset%nz"
      ERROR, "diffusion :: diffuse_3D ", "k = ", k
      ERROR, i, j, k
      ERROR, p%lon0, p%lat0, p%depth0
      ERROR, p%lon1, p%lat1, p%depth1
      ERROR, x0, y0, z0
      call throw_error("diffusion :: diffuse_3D", "k >= fieldset%nz")
    end if
#endif

    ! Calculates the Smagorinsky diffusion coefficient for the horizontal
    ! mixing. If the diffusion coefficient is smaller than the minimum
    ! diffusion coefficient, the minimum diffusion coefficient is used
    ! (diffusion_hor_const set in namelist).
    Ah = max(Ah_Smagorinsky(fieldset, time, i, j, k), diffusion_hor_const)

    select case (vertical_diffusion_method)
    case (DIFF_VARIABLE)
      kv = fieldset%get("KV", time, i, j, k)
    case (DIFF_DEFAULT)
      kv = diffusion_vert_const
    case default
      call throw_error("diffusion :: diffuse_3D", "Unknown vertical diffusion method")
    end select

    debug(Ah); debug(kv)

    x1 = x0 + normal_random() * sqrt(2 * Ah * dt)
    y1 = y0 + normal_random() * sqrt(2 * Ah * dt)
#ifndef NO_DIFFUSE_VERTICAL
    z1 = z0 + normal_random() * sqrt(2 * kv * dt)
#else
    z1 = z0
#endif

    debug(x1); debug(y1); debug(z1)

#ifdef DEBUG
    if ((z1 < -150._rk) .or. (z1 > 15._rk)) then
      ERROR, "diffusion :: diffuse_3D", "z1 = ", z1
      ERROR, "diffusion :: diffuse_3D", "z0 = ", z0
      ERROR, "diffusion :: diffuse_3D", "kv = ", kv
      ERROR, "diffusion :: diffuse_3D", "Ah = ", Ah
      ERROR, "diffusion :: diffuse_3D", "vertical_diffusion_method = ", vertical_diffusion_method
      call throw_error("diffusion :: diffuse_3D", "z1 < -150._rk .or. z1 > 15._rk")
    end if
#endif

    call fieldset%domain%xy2lonlat(x1, y1, p%lon1, p%lat1)
    call fieldset%search_indices(t=time, x=p%lon1, y=p%lat1, z=z1, i=p%i1, j=p%j1, k=p%k1, ir=p%ir1, jr=p%jr1, kr=p%kr1)
    p%depth1 = z1

    dbgtail(diffusion :: diffuse_3D)
    return
  end subroutine diffuse_3D
  !===========================================
  subroutine diffuse(p, fieldset, time, dif_3d)

    type(t_particle), intent(inout) :: p
    type(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)            :: time
    logical, intent(in)             :: dif_3d

    if (p%state /= ST_SUSPENDED) return

    select case (dif_3d)
    case (.true.)
      call diffuse_3D(p, fieldset, time)
    case (.false.)
      call diffuse_2D(p, fieldset, time)
    end select

    return
  end subroutine diffuse
end module mod_diffusion
