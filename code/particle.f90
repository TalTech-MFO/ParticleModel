#include "cppdefs.h"
#include "particle.h"
#include "field.h"
#include "file.h"
#include "ompdefs.h"
module mod_particle
  !----------------------------------------------------------------
  ! This is the particle type definition
  !----------------------------------------------------------------
  use mod_common
  ! use mod_domain_vars, only: x0, y0, dx, dy, seamask, depdata, nx, ny
  ! Pass loop vars into functions rather than import?
  use time_vars, only: dt
  use mod_params, only: run_3d, pi
  use mod_vector, only: t_vector
  use mod_domain, only: t_domain
  use mod_fieldset, only: t_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_particle
  !---------------------------------------------
  ! Particle state type, holds the current parameters of the particle
  type t_particle_state
    logical        :: is_active = .true.          !< Skip particle in loop if is_active == .false.
    type(t_vector) :: pos                         !< Position (lon, lat, depth)
    type(t_vector) :: vel                         !< Velocity (u, v, w)
    integer        :: state = ST_SUSPENDED        !< States: active, beached, on boundary, bottom. Enumerated in cppdefs.h
    real(rk)       :: age = ZERO                  !< Age
    real(rk)       :: radius = ZERO               !< Radius
    real(rk)       :: density = ZERO              !< Density
    real(rk)       :: time_on_beach = ZERO        !< Time spent on the beach
    real(rk)       :: trajectory = ZERO           !< Trajectory length
    ! Biofouling variables
    real(rk)       :: growth_biofilm = ZERO       !< Attached algal growth
    real(rk)       :: h_biofilm = ZERO            !< Thickness of biofilm
  end type t_particle_state
  !---------------------------------------------
  ! Particle parameters type, holds the parameters that are constant for the particle
  type t_particle_params
    real(rk) :: id                    !< Origin of particle, number
    real(rk) :: density = ZERO        !< Initial density
    real(rk) :: radius = ZERO         !< Initial radius
    real(rk) :: max_age = ZERO        !< Maximum age
    real(rk) :: beaching_time = ZERO  !< Time to beach
    logical  :: kill_beached          !< Set is_active=.false. if beached or on boundary?
    logical  :: kill_boundary         !< Set is_active=.false. if beached or on boundary?
    logical  :: restarted = .false.   !< Particle was restarted from file
  end type t_particle_params
  !---------------------------------------------
  ! Particle type
  type t_particle
    type(t_particle_state)  :: st     !< Particle state
    type(t_particle_params) :: par    !< Particle parameters
  contains
    private
    procedure, public :: update
    procedure, public :: get_state_array
    procedure, public :: set_state_array
    procedure         :: check_boundaries
    procedure         :: check_age
    procedure, public :: check_depth
    procedure, public :: print_info
    procedure, public :: get_nvars
    procedure, public :: volume
    procedure, public :: surface_area
  end type t_particle

  interface t_particle
    module procedure :: ctor_particle
    module procedure :: ctor_particle_restart
  end interface t_particle

  !===================================================
contains
  !===========================================
  type(t_particle) function ctor_particle(lon, lat, depth, &
                                          id, beaching_time, &
                                          rho, radius, max_age, &
                                          kill_beached, kill_boundary, &
                                          is_active, &
                                          fieldset, time) result(p)
    real(rk), intent(in)         :: lon, lat, depth
    real(rk), intent(in)         :: id
    real(rk), intent(in)         :: beaching_time
    real(rk), intent(in)         :: rho, radius
    real(rk), intent(in)         :: max_age
    logical, intent(in)          :: kill_beached, kill_boundary
    logical, intent(in)          :: is_active
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in)         :: time

    p%par%id = id
    p%st%pos%x = lon
    p%st%pos%y = lat
    p%st%pos%z = depth
    p%st%is_active = is_active

    if (.not. p%st%is_active) then
      p%st%state = ST_INIT_ERROR
      ! No more information is needed if the particle is not active from the start
      return
    end if

    p%par%beaching_time = beaching_time
    p%par%density = rho
    p%st%density = rho
    p%par%radius = radius
    p%st%radius = radius
    p%par%max_age = max_age
    p%par%kill_beached = kill_beached
    p%par%kill_boundary = kill_boundary

    ! call fieldset%search_indices(time, lon, lat, depth, i=p%i0, j=p%j0, k=p%k0, ir=p%ir0, jr=p%jr0, kr=p%kr0)

  end function ctor_particle
  !===========================================
  type(t_particle) function ctor_particle_restart(lon, lat, depth, &
                                                  i0, j0, k0, &
                                                  ir0, jr0, kr0, &
                                                  id, beaching_time, &
                                                  rho, rho0, &
                                                  radius, radius0, &
                                                  h_biofilm, &
                                                  age, max_age, kill_beached, kill_boundary, &
                                                  u0, v0, w0, vel_vertical, &
                                                  traj_len, time_on_beach, is_active, state) result(p)
    real(rk), intent(in)         :: lon, lat, depth
    real(rk), intent(in)         :: id
    real(rk), intent(in)         :: beaching_time
    real(rk), intent(in)         :: rho, rho0
    real(rk), intent(in)         :: radius, radius0
    real(rk), intent(in)         :: h_biofilm
    real(rk), intent(in)         :: age, max_age
    real(rk), intent(in)         :: ir0, jr0, kr0
    real(rk), intent(in)         :: u0, v0, w0, vel_vertical
    real(rk), intent(in)         :: traj_len, time_on_beach
    integer, intent(in)          :: i0, j0, k0
    integer, intent(in)          :: state
    logical, intent(in)          :: kill_beached, kill_boundary, is_active

    p%st%pos%x = lon
    p%st%pos%y = lat
    p%st%pos%z = depth

    p%st%vel%x = u0
    p%st%vel%y = v0
    p%st%vel%z = w0

    p%st%trajectory = traj_len

    p%par%id = id
    p%st%is_active = is_active
    p%st%state = state

    p%par%beaching_time = beaching_time
    p%st%time_on_beach = time_on_beach

    p%par%density = rho
    p%st%density = rho

    p%par%radius = radius
    p%st%radius = radius

    p%st%h_biofilm = h_biofilm

    p%st%age = age
    p%par%max_age = max_age

    p%par%kill_beached = kill_beached
    p%par%kill_boundary = kill_boundary

    p%par%restarted = .true.

  end function ctor_particle_restart
  !===========================================
  pure real(rk) function volume(this)
    class(t_particle), intent(in) :: this

    volume = 4./3.*pi * this%par%radius**3

  end function volume
  !===========================================
  pure real(rk) function surface_area(this)
    class(t_particle), intent(in) :: this

    surface_area = 4.*pi * this%par%radius**2.

  end function surface_area
  !===========================================
  pure integer function get_nvars(this)
    class(t_particle), intent(in) :: this

    get_nvars = 10

  end function get_nvars
  !===========================================
  function get_state_array(this) result(res)
    class(t_particle), intent(in) :: this
    ! real(rk), allocatable :: res(:)
    real(rk) :: res(this%get_nvars())

    ! allocate (res(this%get_nvars()))
    res(1) = this%st%pos%x
    res(2) = this%st%pos%y
    res(3) = this%st%pos%z
    res(4) = this%st%vel%x
    res(5) = this%st%vel%y
    res(6) = this%st%vel%z
    res(7) = this%st%density
    res(8) = this%st%radius
    res(9) = this%st%h_biofilm
    res(10) = this%st%age

    return
  end function get_state_array
  !===========================================
  subroutine set_state_array(this, sa)
    class(t_particle), intent(inout) :: this
    real(rk), intent(in) :: sa(this%get_nvars())

    this%st%pos%x = sa(1)
    this%st%pos%y = sa(2)
    this%st%pos%z = sa(3)
    this%st%vel%x = sa(4)
    this%st%vel%y = sa(5)
    this%st%vel%z = sa(6)
    this%st%density = sa(7)
    this%st%radius = sa(8)
    this%st%h_biofilm = sa(9)
    this%st%age = sa(10)

    return
  end subroutine set_state_array
  !===========================================
  subroutine update(this, fieldset, time)

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)             :: time
    real(rk)                         :: x0, y0
    real(rk)                         :: x1, y1

    ! this%age = this%age + dt
    call this%check_age()

    ! Only update the age if the particle is beached or otherwise not active (but still alive)
    if (this%st%state < ST_SUSPENDED) return

    call this%check_boundaries(fieldset, time)

    ! ! Trajectory length should be calculated in the kernel
    ! call fieldset%domain%lonlat2xy(this%lon0, this%lat0, x0, y0)
    ! call fieldset%domain%lonlat2xy(this%lon1, this%lat1, x1, y1)
    ! this%traj_len = this%traj_len + &
    !                 sqrt((x1 - x0)**2 + &
    !                      (y1 - y0)**2 + &
    !                      (this%depth1 - this%depth0)**2) ! This will always be 0 if run_3d=.false.

    ! this%i0 = this%i1
    ! this%j0 = this%j1
    ! this%k0 = this%k1
    ! this%ir0 = this%ir1
    ! this%jr0 = this%jr1
    ! this%kr0 = this%kr1

    ! this%lon0 = this%lon1
    ! this%lat0 = this%lat1
    ! this%depth0 = this%depth1 ! This will always stay the same if run_3d=.false.

    ! this%u0 = this%u1
    ! this%v0 = this%v1
    ! this%w0 = this%w1

    return
  end subroutine update
  !===========================================
  subroutine check_age(this)

    class(t_particle), intent(inout) :: this

    if ((this%par%max_age > ZERO) .and. (this%st%age > this%par%max_age)) then
      this%st%is_active = .false.
    end if

  end subroutine check_age
  !===========================================
  subroutine check_depth(this, fieldset, t)
    !---------------------------------------------
    ! TODO: Interpolation for bathymetry?
    !---------------------------------------------

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)             :: t
    real(rk)                         :: dep
    real(rk)                         :: elev

#ifdef PARTICLE_SNAP_SEALVL
    elev = fieldset%sealevel(t, this%ir1, this%jr1)
    if (this%depth1 >= elev) then
      this%depth1 = elev
      return
    end if
#endif

    ! ! TODO: Interpolation for bathymetry?
    ! dep = fieldset%domain%get_bathymetry(this%i1, this%j1)

    ! The particle is past the bottom
    if (this%st%pos%z <= -1.0 * dep) then
      this%st%pos%z = -1.0 * dep
      this%st%state = ST_BOTTOM
      this%st%vel%z = ZERO
      return
    end if

    ! Reset to SUSPENDED if resuspended
    if ((this%st%pos%z > -1.0 * dep) .and. (this%st%state == ST_BOTTOM)) then
      this%st%state = ST_SUSPENDED
      return
    end if

  end subroutine check_depth
  !===========================================
  subroutine check_boundaries(this, fieldset, time)

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)             :: time
    integer                          :: i, j
    integer                          :: seamask_val

    i = this%i1
    j = this%j1

    if (i < 1 .or. i > fieldset%domain%nx .or. &
        j < 1 .or. j > fieldset%domain%ny) then
      this%st%is_active = .false.
      this%st%state = ST_BOUNDARY
      return
    end if

    seamask_val = fieldset%domain%get_seamask(i=i, j=j)

    select case (seamask_val)
    case (DOM_BEACH)
      this%st%time_on_beach = this%st%time_on_beach + dt
      !---------------------------------------------
      ! Change state if beaching time exceeded or on boundary
#ifndef PARTICLE_BEACH_IMMEDIATELY
      if (this%st%time_on_beach >= this%par%beaching_time) then
#endif
        if (this%par%kill_beached) this%st%is_active = .false.
        this%st%state = ST_BEACHED
#ifndef PARTICLE_BEACH_IMMEDIATELY
      end if
#endif

    case (DOM_LAND)
! #if defined PARTICLE_BOUNCE ! ! These two methods are deleted and should never be rewritten again.
!       call this%bounce(fieldset)
! #elif defined PARTICLE_REDIRECT
!       call this%redirect(fieldset)
#if defined PARTICLE_BEACH_IMMEDIATELY
      if (this%kill_beached) this%is_active = .false.
      this%state = ST_BEACHED
#else
      this%i1 = this%i0
      this%j1 = this%j0
      this%k1 = this%k0
      this%ir1 = this%ir0
      this%jr1 = this%jr0
      this%kr1 = this%kr0
      this%lon1 = this%lon0
      this%lat1 = this%lat0
      this%depth1 = this%depth0
#endif

      !---------------------------------------------
      ! The bounce can happen only in the beach zone, so add to time on beach
      this%time_on_beach = this%time_on_beach + dt
      !---------------------------------------------
      ! Change state if beaching time exceeded or on boundary
      if (this%time_on_beach >= this%beaching_time) then
        if (this%kill_beached) this%is_active = .false.
        this%state = ST_BEACHED
      end if

    case (DOM_SEA)
      this%time_on_beach = ZERO

    case (DOM_BOUNDARY)
      if (this%kill_boundary) this%is_active = .false.
      this%state = ST_BOUNDARY
    end select

    if (run_3d) call this%check_depth(fieldset, time)

    return
  end subroutine check_boundaries
  !===========================================
  subroutine print_info(this)

    class(t_particle), intent(in) :: this

    ! FMT1, "Indices (i, j, k)"
    ! FMT2, this%i0, this%j0, this%k0
    ! FMT2, this%i1, this%j1, this%k1

    ! FMT1, "Real indices (ir, jr, kr)"
    ! FMT2, this%ir0, this%jr0, this%kr0
    ! FMT2, this%ir1, this%jr1, this%kr1

    FMT1, "Position (lon, lat, depth)"
    FMT2, this%st%pos%x, this%st%pos%y, this%st%pos%z

    FMT1, "Velocity (u, v, w)"
    FMT2, this%st%vel%x, this%st%vel%y, this%st%vel%z

    ! FMT1, "Velocity from buoyancy"
    ! FMT2, this%vel_vertical

    FMT1, "State"
    FMT2, "state: ", this%st%state
    FMT2, "is_active: ", this%st%is_active
    FMT2, "restarted: ", this%par%restarted

    FMT1, "Other characteristics"
    FMT2, "Age: ", this%st%age
    FMT2, "Distance: ", this%st%trajectory
    FMT2, "Radius (R, R0): ", this%st%radius, this%par%radius
    FMT2, "Density (rho, rho0): ", this%st%density, this%par%density
    FMT2, "Time on beach: ", this%st%time_on_beach
    FMT2, "Beaching time: ", this%par%beaching_time
    FMT2, "Biofilm thickness: ", this%st%h_biofilm

  end subroutine print_info

end module mod_particle
