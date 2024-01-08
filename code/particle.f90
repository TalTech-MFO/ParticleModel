#include "cppdefs.h"
#include "particle.h"
#include "field.h"
module mod_particle
  !----------------------------------------------------------------
  ! This is the particle type definition
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  ! use mod_domain_vars, only: x0, y0, dx, dy, seamask, depdata, nx, ny
  ! Pass loop vars into functions rather than import?
  use mod_domain, only: t_domain
  use mod_fieldset, only: t_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_particle
  !---------------------------------------------
  ! Particle type
  type t_particle
    !---------------------------------------------
    logical  :: is_active = .true.          ! Skip particle in loop if is_active == .false.
    logical  :: kill_beached, kill_boundary ! Set is_active=.false. if beached or on boundary?
    integer  :: warnings = 0
    integer  :: state = ST_SUSPENDED        ! States: active, beached, on boundary, bottom. Enumerated in cppdefs.h
    real(rk) :: id                          ! Origin of particle, number
    logical  :: restarted = .false.         ! Particle was restarted from file
    !---------------------------------------------
    ! Indices
    integer  :: i0, j0, k0                  ! Position (grid cell indices, original)
    real(rk) :: ir0, jr0, kr0               ! Position (real indices, original)
    integer  :: i1, j1, k1                  ! Position (grid cell indices, t + dt)
    real(rk) :: ir1, jr1, kr1               ! Position (real indices, t + dt)
    !---------------------------------------------
    ! Coordinates
    real(rk) :: lon0 = ZERO                 ! Position (original)
    real(rk) :: lat0 = ZERO                 ! Position (original)
    real(rk) :: depth0 = ZERO               ! Position (original)
    real(rk) :: lon1 = ZERO                 ! Position (t + dt)
    real(rk) :: lat1 = ZERO                 ! Position (t + dt)
    real(rk) :: depth1 = ZERO               ! Position (t + dt)
    !---------------------------------------------
    ! Velocity
    real(rk) :: u0 = ZERO                   ! Velocity (original)
    real(rk) :: v0 = ZERO                   ! Velocity (original)
    real(rk) :: w0 = ZERO                   ! Velocity (original)
    real(rk) :: u1 = ZERO                   ! Velocity (t + dt)
    real(rk) :: v1 = ZERO                   ! Velocity (t + dt)
    real(rk) :: w1 = ZERO                   ! Velocity (t + dt)
    real(rk) :: vel_vertical = ZERO         ! Settling velocity (Kooi)
    real(rk) :: u_diff = ZERO               ! Diffusion velocity
    real(rk) :: v_diff = ZERO               ! Diffusion velocity
    real(rk) :: w_diff = ZERO               ! Diffusion velocity
    !---------------------------------------------
    real(rk) :: rho = ZERO                  ! Density
    real(rk) :: rho0 = ZERO                 ! Initial density
    real(rk) :: radius = ZERO               ! Radius
    real(rk) :: radius0 = ZERO              ! Initial radius
    real(rk) :: age = ZERO                  ! Age
    real(rk) :: max_age = ZERO              ! Maximum age
    real(rk) :: traj_len = ZERO             ! Particle trajectory length
    real(rk) :: time_on_beach = ZERO        ! Time spent in the beach area
    real(rk) :: beaching_time               ! Different particles may essentialy have different beaching times
    !---------------------------------------------
    ! Environment variables
    real(rk) :: delta_rho = ZERO            ! Density difference between particle and surrounding water
    real(rk) :: kin_visc = ZERO              ! Kinematic viscosity surrounding particle
    real(rk) :: u_star = ZERO               ! Friction velocity surrounding particle
    !---------------------------------------------
    ! Biofouling variables
    real(rk) :: growth_biofilm = ZERO       ! Attached algal growth
    real(rk) :: h_biofilm = ZERO            ! Thickness of biofilm

  contains
    private
    procedure, public :: update
    procedure         :: check_boundaries
    procedure         :: check_age
    procedure, public :: check_depth
    procedure, public :: print_info
    procedure, public :: get_variable
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

    p%id = id
    p%lon0 = lon
    p%lat0 = lat
    p%depth0 = depth
    p%lon1 = lon
    p%lat1 = lat
    p%depth1 = depth
    p%is_active = is_active

    if (.not. p%is_active) then
      p%state = ST_INIT_ERROR
      ! No more information is needed if the particle is not active from the start
      return
    end if

    p%beaching_time = beaching_time
    p%rho = rho
    p%rho0 = rho
    p%radius = radius
    p%radius0 = radius
    p%max_age = max_age
    p%kill_beached = kill_beached
    p%kill_boundary = kill_boundary

    call fieldset%search_indices(time, lon, lat, depth, i=p%i0, j=p%j0, k=p%k0, ir=p%ir0, jr=p%jr0, kr=p%kr0)
    p%i1 = p%i0
    p%j1 = p%j0
    p%k1 = p%k0
    p%ir1 = p%ir0
    p%jr1 = p%jr0
    p%kr1 = p%kr0

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

    p%lon0 = lon
    p%lat0 = lat
    p%depth0 = depth

    p%lon1 = lon
    p%lat1 = lat
    p%depth1 = depth

    p%i0 = i0
    p%j0 = j0
    p%k0 = k0
    p%ir0 = ir0
    p%jr0 = jr0
    p%kr0 = kr0

    p%i1 = i0
    p%j1 = j0
    p%k1 = k0
    p%ir1 = ir0
    p%jr1 = jr0
    p%kr1 = kr0

    p%u0 = u0
    p%v0 = v0
    p%w0 = w0

    p%u1 = u0
    p%v1 = v0
    p%w1 = w0

    p%vel_vertical = vel_vertical

    p%traj_len = traj_len

    p%id = id
    p%is_active = is_active
    p%state = state

    p%beaching_time = beaching_time
    p%time_on_beach = time_on_beach

    p%rho = rho
    p%rho0 = rho0

    p%radius = radius
    p%radius0 = radius0
    p%h_biofilm = h_biofilm

    p%age = age
    p%max_age = max_age

    p%kill_beached = kill_beached
    p%kill_boundary = kill_boundary

    p%restarted = .true.

  end function ctor_particle_restart
  !===========================================
  real(rk) function get_variable(this, name) result(res)
    class(t_particle), intent(in) :: this
    character(len=LEN_CHAR_S), intent(in) :: name

    select case (name)
    case ("u")
      res = this%u0
    case ("v")
      res = this%v0
    case ("w")
      res = this%w0
    case ("u_diff")
      res = this%u_diff
    case ("v_diff")
      res = this%v_diff
    case ("w_diff")
      res = this%w_diff
    case ("vel_vertical")
      res = this%vel_vertical
    case ("rho")
      res = this%rho
    case ("radius")
      res = this%radius
    case ("age")
      res = this%age
    case ("traj_len")
      res = this%traj_len
    case ("time_on_beach")
      res = this%time_on_beach
    case ("h_biofilm")
      res = this%h_biofilm
    case ("growth_biofilm")
      res = this%growth_biofilm
    case ("volume")
      res = this%volume()
    case ("surface_area")
      res = this%surface_area()
    end select

  end function get_variable
  !===========================================
  elemental real(rk) function volume(this)
    class(t_particle), intent(in) :: this
    real(rk), parameter :: pi = 4.*atan(1.)

    volume = 4./3.*pi * this%radius0**3

  end function volume
  !===========================================
  elemental real(rk) function surface_area(this)
    class(t_particle), intent(in) :: this
    real(rk), parameter :: pi = 4.*atan(1.)

    surface_area = 4.*pi * this%radius0**2.

  end function surface_area
  !===========================================
  subroutine update(this, fieldset, time, dt)

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)             :: time
    real(rk), intent(in)             :: dt
    real(rk)                         :: x0, y0
    real(rk)                         :: x1, y1

    this%age = this%age + dt
    call this%check_age()

    ! Only update the age if the particle is beached or otherwise not active (but still alive)
    if (this%state < ST_SUSPENDED) return

    call this%check_boundaries(fieldset, time, dt)

    call fieldset%domain%lonlat2xy(this%lon0, this%lat0, x0, y0)
    call fieldset%domain%lonlat2xy(this%lon1, this%lat1, x1, y1)
    this%traj_len = this%traj_len + &
                    sqrt((x1 - x0)**2 + &
                         (y1 - y0)**2 + &
                         (this%depth1 - this%depth0)**2) ! This will always be 0 if run_3d=.false.

    this%i0 = this%i1
    this%j0 = this%j1
    this%k0 = this%k1
    this%ir0 = this%ir1
    this%jr0 = this%jr1
    this%kr0 = this%kr1

    this%lon0 = this%lon1
    this%lat0 = this%lat1
    this%depth0 = this%depth1 ! This will always stay the same if run_3d=.false.

    this%u0 = this%u1
    this%v0 = this%v1
    this%w0 = this%w1

    return
  end subroutine update
  !===========================================
  subroutine check_age(this)

    class(t_particle), intent(inout) :: this

    if ((this%max_age > ZERO) .and. (this%age > this%max_age)) then
      this%is_active = .false.
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

    dbghead(particle :: check_depth)

    debug(this%ir1)
    debug(this%jr1)

#ifdef PARTICLE_SNAP_SEALVL
    elev = fieldset%sealevel(t, this%ir1, this%jr1)
    if (this%depth1 >= elev) then
      this%depth1 = elev
      dbgtail(particle :: check_depth)
      return
    end if
#endif

    ! dep = fieldset%domain%get_bathymetry(this%i1, this%j1)
    dep = -ONE * fieldset%domain%get_bathymetry(this%i1, this%j1)

    ! The particle is past the bottom
    ! if (this%depth1 <= -1.0 * dep) then
    if (this%depth1 <= dep) then
      ! ERROR, "Particle is settling:"
      ! ERROR, "  dep    = ", dep
      ! ERROR, "  lon0   = ", this%lon0
      ! ERROR, "  lat0   = ", this%lat0
      ! ERROR, "  depth0 = ", this%depth0
      ! ERROR, "  i0     = ", this%i0
      ! ERROR, "  j0     = ", this%j0
      ! ERROR, "  k0     = ", this%k0
      ! ERROR, "  ir0    = ", this%ir0
      ! ERROR, "  jr0    = ", this%jr0
      ! ERROR, "  kr0    = ", this%kr0
      ! ERROR, "  lon1   = ", this%lon1
      ! ERROR, "  lat1   = ", this%lat1
      ! ERROR, "  depth1 = ", this%depth1
      ! ERROR, "  i1     = ", this%i1
      ! ERROR, "  j1     = ", this%j1
      ! ERROR, "  k1     = ", this%k1
      ! ERROR, "  ir1    = ", this%ir1
      ! ERROR, "  jr1    = ", this%jr1
      ! ERROR, "  kr1    = ", this%kr1
      ! ERROR, "  u1     = ", this%u1
      ! ERROR, "  v1     = ", this%v1
      ! ERROR, "  w1     = ", this%w1
      ! ERROR, "  du     = ", this%u_diff
      ! ERROR, "  dv     = ", this%v_diff
      ! ERROR, "  dw     = ", this%w_diff
      ! ERROR, "  vs     = ", this%vel_vertical
      ! ERROR, "  state  = ", this%state
      ! call throw_error("particle :: check_depth", "Particle settled.")

      this%state = ST_BOTTOM
      this%depth1 = dep
      this%w1 = ZERO
      dbgtail(particle :: check_depth)
      return
    end if

    ! Reset to SUSPENDED if resuspended
    if ((this%depth1 > dep) .and. (this%state == ST_BOTTOM)) then
      this%state = ST_SUSPENDED
      dbgtail(particle :: check_depth)
      return
    end if

    dbgtail(particle :: check_depth)
    return
  end subroutine check_depth
  !===========================================
  subroutine check_boundaries(this, fieldset, time, dt)

    class(t_particle), intent(inout) :: this
    class(t_fieldset), intent(in)    :: fieldset
    real(rk), intent(in)             :: time
    real(rk), intent(in)             :: dt
    integer                          :: i, j
    integer                          :: seamask_val

    i = this%i1
    j = this%j1

    if (i < 1 .or. i > fieldset%domain%nx .or. &
        j < 1 .or. j > fieldset%domain%ny) then
      this%is_active = .false.
      this%state = ST_BOUNDARY
      return
    end if

    seamask_val = fieldset%domain%get_seamask(i=i, j=j)

    select case (seamask_val)
    case (DOM_BEACH)
      this%time_on_beach = this%time_on_beach + dt
      !---------------------------------------------
      ! Change state if beaching time exceeded or on boundary
#ifndef PARTICLE_BEACH_IMMEDIATELY
      if (this%time_on_beach >= this%beaching_time) then
#endif
        if (this%kill_beached) this%is_active = .false.
        this%state = ST_BEACHED
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

    call this%check_depth(fieldset, time)

    return
  end subroutine check_boundaries
  !===========================================
  subroutine print_info(this)

    class(t_particle), intent(in) :: this

    FMT1, "Indices (i, j, k)"
    FMT2, this%i0, this%j0, this%k0
    FMT2, this%i1, this%j1, this%k1

    FMT1, "Real indices (ir, jr, kr)"
    FMT2, this%ir0, this%jr0, this%kr0
    FMT2, this%ir1, this%jr1, this%kr1

    FMT1, "Position (lon, lat, depth)"
    FMT2, this%lon0, this%lat0, this%depth0
    FMT2, this%lon1, this%lat1, this%depth1

    FMT1, "Velocity (u, v, w)"
    FMT2, this%u0, this%v0, this%w0
    FMT2, this%u1, this%v1, this%w1

    FMT1, "Velocity from buoyancy"
    FMT2, this%vel_vertical

    FMT1, "State"
    FMT2, "state: ", this%state
    FMT2, "is_active: ", this%is_active
    FMT2, "restarted: ", this%restarted

    FMT1, "Other characteristics"
    FMT2, "Age: ", this%age
    FMT2, "Distance: ", this%traj_len
    FMT2, "Radius (R, R0): ", this%radius, this%radius0
    FMT2, "Density (rho, rho0): ", this%rho, this%rho0
    FMT2, "Time on beach: ", this%time_on_beach
    FMT2, "Beaching time: ", this%beaching_time
    FMT2, "Biofilm thickness: ", this%h_biofilm

  end subroutine print_info

end module mod_particle
