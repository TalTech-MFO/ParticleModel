!===================================================
module mod_particle_vars
  !----------------------------------------------------------------
  ! This module includes variables related to particles:
  ! - number of particles, initial locations or something (maybe)...
  ! - anything else?
  !----------------------------------------------------------------
#ifdef USE_OMP
  use omp_lib
#endif
  use mod_precdefs
  use mod_errors
  use mod_particle
  use nc_manager, only: nc_read_real_1d, nc_read_real_2d, nc_get_dim_len, nc_var_exists
  use mod_domain, only: t_domain
  use time_vars, only: nTimes, run_start_dt, run_end_dt, dt
  use run_params, only: runid, restart, restart_path
  use mod_datetime, only: t_datetime, datetime_from_netcdf
  use mod_fieldset, only: t_fieldset
  implicit none
  private
  !===================================================
  public :: inputstep, particle_init_method, coordfile, &
            max_age, kill_beached, kill_boundary, runparts, n_total_particles
  public :: particles
  public :: initialise_particles, release_particles
  !---------------------------------------------
  logical                       :: kill_beached, kill_boundary ! Set is_active=.false. if beached or on boundary?
  integer                       :: inputstep, &                ! How often are particles released?
                                   particle_init_method, &     ! Read initial positions (1 - txt, 2 - .nc)
                                   n_particles, &              ! Number of particles
                                   n_restart_particles, &
                                   n_total_particles, &        ! Total number of particles (including restart)
                                   n_init_times, &
                                   runparts = 0, &             ! Number of particles to loop over
                                   i_release = 1
  real(rk)                      :: max_age                     ! Lifetime (for all particles) in timesteps
  character(len=LEN_CHAR_L)     :: coordfile                   ! File containing particle locations at init.
  type(t_particle), allocatable :: particles(:)                ! Array of particles
  !---------------------------------------------
  type, private :: t_initial_position
    type(t_datetime)                    :: release_date
    integer                             :: next_idx = 1, &
                                           time_idx = 1, &
                                           n_particles
#ifdef IGNORE_BAD_PARTICLES
    integer                             :: n_good_particles = 0
#endif
    real(rk), allocatable, dimension(:) :: x, y, z, &
                                           rho, radius, &
                                           beaching_time, &
                                           id
    logical, allocatable, dimension(:) :: is_active
  contains
    procedure :: allocate_n_init_particles
    procedure :: check_initial_coordinates
  end type t_initial_position
  !---------------------------------------------
  type(t_initial_position), allocatable :: init_coords(:)
  !---------------------------------------------
  interface initialise_particles
    module procedure init_particles
  end interface initialise_particles
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  subroutine init_particles(fieldset)
    type(t_fieldset), intent(in) :: fieldset

    select case (particle_init_method)
    case (TXT_FILE)
      call init_particles_from_coordfile(fieldset)
    case (NC_FILE)
      call init_particles_from_netcdf(fieldset)
    end select

    if (restart) then
      call check_restart_file
      allocate (particles(n_particles + n_restart_particles))
      FMT2, "Allocated array for", n_particles + n_restart_particles, "particles"
      call read_restart_file
      n_total_particles = n_particles + n_restart_particles
    else
      allocate (particles(n_particles))
      FMT2, "Allocated array for", n_particles, "particles"
      n_total_particles = n_particles
    end if

  end subroutine init_particles
  !===========================================
  subroutine check_restart_file()

    character(len=LEN_CHAR_L) :: restart_filename
    character(len=14) :: time_str
    logical  :: file_exists

    write (time_str, '(i0.14)') run_start_dt%short_format(include_time=.true.)

    restart_filename = trim(restart_path)//"/"//trim(runid)//"."//trim(time_str)//".restart.dat"

    inquire (file=trim(restart_filename), exist=file_exists)
    if (file_exists) then
      open (RESTARTFILE, file=trim(restart_filename), action='read', iostat=ierr)
      read (RESTARTFILE, *) n_restart_particles
      close (RESTARTFILE)
    else
      call throw_error("particle :: read_restart_file", "No restart file found.")
    end if

  end subroutine check_restart_file
  !===========================================
  subroutine read_restart_file()
    !---------------------------------------------
    ! Read restart file (latest position)
    ! Initialise in particles array
    ! Update runparts
    !---------------------------------------------
    character(len=LEN_CHAR_L) :: restart_filename
    character(len=14) :: time_str
    integer  :: i
    real(rk) :: lon, lat, depth
    real(rk) :: id
    real(rk) :: beaching_time
    real(rk) :: rho, rho0
    real(rk) :: radius, radius0
    real(rk) :: h_biofilm
    real(rk) :: age, max_age
    real(rk) :: ir0, jr0, kr0
    real(rk) :: u0, v0, w0, vel_vertical
    real(rk) :: traj_len, time_on_beach
    integer  :: i0, j0, k0
    integer  :: state
    logical  :: kill_beached, kill_boundary, is_active

    write (time_str, '(i0.14)') run_start_dt%short_format(include_time=.true.)

    restart_filename = trim(restart_path)//"/"//trim(runid)//"."//trim(time_str)//".restart.dat"

    open (RESTARTFILE, file=trim(restart_filename), action='read', status='old', iostat=ierr)
    read (RESTARTFILE, *)
    do i = 1, n_restart_particles
      read (RESTARTFILE, *, iostat=ierr) lon, lat, depth, &
        i0, j0, k0, &
        ir0, jr0, kr0, &
        id, beaching_time, &
        rho, rho0, &
        radius, radius0, &
        h_biofilm, &
        age, max_age, kill_beached, kill_boundary, &
        u0, v0, w0, vel_vertical, &
        traj_len, time_on_beach, is_active, state

      if (ierr /= 0) then
        call throw_error("particle :: read_restart_file", "Error reading restart file.")
      end if

      particles(i) = t_particle(lon, lat, depth, &
                                i0, j0, k0, &
                                ir0, jr0, kr0, &
                                id, beaching_time, &
                                rho, rho0, &
                                radius, radius0, &
                                h_biofilm, &
                                age, max_age, kill_beached, kill_boundary, &
                                u0, v0, w0, vel_vertical, &
                                traj_len, time_on_beach, is_active, state)
    end do
    close (RESTARTFILE)
    runparts = n_restart_particles

    return
  end subroutine read_restart_file
  !===========================================
  subroutine allocate_n_init_particles(this)
    class(t_initial_position), intent(inout) :: this
    integer :: n_init_p

    n_init_p = this%n_particles

    allocate (this%x(n_init_p), this%y(n_init_p), this%z(n_init_p), &
              this%rho(n_init_p), this%radius(n_init_p), &
              this%beaching_time(n_init_p), this%id(n_init_p), &
              this%is_active(n_init_p))

    this%x = ZERO; this%y = ZERO; this%z = ZERO; 
    this%rho = ZERO; this%radius = ZERO; 
    this%beaching_time = ZERO; this%id = ZERO
    this%is_active = .true.

  end subroutine allocate_n_init_particles
  !===========================================
  subroutine check_initial_coordinates(this, domain)
    class(t_initial_position), intent(inout) :: this
    class(t_domain), intent(in) :: domain
    integer :: ipart, i, j, on_land
    integer :: n_good_particles, n_bad_particles
    integer, allocatable, dimension(:, :) :: seamask
    real(rk) :: lon_l, lon_u, lat_l, lat_u ! domain limits

    n_good_particles = 0
    n_bad_particles = 0

    allocate (seamask(domain%nx, domain%ny))
    seamask = domain%get_seamask()

    lon_l = domain%get_lons(1)
    lon_u = domain%get_lons(domain%nx)
    lat_l = domain%get_lats(1)
    lat_u = domain%get_lats(domain%ny)

    on_land = 0

    START_OMP_DO private(i, j) shared(this, domain) reduction(+:on_land, n_good_particles, n_bad_particles)
    do ipart = 1, this%n_particles
      if ((this%x(ipart) < lon_l) .or. (this%x(ipart) > lon_u) .or. &
          (this%y(ipart) < lat_l) .or. (this%y(ipart) > lat_u)) then
        ERROR, "Particle", ipart, ":", &
          this%x(ipart), this%y(ipart)
        call throw_error("particle_vars :: check_initial_coordinates", "Particle initialised outside of domain")
      end if
      call domain%get_index(this%x(ipart), i=i, dim=1)
      call domain%get_index(this%y(ipart), i=j, dim=2)
      if (seamask(i, j) == DOM_LAND) then
        this%is_active(ipart) = .false.
        on_land = on_land + 1
#ifdef IGNORE_BAD_PARTICLES
        n_bad_particles = n_bad_particles + 1
      else
        n_good_particles = n_good_particles + 1
#endif
      end if
    end do
    END_OMP_DO
    if (on_land > 0) then
      call throw_warning("particle_vars :: check_initial_coordinates", "Particles initialised on land")
      WARNING, on_land, " of ", this%n_particles, " particles on land, time idx = ", this%time_idx
    end if
#ifdef IGNORE_BAD_PARTICLES
    ! Check that n_good_particles + n_bad_particles = n_particles
    if (n_good_particles + n_bad_particles .ne. this%n_particles) then
      call throw_error("particle_vars :: check_initial_coordinates", "n_good_particles + n_bad_particles .ne. n_particles")
    end if
    this%n_good_particles = n_good_particles
#endif

    return
  end subroutine check_initial_coordinates
  !===========================================
  subroutine init_particles_from_coordfile(fieldset)
    !---------------------------------------------
    ! Allocate array for estimated amount of particles
    !---------------------------------------------
    type(t_fieldset), intent(in) :: fieldset
    integer :: ipart

    FMT1, "======== Init particles ========"

    allocate (init_coords(1))

    open (COORDFILE, file=trim(coordfile), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("particle_vars :: init_particles_from_coordfile", "Failed to open "//trim(coordfile), ierr)
    read (COORDFILE, *) init_coords(1)%n_particles
    call init_coords(1)%allocate_n_init_particles
    do ipart = 1, init_coords(1)%n_particles
      read (COORDFILE, *, iostat=ierr) init_coords(1)%x(ipart), init_coords(1)%y(ipart), &
        init_coords(1)%z(ipart), init_coords(1)%id(ipart), &
        init_coords(1)%beaching_time(ipart), init_coords(1)%rho(ipart), init_coords(1)%radius(ipart)
  if (ierr .ne. 0) call throw_error("particle_vars :: init_particles_from_coordfile", "Failed to read from "//trim(coordfile), ierr)
    end do
    close (COORDFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("particle_vars :: init_particles_from_coordfile", "Failed to close "//trim(coordfile), ierr)

    call init_coords(1)%check_initial_coordinates(fieldset%domain)

    ! If inputstep is < 0, it means no additional particles will be released (for restart)
    if (inputstep > 0) then
#ifdef IGNORE_BAD_PARTICLES
      n_particles = init_coords(1)%n_good_particles * (nTimes / inputstep) + init_coords(1)%n_good_particles
#else
      n_particles = init_coords(1)%n_particles * (nTimes / inputstep) + init_coords(1)%n_particles
#endif
    else
      n_particles = 0
    end if

    FMT2, "Finished init particles"

  end subroutine init_particles_from_coordfile
  !===========================================
  subroutine init_particles_from_netcdf(fieldset)
    !---------------------------------------------
    ! Allocate array for estimated amount of particles
    !---------------------------------------------
    type(t_fieldset), intent(in) :: fieldset
    integer :: itime, i_end
    real(rk), allocatable :: nInitParticles(:)

    FMT1, "======== Init particles ========"

    call nc_get_dim_len(trim(coordfile), "time", n_init_times)
    allocate (init_coords(n_init_times), nInitParticles(n_init_times))

    call nc_read_real_1d(trim(coordfile), "n_particles", n_init_times, nInitParticles)

    do itime = 1, n_init_times
      init_coords(itime)%time_idx = itime
      if (itime < n_init_times) then
        init_coords(itime)%next_idx = itime + 1
      else
        ! Could also be periodic (last next_idx = 1)
        ! will stop releasing particles when the init file runs out
        init_coords(itime)%next_idx = 1
        ! init_coords(itime)%next_idx = itime
      end if

      init_coords(itime)%n_particles = int(nInitParticles(itime))
      call init_coords(itime)%allocate_n_init_particles

      init_coords(itime)%release_date = datetime_from_netcdf(trim(coordfile), itime)

      if (nc_var_exists(trim(coordfile), "x")) then
        call nc_read_real_2d(trim(coordfile), "x", 1, int(nInitParticles(itime)), init_coords(itime)%x, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%x = ZERO
      end if

      if (nc_var_exists(trim(coordfile), "y")) then
        call nc_read_real_2d(trim(coordfile), "y", 1, int(nInitParticles(itime)), init_coords(itime)%y, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%y = ZERO
      end if

      if (nc_var_exists(trim(coordfile), "z")) then
        call nc_read_real_2d(trim(coordfile), "z", 1, int(nInitParticles(itime)), init_coords(itime)%z, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%z = ZERO
      end if

      if (nc_var_exists(trim(coordfile), "id")) then
        call nc_read_real_2d(trim(coordfile), "id", 1, int(nInitParticles(itime)), init_coords(itime)%id, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%id = ZERO
      end if

      if (nc_var_exists(trim(coordfile), "beaching_time")) then
        call nc_read_real_2d(trim(coordfile), "beaching_time", 1, int(nInitParticles(itime)), init_coords(itime)%beaching_time, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%beaching_time = 86400.0d0
      end if

      if (nc_var_exists(trim(coordfile), "rho")) then
        call nc_read_real_2d(trim(coordfile), "rho", 1, int(nInitParticles(itime)), init_coords(itime)%rho, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%rho = 30.0d0
      end if

      if (nc_var_exists(trim(coordfile), "radius")) then
        call nc_read_real_2d(trim(coordfile), "radius", 1, int(nInitParticles(itime)), init_coords(itime)%radius, &
                             start=[1, itime], count=[int(nInitParticles(itime)), 1])
      else
        init_coords(itime)%radius = 0.001
      end if

      call init_coords(itime)%check_initial_coordinates(fieldset%domain)

    end do

    do itime = 1, n_init_times
      if (run_start_dt <= init_coords(itime)%release_date) then
        i_release = itime
        FMT2, "First release time:"
        call init_coords(itime)%release_date%print_short_date()
        exit
      end if
    end do
    do itime = 1, n_init_times
      if (run_end_dt <= init_coords(itime)%release_date) then
        i_end = itime
        FMT2, "Final release time:"
        call init_coords(itime)%release_date%print_short_date()
        exit
      end if
    end do

    if (inputstep > 0) then
      if (n_init_times > 1) then
#ifdef IGNORE_BAD_PARTICLES
        ! Count the number of good particles, i.e. particles that are not outside the space or time domain
        do itime = i_release, i_end
          n_particles = n_particles + int(init_coords(itime)%n_good_particles)
        end do
#else
        n_particles = int(sum(nInitParticles))
#endif
      else
#ifdef IGNORE_BAD_PARTICLES
        n_particles = int(init_coords(1)%n_good_particles) * (nTimes / inputstep) + int(nInitParticles(1))
#else
        n_particles = int(nInitParticles(1)) * (nTimes / inputstep) + int(nInitParticles(1))
#endif
      end if
    else
      n_particles = 0
    end if

    FMT2, "Finished init particles"

  end subroutine init_particles_from_netcdf
  !===========================================
  subroutine release_particles(itime, date, fieldset, fieldset_time)
    integer, intent(in) :: itime
    type(t_datetime), intent(in) :: date
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: fieldset_time
    integer :: ipart
    character(len=LEN_CHAR_L) :: info
#ifdef IGNORE_BAD_PARTICLES
    integer :: i_good_part

    i_good_part = 0
#endif

    if (inputstep <= 0) return

    select case (particle_init_method)
    case (TXT_FILE)
      if (mod(itime, inputstep) /= 0) then
        return
      end if
    case (NC_FILE)
      if (date < init_coords(i_release)%release_date) then
        return
      end if
    end select

#ifdef IGNORE_BAD_PARTICLES
    ! FMT2, "Releasing ", init_coords(i_release)%n_good_particles, " new particles at itime = ", itime
    do ipart = 1, init_coords(i_release)%n_particles
      if (init_coords(i_release)%is_active(ipart)) then
        i_good_part = i_good_part + 1
        particles(i_good_part + runparts) = t_particle(lon=init_coords(i_release)%x(ipart), &
                                                       lat=init_coords(i_release)%y(ipart), &
                                                       depth=init_coords(i_release)%z(ipart), &
                                                       id=init_coords(i_release)%id(ipart), &
                                                       beaching_time=init_coords(i_release)%beaching_time(ipart), &
                                                       rho=init_coords(i_release)%rho(ipart), &
                                                       radius=init_coords(i_release)%radius(ipart), &
                                                       max_age=max_age, &
                                                       kill_beached=kill_beached, &
                                                       kill_boundary=kill_boundary, &
                                                       is_active=init_coords(i_release)%is_active(ipart), &
                                                       fieldset=fieldset, &
                                                       time=fieldset_time)
      end if
    end do
    if (i_good_part .ne. init_coords(i_release)%n_good_particles) call throw_error("particle :: release_particles", "i_good_part .ne. init_coords(i_release)%n_good_particles")
    runparts = runparts + init_coords(i_release)%n_good_particles
    write (info, '(a,i9,a,i9,a)') "| "//date%nice_format()//" | Released ", i_good_part, "  particles  | ", runparts, " particles |"
    FMT2, LINE, LINE, LINE
    FMT2, trim(info)
    ! FMT2, runparts, "particles"
    i_release = init_coords(i_release)%next_idx
#else
    FMT2, "Releasing ", init_coords(i_release)%n_particles, " new particles at itime = ", itime
    do ipart = 1, init_coords(i_release)%n_particles
      particles(ipart + runparts) = t_particle(lon=init_coords(i_release)%x(ipart), &
                                               lat=init_coords(i_release)%y(ipart), &
                                               depth=init_coords(i_release)%z(ipart), &
                                               id=init_coords(i_release)%id(ipart), &
                                               beaching_time=init_coords(i_release)%beaching_time(ipart), &
                                               rho=init_coords(i_release)%rho(ipart), &
                                               radius=init_coords(i_release)%radius(ipart), &
                                               max_age=max_age, &
                                               kill_beached=kill_beached, &
                                               kill_boundary=kill_boundary, &
                                               is_active=init_coords(i_release)%is_active(ipart), &
                                               fieldset=fieldset, &
                                               time=fieldset_time)
    end do
    runparts = runparts + init_coords(i_release)%n_particles
    FMT2, runparts, "particles"
    i_release = init_coords(i_release)%next_idx
#endif

    return
  end subroutine release_particles

end module mod_particle_vars
