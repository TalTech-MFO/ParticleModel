#include "cppdefs.h"
module mod_postprocess
  !----------------------------------------------------------------
  ! [module description]
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use mod_variable
  use mod_list
  use mod_domain
  use mod_particle
  use nc_manager
  use netcdf
  use mod_datetime
  use mod_measure
  use mod_accumulator
  use mod_snapshot
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_postprocessor
  !---------------------------------------------
  real(rk), parameter :: pi = 4.*atan(1.)
  !---------------------------------------------
  type t_postprocessor
    private
    real(rk), allocatable :: bin_lon(:)
    real(rk), allocatable :: bin_lat(:)
    real(rk) :: dlon, dlat
    real(rk) :: bin_size_m
    real(rk) :: lonmin, lonmax, latmin, latmax
    integer :: nlon, nlat
    type(t_list) :: measures
    integer :: outputstep = 1
    character(len=LEN_CHAR_L) :: outfile, restart_outfile
    integer :: write_idx = 1
    logical :: should_write_restart = .false.
  contains
    private
    generic, public   :: add_measure => add_basic_measure, add_property_measure
    procedure         :: add_basic_measure, add_property_measure
    procedure, public :: reset_measures
    procedure, public :: after_timestep
    procedure, public :: init_output, init_restart
    procedure, public, nopass :: from_restart
    procedure, public :: save, write_restart
    procedure :: get_index
  end type t_postprocessor
  !---------------------------------------------
  interface t_postprocessor
    module procedure :: ctor_postprocessor
  end interface t_postprocessor
  !===================================================
contains
  !===========================================
  type(t_postprocessor) function ctor_postprocessor(domain, bin_size_m) result(this)
    !---------------------------------------------
    ! Creating the postprocessor object
    ! TODO: Domain resolution
    !---------------------------------------------
    type(t_domain), target, intent(in) :: domain
    real(rk), intent(in) :: bin_size_m
    real(rk) :: meanlat
    integer :: i

    dbghead(postprocess :: ctor_postprocessor)
    FMT1, "======== Init postprocessor ========"

    this%lonmin = domain%get_lons(1); debug(this%lonmin)
    this%lonmax = domain%get_lons(domain%nx); debug(this%lonmax)
    this%latmin = domain%get_lats(1); debug(this%latmin)
    this%latmax = domain%get_lats(domain%ny); debug(this%latmax)
    meanlat = (this%latmin + this%latmax) / 2.0_rk; debug(meanlat)

    this%bin_size_m = bin_size_m
    call domain%xy2lonlat(bin_size_m, bin_size_m, meanlat, this%dlon, this%dlat)

    debug(this%dlat)
    debug(this%dlon)

    this%nlon = int((this%lonmax - this%lonmin) / this%dlon); debug(this%nlon)
    this%nlat = int((this%latmax - this%latmin) / this%dlat); debug(this%nlat)

    ! Creating bins centers
    this%bin_lon = this%lonmin + this%dlon*[(real(i - HALF, rk), i=1, this%nlon + 1)]
    this%bin_lat = this%latmin + this%dlat*[(real(i - HALF, rk), i=1, this%nlat + 1)]

    debug(shape(this%bin_lon))
    debug(minval(this%bin_lon))
    debug(maxval(this%bin_lon))
    debug(shape(this%bin_lat))
    debug(minval(this%bin_lat))
    debug(maxval(this%bin_lat))

    dbgtail(postprocess :: ctor_postprocessor)
  end function ctor_postprocessor
  !===========================================
  type(t_postprocessor) function from_restart(domain, restart_filename) result(this)
    type(t_domain), target, intent(in) :: domain
    character(len=*), intent(in) :: restart_filename
    integer :: nvars
    character(len=LEN_CHAR_L), allocatable :: varnames(:)
    integer :: i
    class(t_variable), pointer :: measure
    character(len=LEN_CHAR_L) :: measure_type, measure_units, measure_var_name
    real(rk), allocatable :: initial_state(:, :)
    real(rk) :: bin_size_m

    dbghead(postprocess :: from_restart)
    FMT1, "======== Init postprocessor (restart) ========"
    FMT1, "Reading postprocessor from restart file: "//trim(restart_filename)

    call nc_get_attr(restart_filename, "bin_size", bin_size_m)
    this = t_postprocessor(domain, bin_size_m)
    allocate (initial_state(this%nlon, this%nlat))
    call nc_get_file_vars(restart_filename, nvars, varnames)
    do i = 1, nvars
      if (varnames(i) == "lon" .or. varnames(i) == "lat" .or. varnames(i) == "time") cycle
      call nc_get_attr(restart_filename, trim(varnames(i)), "measure_type", measure_type)
      call nc_get_attr(restart_filename, trim(varnames(i)), "units", measure_units)
      call nc_get_attr(restart_filename, trim(varnames(i)), "variable_name", measure_var_name)
      if (measure_var_name == "") then
        call this%add_basic_measure(name=trim(varnames(i)), unit=trim(measure_units), measure_type=trim(measure_type)) ! Adding a basic measure
      else
        call this%add_property_measure(name=trim(varnames(i)), unit=trim(measure_units), &
                                       measure_type=trim(measure_type), variable_name=trim(measure_var_name)) ! Adding a property measure
      end if
      call nc_read_real_2d(restart_filename, trim(varnames(i)), this%nlon, this%nlat, initial_state)
      call this%measures%get_item(trim(varnames(i)), measure)
      select type (measure)
      class is (t_measure)
        call measure%set(initial_state)
      end select
    end do

    dbgtail(postprocess :: from_restart)
  end function from_restart
  !===========================================
  subroutine add_basic_measure(this, name, unit, measure_type)
    class(t_postprocessor), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    character(len=*), intent(in) :: measure_type

    select case (measure_type)
    case ("snapshot")
      call this%measures%add_node(name, t_basic_snapshot(name, this%nlon, this%nlat, unit))
    case ("accumulator")
      call this%measures%add_node(name, t_basic_accumulator(name, this%nlon, this%nlat, unit))
    case default
      call throw_error("postprocess :: add_basic_measure", "Unknown measure type: "//measure_type)
    end select

  end subroutine add_basic_measure
  !===========================================
  subroutine add_property_measure(this, name, variable_name, unit, measure_type)
    class(t_postprocessor), intent(inout) :: this
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: variable_name
    character(len=*), intent(in) :: unit
    character(len=*), intent(in) :: measure_type

    select case (measure_type)
    case ("snapshot")
      call this%measures%add_node(name, t_property_snapshot(name, this%nlon, this%nlat, unit, variable_name))
    case ("accumulator")
      call this%measures%add_node(name, t_property_accumulator(name, this%nlon, this%nlat, unit, variable_name))
    case default
      call throw_error("postprocess :: add_property_measure", "Unknown measure type: "//measure_type)
    end select

  end subroutine add_property_measure
  !===========================================
  subroutine reset_measures(this)
    class(t_postprocessor), intent(inout) :: this
    integer :: imeasure
    class(t_variable), pointer :: measure

    dbghead(postprocess :: reset_measures)

    do imeasure = 1, this%measures%size()
      call this%measures%get_item(imeasure, measure)
      select type (measure)
      class is (t_measure)
        call measure%reset()
      end select
    end do

    dbgtail(postprocess :: reset_measures)
  end subroutine reset_measures
  !===========================================
  subroutine init_output(this, outfile, outputstep)
    class(t_postprocessor), intent(inout) :: this
    character(len=LEN_CHAR_L), intent(in) :: outfile
    integer, intent(in) :: outputstep
    integer :: nc_t_dimid, nc_lon_dimid, nc_lat_dimid
    integer :: imeasure
    class(t_variable), pointer :: measure

    dbghead(postprocess :: init_output)
    FMT1, "======== Init postprocessor output ========"

    this%outfile = trim(outfile)
    this%outputstep = outputstep

    FMT2, "Postprocessing output file: ", trim(this%outfile)
    FMT2, "Output step: ", this%outputstep

    call nc_initialise(outfile)

    call nc_add_dimension(outfile, "time", nc_t_dimid)
    call nc_add_dimension(outfile, "lon", nc_lon_dimid, size(this%bin_lon) - 1)
    call nc_add_dimension(outfile, "lat", nc_lat_dimid, size(this%bin_lat) - 1)

    call nc_add_variable(outfile, "time", "float", 1, [nc_t_dimid])
    call nc_add_attr(outfile, "time", "units", "seconds since 1900-01-01 00:00:00")

    call nc_add_variable(outfile, "lon", "float", 1, [nc_lon_dimid], FILLVALUE_BIG)
    call nc_add_attr(outfile, "lon", "units", "degrees_east")
    call nc_add_attr(outfile, "lon", "minval", this%lonmin)
    call nc_add_attr(outfile, "lon", "maxval", this%lonmax)
    call nc_write(outfile, this%bin_lon(:size(this%bin_lon) - 1), "lon", size(this%bin_lon) - 1)

    call nc_add_variable(outfile, "lat", "float", 1, [nc_lat_dimid], FILLVALUE_BIG)
    call nc_add_attr(outfile, "lat", "units", "degrees_north")
    call nc_add_attr(outfile, "lat", "minval", this%latmin)
    call nc_add_attr(outfile, "lat", "maxval", this%latmax)
    call nc_write(outfile, this%bin_lat(:size(this%bin_lat) - 1), "lat", size(this%bin_lat) - 1)

    do imeasure = 1, this%measures%size()
      call this%measures%get_item(imeasure, measure)
      select type (measure)
      class is (t_measure)
        call nc_add_variable(outfile, measure%get_name(), "float", 3, [nc_lon_dimid, nc_lat_dimid, nc_t_dimid], FILLVALUE_BIG)
        call nc_add_attr(outfile, measure%get_name(), "description", measure%get_measure_description())
        call nc_add_attr(outfile, measure%get_name(), "variable_name", measure%get_variable_name())
        call nc_add_attr(outfile, measure%get_name(), "measure_type", measure%get_measure_type())
        call nc_add_attr(outfile, measure%get_name(), "units", measure%get_units())
      end select
    end do

    call nc_add_attr(outfile, "bin_size", this%bin_size_m)

    dbgtail(postprocess :: init_output)
  end subroutine init_output
  !===========================================
  subroutine init_restart(this, outfile)
    class(t_postprocessor), intent(inout) :: this
    character(len=LEN_CHAR_L), intent(in) :: outfile
    integer :: nc_t_dimid, nc_lon_dimid, nc_lat_dimid
    integer :: imeasure
    class(t_variable), pointer :: measure

    dbghead(postprocess :: init_restart)

    FMT1, "======== Init postprocessor restart output ========"

    this%restart_outfile = trim(outfile)
    this%should_write_restart = .true.

    FMT2, "Postprocessing restart file: ", trim(this%restart_outfile)

    call nc_initialise(outfile)

    call nc_add_dimension(outfile, "time", nc_t_dimid, 1)
    call nc_add_dimension(outfile, "lon", nc_lon_dimid, size(this%bin_lon) - 1)
    call nc_add_dimension(outfile, "lat", nc_lat_dimid, size(this%bin_lat) - 1)

    call nc_add_variable(outfile, "time", "float", 1, [nc_t_dimid])
    call nc_add_attr(outfile, "time", "units", "seconds since 1900-01-01 00:00:00")

    call nc_add_variable(outfile, "lon", "float", 1, [nc_lon_dimid], FILLVALUE_BIG)
    call nc_add_attr(outfile, "lon", "units", "degrees_east")
    call nc_add_attr(outfile, "lon", "minval", this%lonmin)
    call nc_add_attr(outfile, "lon", "maxval", this%lonmax)
    call nc_write(outfile, this%bin_lon(:size(this%bin_lon) - 1), "lon", size(this%bin_lon) - 1)

    call nc_add_variable(outfile, "lat", "float", 1, [nc_lat_dimid], FILLVALUE_BIG)
    call nc_add_attr(outfile, "lat", "units", "degrees_north")
    call nc_add_attr(outfile, "lat", "minval", this%latmin)
    call nc_add_attr(outfile, "lat", "maxval", this%latmax)
    call nc_write(outfile, this%bin_lat(:size(this%bin_lat) - 1), "lat", size(this%bin_lat) - 1)

    do imeasure = 1, this%measures%size()
      call this%measures%get_item(imeasure, measure)
      select type (measure)
      class is (t_measure)
        call nc_add_variable(outfile, measure%get_name(), "float", 3, [nc_lon_dimid, nc_lat_dimid, nc_t_dimid], FILLVALUE_BIG)
        call nc_add_attr(outfile, measure%get_name(), "description", measure%get_measure_description())
        call nc_add_attr(outfile, measure%get_name(), "variable_name", measure%get_variable_name())
        call nc_add_attr(outfile, measure%get_name(), "measure_type", measure%get_measure_type())
        call nc_add_attr(outfile, measure%get_name(), "units", measure%get_units())
      end select
    end do

    call nc_add_attr(outfile, "bin_size", this%bin_size_m)

    dbgtail(postprocess :: init_restart)
  end subroutine init_restart
  !===========================================
  subroutine after_timestep(this, p)
    class(t_postprocessor), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer :: i, j, imeasure
    ! real(rk) :: p_age, p_distance
    class(t_variable), pointer :: measure

    dbghead(postprocessor :: after_timestep)

    call this%get_index(p%lon0, p%lat0, i, j)

    debug(i); debug(j)

    do imeasure = 1, this%measures%size()
      call this%measures%get_item(imeasure, measure)
      select type (measure)
      class is (t_measure)
        call measure%run(p, i, j)
      end select
    end do

    dbgtail(postprocessor :: after_timestep)
  end subroutine after_timestep
  !===========================================
  subroutine get_index(this, lon, lat, i, j)
    class(t_postprocessor), intent(in) :: this
    real(rk), intent(in) :: lon, lat
    integer, intent(out) :: i, j

    dbghead(postprocess :: get_index)

    debug(lon); debug(lat)

    i = floor((lon - this%lonmin) / this%dlon) + 1; debug(i)
    j = floor((lat - this%latmin) / this%dlat) + 1; debug(j)

    if (i < 1) then
      DBG, "i < 1: i = ", i
      i = 1
    elseif (i > size(this%bin_lon) - 1) then
      DBG, "i > bin_lon: i = ", i, size(this%bin_lon) - 1
      i = size(this%bin_lon) - 1
    end if

    if (j < 1) then
      DBG, "j < 1: j = ", j
      j = 1
    elseif (j > size(this%bin_lat) - 1) then
      DBG, "j > bin_lat: i = ", i, size(this%bin_lat) - 1
      j = size(this%bin_lat) - 1
    end if

    dbgtail(postprocess :: get_index)
  end subroutine get_index
  !===========================================
  subroutine save(this, itime, date)
    class(t_postprocessor), intent(inout) :: this
    integer, intent(in) :: itime
    type(t_datetime), intent(in) :: date
    real(rk) :: dateval(1)
    integer :: ncid, varid, imeasure
    class(t_variable), pointer :: measure
    ! character(len=LEN_CHAR_L) :: info

    if (mod(itime, this%outputstep) /= 0) return

    ! write (info, '(a,i11,a,i9,a)') "| "//date%nice_format()//" | Writing post output ", itime, " | ", this%measures%size(), " measures |"
!       FMT2, LINE, LINE, LINE
    !       FMT2, trim(info)

    call nc_check(trim(this%outfile), nf90_open(trim(this%outfile), nf90_write, ncid), "postprocessor :: save :: open")
    call nc_check(trim(this%outfile), nf90_inq_varid(ncid, "time", varid), "postprocessor :: save :: inq_varid")
    dateval = date%date2num()
    call nc_check(trim(this%outfile), &
                  nf90_put_var(ncid, varid, dateval, start=[this%write_idx], count=[1]), "postprocessor :: save :: put_var")

    do imeasure = 1, this%measures%size()
      call this%measures%get_item(imeasure, measure)
      select type (measure)
      class is (t_measure)
        call nc_check(trim(this%outfile), nf90_inq_varid(ncid, measure%get_name(), varid), "postprocessor :: save :: inq_varid")
        call nc_check(trim(this%outfile), &
                      nf90_put_var(ncid, varid, measure%get(), &
                                   start=[1, 1, this%write_idx], &
                                   count=[size(this%bin_lon) - 1, size(this%bin_lat) - 1, 1]), &
                      "postprocessor :: save :: put_var")
      end select
    end do
    call nc_check(trim(this%outfile), nf90_close(ncid), "postprocessor :: save :: close")
    this%write_idx = this%write_idx + 1

  end subroutine save
  !===========================================
  subroutine write_restart(this, date)
    class(t_postprocessor), intent(inout) :: this
    type(t_datetime), intent(in) :: date
    real(rk) :: dateval(1)
    integer :: ncid, varid, imeasure
    class(t_variable), pointer :: measure

    if (.not. this%should_write_restart) return

 call nc_check(trim(this%restart_outfile), nf90_open(trim(this%restart_outfile), nf90_write, ncid), "postprocessor :: write_restart :: open")
    call nc_check(trim(this%restart_outfile), nf90_inq_varid(ncid, "time", varid), "postprocessor :: write_restart :: inq_varid")
    dateval = date%date2num()
    call nc_check(trim(this%restart_outfile), &
                  nf90_put_var(ncid, varid, dateval, start=[1], count=[1]), "postprocessor :: write_restart :: put_var")

    do imeasure = 1, this%measures%size()
      call this%measures%get_item(imeasure, measure)
      select type (measure)
      class is (t_measure)
    call nc_check(trim(this%restart_outfile), nf90_inq_varid(ncid, measure%get_name(), varid), "postprocessor :: write_restart :: inq_varid")
        call nc_check(trim(this%restart_outfile), &
                      nf90_put_var(ncid, varid, measure%get(), &
                                   start=[1, 1, 1], &
                                   count=[size(this%bin_lon) - 1, size(this%bin_lat) - 1, 1]), &
                      "postprocessor :: write_restart :: put_var")
      end select
    end do
    call nc_check(trim(this%restart_outfile), nf90_close(ncid), "postprocessor :: write_restart :: close")

  end subroutine write_restart

end module mod_postprocess
