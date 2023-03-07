#include "cppdefs.h"
#include "field.h"
module mod_domain
  !----------------------------------------------------------------
  ! Initialise domain and seamask
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use nc_manager, only: nc_read_real_1d, nc_read_real_2d, nc_get_var_fillvalue
  use mod_field
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_domain
  !---------------------------------------------
  real(rk), parameter :: pi = 4.*atan(1.)
  real(rk), parameter :: radius_earth = 6371.0e3_rk
  !---------------------------------------------
  type t_domain
    private
    integer, public                 :: nx, ny
    type(t_field_static_1d)         :: lons, lats
    type(t_field_static_2d)         :: depdata
    integer, allocatable            :: seamask(:, :) ! field type does not take integer
    real(rk)                        :: lboundx, lboundy, uboundx, uboundy
    type(t_field_static_2d), public :: dlon, dlat, dx, dy
  contains
    private
    procedure, public :: lonlat2xy, xy2lonlat
    procedure, public :: get_index
    generic, public   :: get_bathymetry => get_bathymetry_whole, get_bathymetry_idx, get_bathymetry_idx_interp
    procedure         :: get_bathymetry_whole, get_bathymetry_idx, get_bathymetry_idx_interp
    generic, public   :: get_seamask => get_seamask_whole, get_seamask_idx
    procedure         :: get_seamask_whole, get_seamask_idx
    generic, public   :: get_lons => get_lons_whole, get_lons_idx, get_lons_interp
    procedure         :: get_lons_whole, get_lons_idx, get_lons_interp
    generic, public   :: get_lats => get_lats_whole, get_lats_idx, get_lats_interp
    procedure         :: get_lats_whole, get_lats_idx, get_lats_interp

  end type t_domain
  !---------------------------------------------
  interface t_domain
    module procedure :: ctor_domain
  end interface t_domain
  !===================================================
contains
  !===========================================
  type(t_domain) function ctor_domain(nx, ny, topofile, lon, lat, bathy) result(d)
    !---------------------------------------------
    ! Initialize the global longitude/latitude
    ! and seamask
    ! TODO: dx, dy, dlon, dlat should be 2D arrays
    ! TODO: Should get dimensions from netCDF file
    ! TODO: Should check topo file for seamask (with the same values for land, sea etc.)
    !---------------------------------------------
#ifdef DEBUG
    use nc_manager
    integer :: nc_x_dimid, nc_y_dimid
#endif
    integer, intent(in) :: nx, ny
    character(len=*), intent(in) :: topofile
    character(len=*), intent(in) :: lon, lat, bathy ! Variable names
    integer :: ii, jj
    real(rk) :: fill_value
    real(rk), allocatable :: lons(:), lats(:)
    real(rk), dimension(:, :), allocatable :: dlon, dlat, dx, dy, depdata
    real(rk) :: dx_min, dx_max, dy_min, dy_max
    real(rk) :: dlon_min, dlon_max, dlat_min, dlat_max

    dbghead(domain :: domain)

    d%nx = nx
    d%ny = ny

    FMT2, "Reading coordinates"
    allocate (lons(nx), lats(ny))

    call nc_read_real_1d(trim(topofile), trim(lon), nx, lons)
    call d%lons%init(n=nx, name=lon, value=lons, units="degrees_east")

    call nc_read_real_1d(trim(topofile), trim(lat), ny, lats)
    call d%lats%init(n=ny, name=lat, value=lats, units="degrees_north")

    d%lboundx = lons(1)
    d%lboundy = lats(1)
    d%uboundx = lons(nx)
    d%uboundy = lats(ny)

    allocate (dlon(nx, ny), dlat(nx, ny), dx(nx, ny), dy(nx, ny))

    do jj = 1, ny - 1
      do ii = 1, nx - 1
        dlon(ii, jj) = lons(ii + 1) - lons(ii)
        dlat(ii, jj) = lats(jj + 1) - lats(jj)
        dx(ii, jj) = dlon(ii, jj) * ((radius_earth * pi / 180._rk) * cos(lats(jj) * pi / 180._rk))
        dy(ii, jj) = dlat(ii, jj) * (radius_earth * pi / 180._rk)
      end do
    end do
    dlon(:, ny) = dlon(:, ny - 1)
    dlat(:, ny) = dlat(:, ny - 1)
    dx(:, ny) = dx(:, ny - 1)
    dy(:, ny) = dy(:, ny - 1)
    dlon(nx, :) = dlon(nx - 1, :)
    dlat(nx, :) = dlat(nx - 1, :)
    dx(nx, :) = dx(nx - 1, :)
    dy(nx, :) = dy(nx - 1, :)

    call d%dlon%init(n1=nx, n2=ny, name="dlon", value=dlon, units="degrees_east")
    call d%dlat%init(n1=nx, n2=ny, name="dlat", value=dlat, units="degrees_north")
    call d%dx%init(n1=nx, n2=ny, name="dx", value=dx, units="meters")
    call d%dy%init(n1=nx, n2=ny, name="dy", value=dy, units="meters")

    FMT2, LINE
    FMT2, "Coordinates:"
    FMT3, var2val(d%lboundy), "[deg N], ", var2val(d%uboundy), "[deg N]"
    FMT3, var2val(d%lboundx), "[deg E], ", var2val(d%uboundx), "[deg E]"
    FMT2, "Cell size:"
    dy_min = minval(dy)
    dy_max = maxval(dy)
    dlat_min = minval(dlat)
    dlat_max = maxval(dlat)
    dx_min = minval(dx)
    dx_max = maxval(dx)
    dlon_min = minval(dlon)
    dlon_max = maxval(dlon)
    FMT3, "min: ", dlat_min, ", max: ", dlat_max, "[deg N]; min: ", dy_min, ", max: ", dy_max, "[m]"
    FMT3, "min: ", dlon_min, ", max: ", dlon_max, "[deg E]; min: ", dx_min, ", max: ", dx_max, "[m]"

    allocate (depdata(nx, ny), d%seamask(nx, ny))

    FMT2, "Reading bathymetry"
    call nc_read_real_2d(trim(topofile), trim(bathy), nx, ny, depdata)
    if (any(isnan(depdata))) then
      FMT3, "Bathymetry contains NaNs. Not trying to find fill value."
    else
      if (.not. nc_get_var_fillvalue(trim(topofile), trim(bathy), fill_value)) then
        ! Assuming that the fill value is a large negative number
        fill_value = minval(depdata)
        call throw_warning("domain :: domain", "Could not find fill value for bathymetry. Assuming that the fill value is a large negative number.")
      end if
      FMT3, "Fill value for bathymetry is assumed to be ", fill_value
      where (depdata == fill_value) depdata = ieee_value(ZERO, ieee_quiet_nan)
    end if
    call d%depdata%init(n1=nx, n2=ny, name=bathy, value=depdata, units="meters")

    !---------------------------------------------
    ! TODO: Seamask could have another value (4) to represent boundaries.
    !       Boundary should have a thickness!
    FMT2, "Making seamask"

    do ii = 2, nx - 1
      do jj = 2, ny - 1
        if (.not. isnan(depdata(ii, jj))) then
          if ((isnan(depdata(ii + 1, jj))) .or. (isnan(depdata(ii - 1, jj))) .or. &
              (isnan(depdata(ii, jj + 1))) .or. (isnan(depdata(ii, jj - 1))) .or. &
              (isnan(depdata(ii + 1, jj + 1))) .or. (isnan(depdata(ii + 1, jj - 1))) .or. &
              (isnan(depdata(ii - 1, jj - 1))) .or. (isnan(depdata(ii - 1, jj + 1)))) then
            d%seamask(ii, jj) = DOM_BEACH
          else
            d%seamask(ii, jj) = DOM_SEA
          end if
        else
          d%seamask(ii, jj) = DOM_LAND
        end if
      end do
    end do
    do ii = 1, nx
      if (.not. isnan(depdata(ii, 1))) then
        d%seamask(ii, 1) = DOM_BOUNDARY
      else
        d%seamask(ii, 1) = DOM_LAND
      end if
      if (.not. isnan(depdata(ii, ny))) then
        d%seamask(ii, ny) = DOM_BOUNDARY
      else
        d%seamask(ii, ny) = DOM_LAND
      end if
    end do
    do jj = 1, ny
      if (.not. isnan(depdata(1, jj))) then
        d%seamask(1, jj) = DOM_BOUNDARY
      else
        d%seamask(1, jj) = DOM_LAND
      end if
      if (.not. isnan(depdata(nx, jj))) then
        d%seamask(nx, jj) = DOM_BOUNDARY
      else
        d%seamask(nx, jj) = DOM_LAND
      end if
    end do

    FMT2, sum(d%seamask, mask=d%seamask == DOM_LAND), " land points"
    FMT2, sum(d%seamask, mask=d%seamask == DOM_SEA) / DOM_SEA, " sea points"
    FMT2, sum(d%seamask, mask=d%seamask == DOM_BEACH) / DOM_BEACH, " beach points"
    FMT2, sum(d%seamask, mask=d%seamask == DOM_BOUNDARY) / DOM_BOUNDARY, " boundary points"
    FMT2, nx * ny, " total points"

#ifdef DEBUG
#define FNAME "seamask.nc"
    call nc_initialise(FNAME)
    call nc_add_dimension(FNAME, "lon", nc_x_dimid, nx)
    call nc_add_dimension(FNAME, "lat", nc_y_dimid, ny)
    call nc_add_variable(FNAME, "seamask", "int", 2, [nc_x_dimid, nc_y_dimid])
    call nc_add_variable(FNAME, "lon", "float", 1, [nc_x_dimid])
    call nc_add_variable(FNAME, "lat", "float", 1, [nc_y_dimid])
    call nc_write(FNAME, d%seamask, "seamask", nx, ny)
    call nc_write(FNAME, d%lons%get(), "lon", nx)
    call nc_write(FNAME, d%lats%get(), "lat", ny)
#undef FNAME
#endif

    deallocate (lons, lats)
    deallocate (dlon, dlat, dx, dy)
    deallocate (depdata)

    FMT2, "Finished init_domain"

    dbgtail(domain :: domain)
    return
  end function ctor_domain
  !===========================================
  function get_lons_whole(this) result(res)
    class(t_domain), intent(in) :: this
    real(rk), dimension(this%nx):: res

    res = this%lons%get()

    return
  end function get_lons_whole
  !===========================================
  real(rk) function get_lons_idx(this, idx) result(res)
    class(t_domain), intent(in) :: this
    integer AIM_INTENT :: idx
    logical :: out_of_bounds

    dbghead(domain :: get_lons_idx)

    debug(idx)

    res = this%lons%get(idx, out_of_bounds)
    debug(res); debug(out_of_bounds)
    if (out_of_bounds) then
      if (idx < 1) then
        res = this%lboundx
#ifdef SNAP_TO_BOUNDS
        idx = 1
#else
        call throw_error("domain :: get_lons_idx", "Index out of bounds! (Less than 1)")
#endif
      else
        res = this%uboundx
#ifdef SNAP_TO_BOUNDS
        ERROR, "Snapping to bound"
        ERROR, "res = ", res
        ERROR, "idx = ", idx
        idx = this%nx
#else
        call throw_error("domain :: get_lons_idx", "Index out of bounds! (Greater than nx)")
#endif

      end if
    end if

    return
  end function get_lons_idx
  !===========================================
  real(rk) function get_lons_interp(this, idx) result(res)
    class(t_domain), intent(in) :: this
    real(rk) AIM_INTENT :: idx
    logical :: out_of_bounds

    res = this%lons%get(idx, out_of_bounds)
    if (out_of_bounds) then
      if (idx < 1) then
        res = this%lboundx
#ifdef SNAP_TO_BOUNDS
        idx = ONE
#else
        call throw_error("domain :: get_lons_interp", "Index out of bounds! (Less than 1)")
#endif
      else
        res = this%uboundx
#ifdef SNAP_TO_BOUNDS
        idx = real(this%nx, rk)
#else
        call throw_error("domain :: get_lons_interp", "Index out of bounds! (Greater than nx)")
#endif
      end if
    end if

    return
  end function get_lons_interp
  !===========================================
  function get_lats_whole(this) result(res)
    class(t_domain), intent(in)  :: this
    real(rk), dimension(this%ny) :: res

    res = this%lats%get()

    return
  end function get_lats_whole
  !===========================================
  real(rk) function get_lats_idx(this, idx) result(res)
    class(t_domain), intent(in) :: this
    integer AIM_INTENT :: idx
    logical :: out_of_bounds

    res = this%lats%get(idx, out_of_bounds)
    if (out_of_bounds) then
      if (idx < 1) then
        res = this%lboundy
#ifdef SNAP_TO_BOUNDS
        idx = 1
#else
        call throw_error("domain :: get_lats_idx", "Index out of bounds! (Less than 1)")
#endif
      else
        res = this%uboundy
#ifdef SNAP_TO_BOUNDS
        idx = this%ny
#else
        call throw_error("domain :: get_lats_idx", "Index out of bounds! (Greater than ny)")
#endif
      end if
    end if

    return
  end function get_lats_idx
  !===========================================
  real(rk) function get_lats_interp(this, idx) result(res)
    class(t_domain), intent(in) :: this
    real(rk) AIM_INTENT :: idx
    real(rk) :: i0, i1
    logical :: out_of_bounds

    res = this%lats%get(idx, out_of_bounds)
    if (out_of_bounds) then
      if (idx < 1) then
        res = this%lboundy
#ifdef SNAP_TO_BOUNDS
        idx = ONE
#else
        call throw_error("domain :: get_lats_interp", "Index out of bounds! (Less than 1)")
#endif
      else
        res = this%uboundy
#ifdef SNAP_TO_BOUNDS
        idx = real(this%ny, rk)
#else
        call throw_error("domain :: get_lats_interp", "Index out of bounds! (Greater than ny)")
#endif
      end if
    end if

    return
  end function get_lats_interp
  !===========================================
  function get_bathymetry_whole(this) result(res)
    class(t_domain), intent(in) :: this
    real(rk), dimension(this%nx, this%ny) :: res

    res = this%depdata%get()

    return
  end function get_bathymetry_whole
  !===========================================
  function get_bathymetry_idx(this, i, j) result(res)
    class(t_domain), intent(in) :: this
    integer AIM_INTENT :: i, j
    real(rk) :: res
    logical :: out_of_bounds

    res = this%depdata%get(i, j, out_of_bounds)
    if (out_of_bounds) then
#ifdef SNAP_TO_BOUNDS
      if (i < 1) i = 1
      if (i > this%nx) i = this%nx
      if (j < 1) j = 1
      if (j > this%ny) j = this%ny
      res = this%depdata%get(i, j)
#else
      call throw_error("domain :: get_bathymetry_idx", "Index out of bounds!")
#endif
    end if

    return
  end function get_bathymetry_idx
  !===========================================
  function get_bathymetry_idx_interp(this, x, y) result(res)
    class(t_domain), intent(in) :: this
    real(rk), intent(in)        :: x, y ! Indices, not coordinates!
    real(rk)                    :: x1, x2, &
                                   y1, y2, &
                                   c11, c12, c21, c22
    integer                     :: i, j
    real(rk)                    :: res
    logical :: out_of_bounds

    res = this%depdata%get(x, y, out_of_bounds)
    if (out_of_bounds) call throw_error("domain :: get_bathymetry_idx_interp", "Index out of bounds!")

    return
  end function get_bathymetry_idx_interp
  !===========================================
  function get_seamask_whole(this) result(res)
    class(t_domain), intent(in) :: this
    integer, dimension(this%nx, this%ny) :: res

    res = this%seamask

  end function get_seamask_whole
  !===========================================
  function get_seamask_idx(this, i, j) result(res)
    class(t_domain), intent(in) :: this
    integer AIM_INTENT          :: i, j
    integer                     :: res

#ifdef SNAP_TO_BOUNDS
    if (i < 1) i = 1
    if (i > this%nx) i = this%nx
    if (j < 1) j = 1
    if (j > this%ny) j = this%ny
#endif

    if (i < 1) call throw_error("domain :: get_seamask_idx", "i is less than 1")
    if (i > this%nx) call throw_error("domain :: get_seamask_idx", "i is greater than nx")
    if (j < 1) call throw_error("domain :: get_seamask_idx", "j is less than 1")
    if (j > this%ny) call throw_error("domain :: get_seamask_idx", "j is greater than ny")

    res = this%seamask(i, j)

    return
  end function get_seamask_idx
  !===========================================
  subroutine lonlat2xy(this, lon, lat, x, y)
    !---------------------------------------------
    ! Convert longitude and latitude to x and y coordinates 
    ! using lboundx and lboundy as the origin
    !---------------------------------------------
    class(t_domain), intent(in) :: this
    real(rk), intent(in)  :: lon, lat
    real(rk), intent(out) :: x, y

    x = (lon - this%lboundx) * ((radius_earth * pi / 180._rk) * cos(lat * pi / 180._rk))
    y = (lat - this%lboundy) * (radius_earth * pi / 180._rk)

    return
  end subroutine lonlat2xy
  !===========================================
  subroutine xy2lonlat(this, x, y, lon, lat)
    !---------------------------------------------
    ! Convert x and y coordinates to longitude and latitude 
    ! using lboundx and lboundy as the origin
    !---------------------------------------------
    class(t_domain), intent(in) :: this
    real(rk), intent(in) :: x, y
    real(rk), intent(out) :: lon, lat

    lat = y / (radius_earth * pi / 180._rk) + this%lboundy
    lon = x / ((radius_earth * pi / 180._rk) * cos(lat * pi / 180._rk)) + this%lboundx

    return
  end subroutine xy2lonlat
  !===========================================
  subroutine get_index(this, loc, dim, i, ir)
    class(t_domain), intent(in) :: this
    real(rk), intent(in) :: loc ! Location in degrees (lon or lat)
    integer, intent(in):: dim
    integer, intent(out), optional :: i
    real(rk), intent(out), optional  :: ir
    integer :: it
    real(rk) :: irt

    dbghead(domain :: get_index)

    select case (dim)
    case (1)
      DBG, "dim = 1"
      ! Longitude
      ! Check bounds
      if (loc < this%lboundx) then
        DBG, "loc < this%lboundx"
#ifdef SNAP_TO_BOUNDS
        if (present(i)) i = 1
        if (present(ir)) ir = ONE
        dbgtail(domain :: get_index)
        return
#else
        call throw_error("domain :: get_index", "lon is less than lboundx")
#endif
      end if
      if (loc > this%uboundx) then
        DBG, "loc > this%uboundx"
#ifdef SNAP_TO_BOUNDS
        if (present(i)) i = this%nx
        if (present(ir)) ir = real(this%nx, rk)
        dbgtail(domain :: get_index)
        return
#else
        call throw_error("domain :: get_index", "lon is greater than uboundx")
#endif
      end if
      ! Get indices
      it = minloc(abs(this%lons%get() - loc), dim=1)
      debug(it)
      if (this%lons%get(it) > loc) then
        DBG, "this%lons(it) > loc"
        it = it - 1
        debug(it)
      end if
      ! Get real indices
      irt = it + (loc - this%lons%get(it)) / (this%lons%get(it + 1) - this%lons%get(it))
      debug(irt)
      if (present(i)) i = it
      if (present(ir)) ir = irt
    case (2)
      DBG, "dim = 2"
      ! Latitude
      ! Check bounds
      if (loc < this%lboundy) then
        DBG, "loc < this%lboundy"
#ifdef SNAP_TO_BOUNDS
        if (present(i)) i = 1
        if (present(ir)) ir = ONE
        dbgtail(domain :: get_index)
        return
#else
        call throw_error("domain :: get_index", "lat is less than lboundy")
#endif
      end if
      if (loc > this%uboundy) then
        DBG, "loc > this%uboundy"
#ifdef SNAP_TO_BOUNDS
        if (present(i)) i = this%ny
        if (present(ir)) ir = real(this%ny, rk)
        dbgtail(domain :: get_index)
        return
#else
        call throw_error("domain :: get_index", "lat is greater than uboundy")
#endif
      end if
      ! Get indices
      it = minloc(abs(this%lats%get() - loc), dim=1)
      debug(it)
      if (this%lats%get(it) > loc) then
        DBG, "this%lats(it) > loc"
        it = it - 1
        debug(it)
      end if
      ! Get real indices
      irt = it + (loc - this%lats%get(it)) / (this%lats%get(it + 1) - this%lats%get(it))
      debug(irt)
      if (present(i)) i = it
      if (present(ir)) ir = irt
    case default
      call throw_error("domain :: get_index", "dim must be 1 or 2")
    end select

    dbgtail(domain :: get_index)
    return
  end subroutine get_index

end module mod_domain
