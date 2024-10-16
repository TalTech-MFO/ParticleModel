#include "cppdefs.h"
#include "field.h"
#include "file.h"
module mod_fieldset
  use mod_precdefs
  use mod_errors
  use mod_datetime
  use nc_manager
  use mod_variable, only: t_variable
  use mod_field
  use mod_domain, only: t_domain
  use mod_list, only: t_list
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_fieldset
  !---------------------------------------------
  character(len=*), parameter     :: dirinfile = 'dirlist.dat' ! List the dirs into this file
  !---------------------------------------------
  type t_fieldset
    private
    !---------------------------------------------
    type(t_list)                    :: fields
    type(t_domain), pointer, public :: domain
    integer, public                 :: num_fields = 0
    !---------------------------------------------
    ! Time variables
    real(rk)                        :: nc_timestep
    type(t_datetime)                :: date_t1, date_t2
    !---------------------------------------------
    ! Dimension variables
    integer, public                        :: nx, ny, nz
    integer, allocatable                   :: dims(:)
    character(len=LEN_CHAR_S), allocatable :: dimnames(:)
    logical                                :: has_vertical
    !---------------------------------------------
    ! Z axis variables
    integer         :: zax_style = DEPTH_VALUES
    integer         :: zax_dir = 1 ! Default: positive up
    integer         :: zax_idx = 0
    integer, public :: zax_bot_idx = 1, zax_top_idx = 1 ! Calling this bottom and top (not surface) to be more general (e.g. for meteo too)
    !---------------------------------------------
    ! Velocity component variables
    integer :: u_idx = 0
    integer :: v_idx = 0
    integer, allocatable :: u_mask(:, :)
    integer, allocatable :: v_mask(:, :)
    !---------------------------------------------
    ! netCDF read variables
    logical                   :: read_first = .true.
    logical                   :: has_more = .true.
    integer                   :: read_idx
    integer                   :: read_idx_increment = 1
    character(len=LEN_CHAR_L) :: current_path   ! Path to current data
    integer                   :: current_ntimes ! Current timestep index
    integer                   :: dirlist_idx
    type(t_datetime)          :: next_read_dt
    !---------------------------------------------
    ! Paths/directories
    character(len=LEN_CHAR_L)              :: PATH, PMAPFILE           ! Path to field data and processor map
    character(len=LEN_CHAR_L)              :: file_prefix, file_suffix ! What comes before and after the proc. number?
    character(len=LEN_CHAR_L), allocatable :: filelist(:)              ! List of files
    integer, allocatable                   :: dirlist(:)               ! List of directories (one or the other will be allocated)
    integer                                :: nentries                 ! Number of directories or files
    !---------------------------------------------
    ! Subdomain variables
    logical              :: has_subdomains = .false.
    integer              :: nproc       ! Number of subdomains
    integer              :: nxp, nyp    ! Size of subdomain
    integer, allocatable :: pmap(:, :)  ! Offsets of subdomains
    integer, allocatable :: pmask(:, :) ! Proc. mask

  contains
    private
    procedure         :: init_dimensions
    procedure         :: init_time
    procedure         :: init_subdomain
    procedure         :: init_dirlist
    procedure         :: init_u_mask
    procedure         :: init_v_mask
    procedure, public :: add_field
    procedure, public :: list_fields
    procedure, public :: get_field_index
    procedure, public :: has_field
    procedure, public :: set_u_component
    procedure, public :: set_v_component
    procedure, public :: set_zax
    procedure         :: get_zax_regular, get_zax_adaptive
    generic, public   :: get_zax => get_zax_regular, get_zax_adaptive
    procedure, public :: top_is_nan
    procedure, public :: bottom_is_nan
    procedure, public :: zax_direction
    procedure, public :: set_start_time
    procedure, public :: set_simulation_timestep
    procedure, public :: get_directory, get_file
    procedure, public :: get_time
    procedure, public :: find_directory, find_file
    procedure, public :: get_pmap
    procedure         :: get_indices_vertical
    procedure, public :: search_indices
    procedure, public :: sealevel
    generic, public   :: set => set_value_key, set_value_idx
    procedure         :: set_value_key, set_value_idx
    generic, public   :: get => get_value_key_real_idx, get_value_key_int_idx, get_value_idx_real_idx, get_value_idx_int_idx
    procedure         :: get_value_key_real_idx, get_value_key_int_idx, get_value_idx_real_idx, get_value_idx_int_idx
    procedure         :: read_field
    procedure         :: read_field_subdomains
    procedure, public :: update
    procedure, public :: read_first_timesteps

  end type t_fieldset
  !---------------------------------------------
  interface t_fieldset
    module procedure :: ctor_fieldset
  end interface t_fieldset
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  ! CONSTRUCTOR
  !===========================================
  type(t_fieldset) function ctor_fieldset(path, domain, dimnames, file_prefix, file_suffix, start, dt, pmap) result(f)
#undef PROC0
#define PROC0 "/"//trim(f%file_prefix)//"0000"//trim(f%file_suffix)//".nc"
    ! integer, intent(in) :: nx, ny, nz
    character(*), intent(in)           :: path
    type(t_domain), target, intent(in) :: domain
    character(*), intent(in)           :: dimnames(:) ! Names of the dimensions, assumed to be in the order [x, y, z]
    character(*), intent(in)           :: file_prefix, file_suffix
    type(t_datetime), intent(in)       :: start
    real(rk), intent(in)               :: dt
    character(*), optional, intent(in) :: pmap
    character(len=LEN_CHAR_L)          :: init_path

    f%PATH = path
    f%domain => domain

    if (present(pmap)) then
      f%PMAPFILE = pmap
      f%has_subdomains = .true.
      call f%init_subdomain()
    end if

    f%file_prefix = file_prefix
    f%file_suffix = file_suffix
    call f%init_dirlist()

    select case (f%has_subdomains)
    case (.true.)
      init_path = f%get_directory(1)
    case (.false.)
      init_path = f%get_file(1)
    end select

    call f%init_time(init_path, start, dt)

    call f%init_dimensions(dimnames)

    return
#undef PROC0
#define PROC0 "/"//trim(this%file_prefix)//"0000"//trim(this%file_suffix)//".nc"
  end function ctor_fieldset
  !===========================================
  ! INITIALIZATION FUNCTIONS
  !===========================================
  subroutine init_dimensions(this, dimnames)
    !---------------------------------------------
    ! This initializes the dimensions of the fieldset.
    ! The dimensions are assumed to be in the order [x, y, z].
    ! Horizontal dimensions are taken from the domain and the vertical dimension
    ! is taken from the netCDF file.
    !---------------------------------------------
    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: dimnames(:) ! Names of the dimensions, assumed to be in the order [x, y, z]
    character(len=LEN_CHAR_L)        :: check_path
    character(len=LEN_CHAR_S), allocatable :: dimname_check(:)
    integer, allocatable              :: dimlen_check(:)
    integer :: ndim
    integer :: idim
    integer :: ndim_check

    FMT1, "======== Init dimensions ========"

    if (.not. associated(this%domain)) call throw_error("fieldset :: init_dimensions", "Domain not associated")

    ndim = size(dimnames)
    allocate (this%dims(ndim))
    this%dims = 0
    allocate (this%dimnames(ndim), source=dimnames)
    this%has_vertical = .false.

    ! Get the dimensions from the netCDF file
    if (this%has_subdomains) then
      write (check_path, '(a)') trim(this%current_path)//PROC0
    else
      write (check_path, '(a)') trim(this%current_path)
    end if
    FMT2, "Checking dimensions in file: "//trim(check_path)
    call nc_get_file_dims(check_path, ndims=ndim_check, dimnames=dimname_check, dimlens=dimlen_check)
    FMT2, "Dimensions: ", ndim_check
    if (ndim_check < 2) call throw_error("fieldset :: init_dimensions", "Dataset should have at least 2 dimensions")
    do idim = 1, ndim_check
      FMT3, trim(dimname_check(idim))//"(", dimlen_check(idim), ")"
    end do
    ! Map the dimensions to the fieldset
    do idim = 1, ndim
      if (any(dimname_check == trim(dimnames(idim)))) then
        this%dims(idim) = dimlen_check(findloc(dimname_check, trim(dimnames(idim)), dim=1))
        if (idim == 3) then
          this%has_vertical = .true.
        end if
      else
        ERROR, "Dimension ", trim(dimnames(idim)), " not found in file ", trim(check_path)
        call throw_error("fieldset :: init_dimensions", "Dimension not found in file")
      end if
    end do

    if (this%has_subdomains) then
      ! TODO: This is a dirty bug fix, because the subdomain dimensions are not
      !       the same as the whole domain. This should be fixed in the future.
      this%dims(1) = this%domain%nx
      this%dims(2) = this%domain%ny
    end if
    this%nx = this%dims(1)
    this%ny = this%dims(2)
    if (this%has_vertical) then
      this%nz = this%dims(3)
    end if

    ! Some sanity checks
    if (this%nx .ne. this%domain%nx) call throw_error("fieldset :: init_dimensions", "Horizontal dimension mismatch (nx: "//trim(dimnames(1))//")")
    if (this%ny .ne. this%domain%ny) call throw_error("fieldset :: init_dimensions", "Horizontal dimension mismatch (ny: "//trim(dimnames(2))//")")

    FMT2, "Fieldset dimensions:"
    FMT3, trim(dimnames(1))//"(", this%nx, ")"
    FMT3, trim(dimnames(2))//"(", this%ny, ")"
    if (.not. this%has_vertical) then
      FMT3, "( No vertical dimension )"
    else
      FMT3, trim(dimnames(3))//"(", this%nz, ")"
    end if

    return
  end subroutine init_dimensions
  !===========================================
  subroutine init_time(this, init_path, start, dt)
    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: init_path
    type(t_datetime), intent(in)     :: start
    real(rk), intent(in)             :: dt
    real(rk)                         :: t1, t2
    character(len=LEN_CHAR_L)        :: info

    FMT1, "======== Init time ========"

    if (this%has_subdomains) then
      t1 = nc_read_time_val(trim(init_path)//PROC0, 1)
      t2 = nc_read_time_val(trim(init_path)//PROC0, 2)
    else
      t1 = nc_read_time_val(trim(init_path), 1)
      t2 = nc_read_time_val(trim(init_path), 2)
    end if

    this%nc_timestep = t2 - t1
    if (this%nc_timestep <= ZERO) then
      ERROR, "Fieldset time step <= 0: ", this%nc_timestep, " [s]"
      call throw_error("fieldset :: fieldset", "Bad time step")
    end if

    FMT2, "Fieldset time step:", this%nc_timestep, " [s]"

    call this%set_start_time(start)
    call this%set_simulation_timestep(dt)

    FMT2, "The starting path is "//trim(this%current_path)
    write(info, "(a,i4,a,i4)") "The starting time is "//this%date_t1%nice_format()//" at time step ", this%read_idx, " of ", this%current_ntimes
    FMT2, trim(info)
    ! FMT2, "The starting time is "
    ! call this%date_t1%print_short_date()
    ! FMT2, "at time step ", this%read_idx, " of ", this%current_ntimes

    return
  end subroutine init_time
  !===========================================
  subroutine init_subdomain(this)
    !---------------------------------------------
    ! This maps the pieces of chunked data using par_setup.
    !---------------------------------------------
    class(t_fieldset), intent(inout) :: this
    integer                          :: imax, jmax ! Total size of the domain
    integer                          :: ioff, joff ! Subdomain offset
    integer                          :: pnum, iend, jend
    integer                          :: i, j, iproc

    FMT1, "======== Init subdomains ========"

    if (.not. associated(this%domain)) call throw_error("fieldset :: init_subdomain", "Domain not associated")

    open (PROCFILE, file=trim(this%PMAPFILE), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: init_subdomain", "Failed to open "//trim(this%PMAPFILE), ierr)
    read (PROCFILE, *) this%nproc
    allocate (this%pmap(this%nproc, 2))
    read (PROCFILE, *) this%nxp, this%nyp, imax, jmax
    if ((imax .ne. this%domain%nx) .or. (jmax .ne. this%domain%ny)) then
      ERROR, "(nx = ", this%domain%nx, "imax = ", imax, "ny = ", this%domain%ny, "jmax = ", jmax, ")"
      call throw_error("fieldset :: init_subdomain", "Given domain size does not match one in parallel setup file!")
    else
      FMT2, "Number of subdomains: ", this%nproc, " (", this%nxp, " x ", this%nyp, ")"
    end if
    allocate (this%pmask(imax, jmax))
    do i = 1, imax
      do j = 1, jmax
        this%pmask(i, j) = -10
      end do
    end do
    iproc = 0
    do while (iproc .lt. this%nproc)
      read (PROCFILE, *) pnum, ioff, joff
      i = ioff + 1; j = joff + 1
      iproc = iproc + 1
      this%pmap(iproc, 1) = ioff; this%pmap(iproc, 2) = joff
      iend = this%nxp; jend = this%nyp
      if (i .lt. 1) then
        !---------------------------------------------
        ! If i is less than 1, then ioff must be negative.
        ! Therefore this processor actually covers less area,
        ! i.e. iend must be reduced ioff amount
        iend = this%nxp + ioff ! Same as nxp - abs(ioff)
        i = 1
      end if
      if (j .lt. 1) then
        !---------------------------------------------
        ! Same logic here
        jend = this%nyp + joff
        j = 1
      end if
      this%pmask(i:i + iend, j:j + jend) = pnum
    end do
    close (PROCFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: init_subdomain", "Failed to close "//trim(this%PMAPFILE), ierr)

    FMT2, "Processor mask initialised"

    return
  end subroutine init_subdomain
  !===========================================
  subroutine init_dirlist(this)
    !---------------------------------------------
    ! Get a list of directories or files that contain the data
    ! The files should be named so that ls command would give them
    ! in the right order. Including the date in the file name should be enough.
    ! TODO: Alternatively, if all the files are in one directory,
    !       it should list all the files in this%PATH.
    !       (separate routine e.g., init_filelist ?)
    ! TODO (later, probably never): sort the files somehow so all of this
    !       would not depend on ls getting it right.
    !---------------------------------------------
    class(t_fieldset), intent(inout) :: this
    logical                          :: dirlist_exists
    integer                          :: idir

    FMT1, "======== Init dirlist ========"
    inquire (file=trim(dirinfile), exist=dirlist_exists)
    if (.not. dirlist_exists) then
#ifndef NOSYSCALLS
      select case (this%has_subdomains)
      case (.true.)
        FMT2, "Getting directory list..."
        call system('( ls -d '//trim(this%PATH)//'/*/ | wc -l && ls -d ' &
                    //trim(this%PATH)//'/*/ | xargs -n 1 basename ) > ' &
                    //trim(dirinfile), status=ierr)
        if (ierr .ne. 0) call throw_error("fieldset :: init_dirlist", "Could not get directory list from "//trim(this%PATH), ierr)
      case (.false.)
        FMT2, "Getting file list..."
        call system('( ls '//trim(this%PATH)//'/'//trim(this%file_prefix)//'*' &
                    //trim(this%file_suffix)//'*.nc | wc -l && ls '//trim(this%PATH) &
                    //'/'//trim(this%file_prefix)//'*'//trim(this%file_suffix) &
                    //'*.nc | xargs -n 1 basename ) > '//trim(dirinfile), status=ierr)
        if (ierr .ne. 0) call throw_error("fieldset :: init_dirlist", "Could not get file list from "//trim(this%PATH), ierr)
      end select
#else
      call throw_error("fieldset :: init_dirlist", "File/directory list ("//trim(dirinfile)//") does not exist!")
#endif
    end if

    open (DIRFILE, file=trim(dirinfile), action='read', iostat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: init_dirlist", "Failed to open "//trim(dirinfile), ierr)
    read (DIRFILE, *) this%nentries
    FMT2, "Found ", this%nentries, "files/directories in "//trim(this%PATH)
    select case (this%has_subdomains)
    case (.true.)
      allocate (this%dirlist(this%nentries))
      do idir = 1, this%nentries
        read (DIRFILE, *) this%dirlist(idir)
      end do
    case (.false.)
      allocate (this%filelist(this%nentries))
      do idir = 1, this%nentries
        read (DIRFILE, '(a)') this%filelist(idir)
      end do
    end select
    close (DIRFILE, iostat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: init_dirlist", "Failed to close "//trim(dirinfile), ierr)

    FMT2, "Dirlist initialized"

  end subroutine init_dirlist
  !===========================================
  subroutine add_field(this, field_name, nc_varname)
    class(t_fieldset), intent(inout)       :: this
    character(*), intent(in)               :: field_name, nc_varname
    integer                                :: fld_ndim
    integer                                :: fld_time_idx
    logical                                :: fld_static
    character(len=LEN_CHAR_S), allocatable :: fld_dimnames(:)
    integer, allocatable                   :: fld_dimlen(:)
    integer                                :: idim
    character(len=LEN_CHAR_L)              :: check_path
    real(rk)                               :: fill_value
    character(len=LEN_CHAR_L)              :: errmsg

    fld_static = .true.
    fld_time_idx = -1

    ! Check if the input arguments match the actual field dimensions
    if (this%has_subdomains) then
      write (check_path, '(a)') trim(this%current_path)//PROC0
    else
      write (check_path, '(a)') trim(this%current_path)
    end if

    if (.not. nc_var_exists(check_path, nc_varname)) call throw_error("fieldset :: add_field", "Variable "//trim(nc_varname)//" does not exist in "//trim(check_path))

    if (this%has_subdomains) then
      call nc_get_var_dims(check_path, nc_varname, ndims=fld_ndim, dimnames=fld_dimnames)
    else
      call nc_get_var_dims(check_path, nc_varname, ndims=fld_ndim, dimnames=fld_dimnames, dimlens=fld_dimlen)
    end if

    if (any(fld_dimnames == "time")) then
      fld_static = .false.
      fld_time_idx = findloc(fld_dimnames, "time", dim=1)
    end if

    do idim = 1, fld_ndim
      if (idim == fld_time_idx) cycle
      if (.not. any(this%dimnames == fld_dimnames(idim))) then
        write (errmsg, '(a)') "Dimension "//trim(fld_dimnames(idim))//" of variable "//trim(nc_varname)//" does not exist in "//trim(check_path)
        call throw_error("fieldset :: add_field", errmsg)
      end if
    end do

    if (.not. fld_static) then
      fld_ndim = fld_ndim - 1 ! Remove the time dimension
    end if

    if (this%has_subdomains) then
      allocate (fld_dimlen(fld_ndim))
      do idim = 1, fld_ndim
        fld_dimlen(idim) = this%dims(idim)
      end do
    end if

    ! Search for _FillValue attribute
    if (.not. nc_get_var_fillvalue(check_path, nc_varname, fill_value)) then
      call throw_warning("fieldset :: add_field", "No _FillValue attribute found for variable "//trim(nc_varname)//" in "//trim(check_path)//". Using default.")
      fill_value = FILLVALUE_BIG
    end if

    ! Add the field to the fieldset
    select case (fld_static)
    case (.true.)
      select case (fld_ndim)
      case (1)
        call this%fields%add_node(field_name, &
                                  t_field_static(n=fld_dimlen(1), &
                                                 name=nc_varname, &
                                                 fill_value=fill_value))
      case (2)
        call this%fields%add_node(field_name, &
                                  t_field_static(n1=fld_dimlen(1), &
                                                 n2=fld_dimlen(2), &
                                                 name=nc_varname, &
                                                 fill_value=fill_value))
      case (3)
        call this%fields%add_node(field_name, &
                                  t_field_static(n1=fld_dimlen(1), &
                                                 n2=fld_dimlen(2), &
                                                 n3=fld_dimlen(3), &
                                                 name=nc_varname, &
                                                 fill_value=fill_value))
      end select
    case (.false.)
      select case (fld_ndim)
      case (1)
        call this%fields%add_node(field_name, &
                                  t_field_dynamic(n=fld_dimlen(1), &
                                                  timestep=this%nc_timestep, &
                                                  name=nc_varname, &
                                                  fill_value=fill_value))
      case (2)
        call this%fields%add_node(field_name, &
                                  t_field_dynamic(n1=fld_dimlen(1), &
                                                  n2=fld_dimlen(2), &
                                                  timestep=this%nc_timestep, &
                                                  name=nc_varname, &
                                                  fill_value=fill_value))
      case (3)
        call this%fields%add_node(field_name, &
                                  t_field_dynamic(n1=fld_dimlen(1), &
                                                  n2=fld_dimlen(2), &
                                                  n3=fld_dimlen(3), &
                                                  timestep=this%nc_timestep, &
                                                  name=nc_varname, &
                                                  fill_value=fill_value))
      end select
    end select

    this%num_fields = this%num_fields + 1

    return
  end subroutine add_field
  !===========================================
  subroutine init_u_mask(this)
    !---------------------------------------------
    ! Make the velocity u component (similarly for v comp. in init_v_mask)
    ! to deal with boundaries.
    ! SEA will be 0, BEACH areas will be 1 or -1 depending on where the land is.
    ! Adjacent land cells will have the same velocity tangential to land (0 gradient).
    ! Normal (cross-boundary) components will be 0.
    !---------------------------------------------
    class(t_fieldset), intent(inout) :: this
    integer                          :: i, j
    integer                          :: seamask(this%nx, this%ny)

    allocate (this%u_mask(this%nx, this%ny), stat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: init_u_mask", "Could not allocate u component mask")

    seamask = this%domain%get_seamask()

    this%u_mask = 0
    do i = 1, this%nx
      do j = 1, this%ny
        if (seamask(i, j) == DOM_BEACH) then
          if (seamask(i, j + 1) == DOM_LAND) then
            this%u_mask(i, j) = -1
          else if (seamask(i, j - 1) == DOM_LAND) then
            this%u_mask(i, j) = 1
          end if
        end if
      end do
    end do

  end subroutine init_u_mask
  !===========================================
  subroutine init_v_mask(this)
    class(t_fieldset), intent(inout) :: this
    integer                          :: i, j
    integer                          :: seamask(this%nx, this%ny)

    allocate (this%v_mask(this%nx, this%ny), stat=ierr)
    if (ierr .ne. 0) call throw_error("fieldset :: init_v_mask", "Could not allocate v component mask")

    seamask = this%domain%get_seamask()

    this%v_mask = 0
    do i = 1, this%nx
      do j = 1, this%ny
        if (seamask(i, j) == DOM_BEACH) then
          if (seamask(i + 1, j) == DOM_LAND) then
            this%v_mask(i, j) = -1
          else if (seamask(i - 1, j) == DOM_LAND) then
            this%v_mask(i, j) = 1
          end if
        end if
      end do
    end do

  end subroutine init_v_mask
  !===========================================
  ! GETTERS
  !===========================================
  real(rk) function get_value_key_real_idx(this, field_name, t, i, j, k) result(res)
    !---------------------------------------------
    ! Get the value of the field at the given time and position
    ! The field can be 2D or 3D, but has to be dynamic (has time dimension).
    !---------------------------------------------
    class(t_fieldset), intent(in)   :: this
    character(*), intent(in)        :: field_name
    real(rk), intent(in)            :: t, i, j
    real(rk), optional, intent(in)  :: k
    class(t_variable), pointer      :: p_field

    res = ZERO

    call this%fields%get_item(field_name, p_field)

    select type (p_field)
    type is (t_field_dynamic_2d)
      res = p_field%get(t=t, x=i, y=j)
    type is (t_field_dynamic_3d)
      if (.not. present(k)) then
        if (this%nz == 1) then
          res = p_field%get(t=t, x=i, y=j, z=ONE)
        else
          call throw_error("fieldset :: get_value_key_real_idx", "k index is missing for 3D field")
        end if
      else
        res = p_field%get(t=t, x=i, y=j, z=k)
      end if
    class default
      call throw_error("fieldset :: get_value_key_real_idx", "The field must be dynamic 2D or 3D field")
    end select

    return
  end function get_value_key_real_idx
  !===========================================
  real(rk) function get_value_key_int_idx(this, field_name, t, i, j, k) result(res)
    class(t_fieldset), intent(in)  :: this
    character(*), intent(in)       :: field_name
    real(rk), intent(in)           :: t
    integer, intent(in)            :: i, j
    integer, optional, intent(in)  :: k
    class(t_variable), pointer     :: p_field

    res = ZERO

    call this%fields%get_item(field_name, p_field)

    select type (p_field)
    type is (t_field_dynamic_2d)
      res = p_field%get(t=t, i=i, j=j)
    type is (t_field_dynamic_3d)
      if (.not. present(k)) then
        if (this%nz == 1) then
          res = p_field%get(t=t, i=i, j=j, k=1)
        else
          call throw_error("fieldset :: get_value_key_int_idx", "k index is missing for 3D field")
        end if
      else
        res = p_field%get(t=t, i=i, j=j, k=k)
      end if
    class default
      call throw_error("fieldset :: get_value_key_int_idx", "The field must be dynamic 2D or 3D field")
    end select

    return
  end function get_value_key_int_idx
  !===========================================
  real(rk) function get_value_idx_real_idx(this, idx, t, i, j, k) result(res)
    class(t_fieldset), intent(in)  :: this
    integer, intent(in)            :: idx
    real(rk), intent(in)           :: t, i, j
    real(rk), optional, intent(in) :: k
    class(t_variable), pointer     :: p_field

    res = ZERO

    call this%fields%get_item(idx, p_field)

    select type (p_field)
    type is (t_field_dynamic_2d)
      res = p_field%get(t=t, x=i, y=j)
    type is (t_field_dynamic_3d)
      if (.not. present(k)) then
        if (this%nz == 1) then
          res = p_field%get(t=t, x=i, y=j, z=ONE)
        else
          call throw_error("fieldset :: get_value_idx_real_idx", "k index is missing for 3D field")
        end if
      else
        res = p_field%get(t=t, x=i, y=j, z=k)
      end if
    class default
      call throw_error("fieldset :: get_value_idx_real_idx", "The field must be dynamic 2D or 3D field")
    end select

    return
  end function get_value_idx_real_idx
  !===========================================
  real(rk) function get_value_idx_int_idx(this, idx, t, i, j, k) result(res)
    class(t_fieldset), intent(in)  :: this
    integer, intent(in)            :: idx
    real(rk), intent(in)           :: t
    integer, intent(in)            :: i, j
    integer, optional, intent(in)  :: k
    class(t_variable), pointer     :: p_field

    res = ZERO

    call this%fields%get_item(idx, p_field)

    select type (p_field)
    type is (t_field_dynamic_2d)
      res = p_field%get(t=t, i=i, j=j)
    type is (t_field_dynamic_3d)
      if (.not. present(k)) then
        if (this%nz == 1) then
          res = p_field%get(t=t, i=i, j=j, k=1)
        else
          call throw_error("fieldset :: get_value_idx_int_idx", "k index is missing for 3D field")
        end if
      else
        res = p_field%get(t=t, i=i, j=j, k=k)
      end if
    class default
      call throw_error("fieldset :: get_value_idx_int_idx", "The field must be dynamic 2D or 3D field")
    end select

    return
  end function get_value_idx_int_idx
  !===========================================
  integer function get_field_index(this, field_name) result(res)
    class(t_fieldset), intent(in) :: this
    character(*), intent(in)      :: field_name

    res = this%fields%node_loc(field_name)

    return
  end function get_field_index
  !===========================================
  real(rk) function get_time(this, date) result(res)
    class(t_fieldset), intent(in) :: this
    type(t_datetime), intent(in)  :: date

    res = date_diff(this%date_t1, date)

  end function get_time
  !===========================================
  character(len=LEN_CHAR_L) function get_directory(this, idx) result(res)
    class(t_fieldset), intent(inout) :: this
    integer, intent(in)           :: idx

    if (idx < 1) call throw_error("fieldset :: get_directory", "Index out of bounds (idx < 1)")
    if (idx > this%nentries) then
       ! call throw_error("fieldset :: get_directory", "Index out of bounds (idx > number of folders)")
       ! If the index is larger than the number of entries then it PROBABLY means that we 
       ! are at the end of the dataset and we should not try to read anymore. If fieldset::update is 
       ! called again, then an error should be thrown.
       this%has_more = .false.
       res = this%current_path
       return
    endif 
    write (res, '(a,i0.8)') trim(this%PATH)//'/', this%dirlist(idx)

    return
  end function get_directory
  !===========================================
  character(len=LEN_CHAR_L) function get_file(this, idx) result(res)
    class(t_fieldset), intent(inout) :: this
    integer, intent(in)           :: idx

    if (idx < 1) call throw_error("fieldset :: get_file", "Index out of bounds (idx < 1)")
    if (idx > this%nentries) then 
      ! call throw_error("fieldset :: get_file", "Index out of bounds (idx > number of files)")
      this%has_more = .false.
      res = this%current_path
      return
    end if
    write (res, '(a)') trim(this%PATH)//'/'//trim(this%filelist(idx))

    return
  end function get_file
  !===========================================
  subroutine get_indices_vertical(this, t, z, i, j, k, kr)
    !---------------------------------------------
    ! Get the vertical indices for a given depth.
    ! Vertical index search is in this module because
    ! adaptive vertical (sigma) coordinates change in time
    ! and need to be read from the netcdf file. The domain class
    ! does not handle this.
    ! TODO: zax direction should be taken into account
    !---------------------------------------------
    class(t_fieldset), intent(in)   :: this
    real(rk), intent(in)            :: t, z
    integer, intent(in)             :: i, j
    integer, optional, intent(out)  :: k
    real(rk), optional, intent(out) :: kr
    real(rk)                        :: zax(this%nz)
    real(rk)                        :: dep
    real(rk)                        :: dz
    integer                         :: ik
    real(rk)                        :: ikr

    dep = this%domain%get_bathymetry(i, j)

    ! If the bathymetry value is negative, that means that this (i,j) point is on land
    ! and we set the index to the top of the water column and return
    if (isnan(dep)) then
      if (present(k)) k = this%zax_top_idx
      if (present(kr)) kr = real(this%zax_top_idx, rk)
#ifdef DEBUG
      call throw_warning("fieldset :: get_indices_vertical", "Particle on ground")
#endif
      return
    end if

    ! Get the Z axis
    select case (this%zax_style)
    case (STATIC_DEPTH_VALUES)
      zax = this%get_zax() ! No input parameters indicates that the array is 1D and we want the whole thing
    case (DEPTH_VALUES, LAYER_THICKNESS)
      zax = this%get_zax(t, i, j) ! Don't want to interpolate the Z axis in space, so we're taking the closest indices
    end select

    ! If the depth is greater than the top of the water column, that means that the particle is
    ! above the top of the water column (has "jumped" out of water) and we set the index to the top of the water column
    ! and return. We must also include the case where the depth is equal to the top of the water column, because
    ! otherwise the index will not be decreased and will cause an error.
    if (z >= zax(this%zax_top_idx)) then
      if (present(k)) k = this%zax_top_idx
      if (present(kr)) kr = real(this%zax_top_idx, rk)
#ifdef DEBUG
      call throw_warning("fieldset :: get_indices_vertical", "Out of bounds (z > top)") ! Better error message
#endif
      return
    end if

    if (z <= zax(this%zax_bot_idx)) then
      if (present(k)) k = this%zax_bot_idx
      if (present(kr)) kr = real(this%zax_bot_idx, rk)
#ifdef DEBUG
      call throw_warning("fieldset :: get_indices_vertical", "Out of bounds (z < bottom)") ! Better error message
#endif
      return
    end if

    ! Because of the previous checks, we know that the particle is in the water column
    ! and actually in the range of the Z axis. We can now find the indices.
    ik = minloc(abs(zax - z), dim=1)
    if (zax(ik) > z) then
      ik = ik - 1
    end if
    dz = (zax(ik + 1) - zax(ik))
    do while (dz <= ZERO) ! It may happen that there are equal values in the Z axis
      ik = ik + 1
      dz = (zax(ik + 1) - zax(ik))
    end do
    ikr = ik + (z - zax(ik)) / dz

    if (present(k)) then
      k = ik
    end if
    if (present(kr)) then
      kr = ikr
    end if

    return
  end subroutine get_indices_vertical
  !===========================================
  function get_zax_regular(this) result(res)
    class(t_fieldset), intent(in) :: this
    class(t_variable), pointer    :: p_field
    real(rk), allocatable         :: res(:)

    if (this%zax_style /= STATIC_DEPTH_VALUES) call throw_error("fieldset :: get_zax_regular", "Z axis is not static")

    res = [ZERO]
    if (this%zax_idx == 0) return

    call this%fields%get_item(this%zax_idx, p_field)
    select type (p_field)
    type is (t_field_static_1d)
      res = p_field%get()
      res = res * real(this%zax_dir, rk)
    class default
      call throw_error("fieldset :: get_zax_regular", "Z axis is not a static 1D field")
    end select

    return
  end function get_zax_regular
  !===========================================
  function get_zax_adaptive(this, t, i, j) result(res)
    class(t_fieldset), intent(in)    :: this
    real(rk), intent(in)             :: t
    integer, intent(in)              :: i, j
    class(t_variable), pointer       :: p_field
    real(rk)                         :: arr_zax(this%nz)
    real(rk)                         :: res(this%nz)
    integer                          :: ik

    res = [ZERO]
    if (this%zax_idx == 0) return

    call this%fields%get_item(this%zax_idx, p_field)
    select type (p_field)
    type is (t_field_dynamic_3d)
      call p_field%slice(slice_dim=3, t=t, idx_other=[i, j], res=arr_zax)
    class default
      call throw_error("fieldset :: get_zax_adaptive", "Z axis is not a dynamic 3D field")
    end select

    select case (this%zax_style)
    case (DEPTH_VALUES)
      res = arr_zax * real(this%zax_dir, rk) ! If zax_dir is positive down (-1), then we want to multiply by -1 because particles' coordinates are positive up
    case (LAYER_THICKNESS)
      res = arr_zax
      res(1) = -ONE * this%domain%get_bathymetry(i, j) ! Bathymetry is positive down
      do ik = 2, this%zax_top_idx
        res(ik) = res(ik - 1) + arr_zax(ik)
      end do
    case default
      call throw_error("fieldset :: get_zax_adaptive", "Z - axis style unknown: zax_style = 1 or 2")
    end select

    return
  end function get_zax_adaptive
  !===========================================
  subroutine search_indices(this, t, x, y, z, i, j, k, ir, jr, kr)
    class(t_fieldset), intent(in)   :: this
    real(rk), optional, intent(in)  :: t           ! time
    real(rk), intent(in)            :: x, y        ! and position (lon, lat, depth)
    real(rk), optional, intent(in)  :: z
    integer, optional, intent(out)  :: i, j, k     ! integer indices out
    real(rk), optional, intent(out) :: ir, jr, kr  ! real indices out
    integer                         :: ii, jj, kk
    real(rk)                        :: iir, jjr, kkr

    call this%domain%get_index(x, i=ii, ir=iir, dim=1)
    call this%domain%get_index(y, i=jj, ir=jjr, dim=2)
    if (present(k) .or. present(kr)) then
      if (.not. (present(t) .and. present(z))) call throw_error("fieldset :: search_indices", "Time or depth missing!")
      if (this%zax_idx < 1) then ! Probably the fastest way to make sure the Z axis exists
        kk = this%nz; kkr = real(this%nz, kind=rk)
      else
        call this%get_indices_vertical(t, z, ii, jj, k=kk, kr=kkr)
      end if
    end if
    if (present(i)) then
      i = ii; 
    end if
    if (present(j)) then
      j = jj; 
    end if
    if (present(k)) then
      k = kk; 
    end if
    if (present(ir)) then
      ir = iir; 
    end if
    if (present(jr)) then
      jr = jjr; 
    end if
    if (present(kr)) then
      kr = kkr; 
    end if

    return
  end subroutine search_indices
  !===========================================
  subroutine find_directory(this, date, dir_name, dir_idx)
    !---------------------------------------------
    ! TODO: might have to check if nc files start at time
    ! [date] 00:00:00 or [date] 00:00:10 (at least in this particular example)
    ! Compare every time?
    ! Some validity check would be nice
    !---------------------------------------------

    class(t_fieldset), intent(in)          :: this
    type(t_datetime), intent(in)           :: date
    character(len=LEN_CHAR_L), intent(out) :: dir_name
    integer, optional, intent(out)         :: dir_idx
    integer                                :: i
    integer(rk)                            :: YYYYMMDD
    character(len=LEN_CHAR_L)              :: testdir

    if (this%nentries == 1) then
      write (testdir, '(a,i0.8)') trim(this%PATH)//'/', this%dirlist(1)
      if (date < datetime_from_netcdf(trim(testdir)//PROC0, 1)) then
        call throw_error("fieldset :: find_directory", "The date is earlier than first date in data file")
      end if
      dir_name = testdir
      if (present(dir_idx)) dir_idx = 1
      return
    end if

    YYYYMMDD = date%short_format(.false.)
    do i = 1, this%nentries
      if (this%dirlist(i) .ge. YYYYMMDD) then
        write (testdir, '(a,i0.8)') trim(this%PATH)//'/', this%dirlist(i)
        if (date < datetime_from_netcdf(trim(testdir)//PROC0, 1)) then
          if (i == 1) call throw_error("fieldset :: find_directory", "The date is earlier than first date in data file")
          write (dir_name, '(a,i0.8)') trim(this%PATH)//'/', this%dirlist(i - 1)
          if (present(dir_idx)) dir_idx = i - 1
        else
          write (dir_name, '(a,i0.8)') trim(this%PATH)//'/', this%dirlist(i)
          if (present(dir_idx)) dir_idx = i
        end if
        return
      end if
    end do

    call throw_error("fieldset :: find_directory", "Did not find folder")

    return
  end subroutine find_directory
  !===========================================
  subroutine find_file(this, date, file_name, file_idx)
    !---------------------------------------------
    ! Assumes a sorted filelist.
    !---------------------------------------------

    class(t_fieldset), intent(in)          :: this
    type(t_datetime), intent(in)           :: date
    character(len=LEN_CHAR_L), intent(out) :: file_name
    integer, optional, intent(out)         :: file_idx
    integer                                :: i, n_times

    if (this%nentries == 1) then
      write (file_name, '(a)') trim(this%PATH)//'/'//trim(this%filelist(1))
      if ((date < datetime_from_netcdf(trim(file_name), 1))) then
        call throw_error("fieldset :: find_file", "The date is earlier than first date in data file")
      end if
      if (present(file_idx)) file_idx = 1
      return
    end if

    do i = 1, this%nentries
      if (date < datetime_from_netcdf(trim(this%PATH)//'/'//trim(this%filelist(i)), 1)) then
        if (i == 1) call throw_error("fieldset :: find_file", "The date is earlier than first date in data file")
        write (file_name, '(a)') trim(this%PATH)//'/'//trim(this%filelist(i - 1))
        if (present(file_idx)) file_idx = i - 1
        return
      end if
    end do
    ! If we reach here, just return the last file and hope for the best
    write (file_name, '(a)') trim(this%PATH)//'/'//trim(this%filelist(this%nentries))
    if (date > datetime_from_netcdf(trim(file_name), n_times)) call throw_error("fieldset :: find_file", "Did not find file")
    if (present(file_idx)) file_idx = this%nentries

    return
  end subroutine find_file
  !===========================================
  subroutine get_pmap(this, pmapout)
    !---------------------------------------------
    ! Just in case I want to use pmap somewhere else
    ! (it's private otherwise)
    !---------------------------------------------
    class(t_fieldset), intent(in) :: this
    integer, intent(inout)        :: pmapout(this%nproc, 2)

    pmapout = this%pmap

    return
  end subroutine get_pmap
  !===========================================
  ! SETTERS
  !===========================================
  subroutine set_value_key(this, field_name, data)
    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: field_name
    real(rk), intent(in)             :: data(this%nx, this%ny, this%nz)
    class(t_variable), pointer       :: p_field

    call this%fields%get_item(field_name, p_field)
    select type (p_field)
    type is (t_field_static_3d)
      call p_field%set(data)
    type is (t_field_dynamic_3d)
      call p_field%set(data)
    class default
      call throw_error("fieldset :: set_value_key", "The field must be 3D")
    end select

    return
  end subroutine set_value_key
  !===========================================
  subroutine set_value_idx(this, idx, data)
    class(t_fieldset), intent(inout) :: this
    integer, intent(in)              :: idx
    real(rk), intent(in)             :: data(this%nx, this%ny, this%nz)
    class(t_variable), pointer       :: p_field

    call this%fields%get_item(idx, p_field)
    select type (p_field)
    type is (t_field_static_3d)
      call p_field%set(data)
    type is (t_field_dynamic_3d)
      call p_field%set(data)
    class default
      call throw_error("fieldset :: set_value_idx", "The field must be 3D")
    end select

    return
  end subroutine set_value_idx
  !===========================================
  subroutine set_u_component(this, u_comp_name)
    !---------------------------------------------
    ! Set the index of the velocity u component and create
    ! the mask (similarly for v component in set_v_component).
    !---------------------------------------------

    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: u_comp_name

    if (.not. this%has_field(trim(u_comp_name))) call throw_error("fieldset :: set_u_component", "Did not find "//trim(u_comp_name)//" in fieldset")
    call this%init_u_mask()
    this%u_idx = this%fields%node_loc(u_comp_name)

  end subroutine set_u_component
  !===========================================
  subroutine set_v_component(this, v_comp_name)
    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: v_comp_name

    if (.not. this%has_field(trim(v_comp_name))) call throw_error("fieldset :: set_v_component", "Did not find "//trim(v_comp_name)//" in fieldset")
    call this%init_v_mask()
    this%v_idx = this%fields%node_loc(v_comp_name)

  end subroutine set_v_component
  !===========================================
  subroutine set_zax(this, zax_name, zax_style, zax_direction)
    class(t_fieldset), intent(inout) :: this
    character(*), intent(in)         :: zax_name
    integer, intent(in)              :: zax_style
    integer, intent(in)              :: zax_direction
    class(t_variable), pointer       :: p_field
    real(rk), allocatable            :: buffer(:)
    character(len=LEN_CHAR_L)        :: read_path

 if (.not. this%has_vertical) call throw_error("fieldset :: set_zax", "The fieldset does not have a vertical axis. Cannot set zax.")

  if (.not. this%has_field(trim(zax_name))) call throw_error("fieldset :: set_zax", "Did not find "//trim(zax_name)//" in fieldset")

    this%zax_bot_idx = 1        ! This is probably not necessary as the layers in the fieldset should always be ordered from bottom to top
    this%zax_top_idx = this%nz
    this%zax_style = zax_style
    if (zax_direction > 0) then
      this%zax_dir = 1
    else
      this%zax_dir = -1
    end if
    this%zax_idx = this%fields%node_loc(trim(zax_name))

    call this%fields%get_item(this%zax_idx, p_field)
    select type (p_field)
    type is (t_field_static_1d)
      if (this%has_subdomains) then
        write (read_path, '(a)') trim(this%current_path)//PROC0
      else
        write (read_path, '(a)') trim(this%current_path)
      end if
      allocate (buffer(this%nz))
      call nc_read_real_1d(read_path, p_field%get_name(), this%nz, buffer)
      call p_field%set(buffer)
    type is (t_field_dynamic_3d)
      ! ! Read later
    class default
      call throw_error("fieldset :: set_zax", "The zax field must be 1D static or 3D dynamic")
    end select

    return
  end subroutine set_zax
  !===========================================
  subroutine set_start_time(this, date)
    !---------------------------------------------
    ! Set the start time of the fieldset.
    ! date0 is the time of the first time step in the first file.
    ! read_idx is the index of the time step closest to date (before the current date).
    ! date1 is the next time step after date0 in the file.
    ! next_read_dt is also the date of the next time step, this is when the next time step will be read.
    !---------------------------------------------
    class(t_fieldset), intent(inout) :: this
    type(t_datetime), intent(in)     :: date
    type(t_datetime)                 :: date0, date1

    if (this%has_subdomains) then
      call this%find_directory(date, this%current_path, this%dirlist_idx)
      call nc_get_dim_len(trim(this%current_path)//PROC0, 'time', this%current_ntimes)
      date0 = datetime_from_netcdf(trim(this%current_path)//PROC0, n=1)
      this%read_idx = int(date_diff(date0, date) / this%nc_timestep) + 1 ! +1 because index starts at 1
      date1 = datetime_from_netcdf(trim(this%current_path)//PROC0, n=this%read_idx)
    else
      call this%find_file(date, this%current_path, this%dirlist_idx)
      call nc_get_dim_len(trim(this%current_path), 'time', this%current_ntimes)
      date0 = datetime_from_netcdf(trim(this%current_path), n=1)
      this%read_idx = int(date_diff(date0, date) / this%nc_timestep) + 1
      date1 = datetime_from_netcdf(trim(this%current_path), n=this%read_idx)
    end if
    this%date_t1 = date1
    this%date_t2 = date1%nextDate(this%nc_timestep)
    this%next_read_dt = this%date_t2

    return
  end subroutine set_start_time
  !===========================================
  subroutine set_simulation_timestep(this, dt)
    class(t_fieldset), intent(inout) :: this
    real(rk), intent(in)             :: dt

    if (dt <= this%nc_timestep) then
      FMT2, "Simulation time step smaller than netCDF time step. Reading every time step."
      return
    else
      if (mod(dt, this%nc_timestep) .ne. 0) then
        ERROR, "Time step ", dt, " [s] not divisible by netCDF time step ", this%nc_timestep, " [s]"
        call throw_error("fieldset :: set_simulation_timestep", "Time step should be divisible by netCDF time step")
      end if
      this%read_idx_increment = int(dt / this%nc_timestep)
      FMT2, "Reading every ", this%read_idx_increment, " time steps"
    end if

    return
  end subroutine set_simulation_timestep
  !===========================================
  ! READ AND UPDATE
  !===========================================
  subroutine update(this, date, ignore_check, update_dates)
    class(t_fieldset), intent(inout) :: this
    class(t_datetime), intent(in)    :: date
    logical, optional, intent(in)    :: ignore_check, update_dates
    class(t_variable), pointer       :: p_field
    integer                          :: i_field
    logical                          :: ign_chk, ud

    if (present(ignore_check)) then
      ign_chk = ignore_check
    else
      ign_chk = .false.
    end if

    ! Check if it's even time to readk
    if ((.not. ign_chk) .and. (date < this%next_read_dt)) then
      return
    end if

    if (.not. this%has_more) then
      call throw_error("fieldset :: update", "Trying to update when there are no more data to be read. Ensure that the dataset covers the simulation period.")
    end if

#ifdef DEBUG
    DBG, "Updating fieldset"
    call date%print_short_date()
#endif

    do i_field = 1, this%num_fields
      call this%fields%get_item(i_field, p_field)
      select type (p_field)
      class is (t_field_dynamic)
        call p_field%swap()
        if (this%has_subdomains) then
          call this%read_field_subdomains(p_field, i_field)
        else
          call this%read_field(p_field, i_field)
        end if
      class default
        ! Do nothing if the field is static
        ! call throw_error("fieldset :: update", "Field "//trim(p_field%get_name())//" is not dynamic")
      end select
    end do

    this%read_idx = this%read_idx + this%read_idx_increment
    if (this%read_idx > this%current_ntimes) then
      this%read_idx = this%read_idx - this%current_ntimes
      ! Hopefully noone will be skipping whole files
      this%dirlist_idx = this%dirlist_idx + 1
      if (this%has_subdomains) then
        this%current_path = this%get_directory(this%dirlist_idx)
        call nc_get_dim_len(trim(this%current_path)//PROC0, 'time', this%current_ntimes)
      else
        this%current_path = this%get_file(this%dirlist_idx)
        call nc_get_dim_len(trim(this%current_path), 'time', this%current_ntimes)
      end if
    end if

    if (present(update_dates)) then
      ud = update_dates
    else
      ud = .true.
    end if

    if (.not. ign_chk) call this%next_read_dt%update(this%nc_timestep)
    if (ud) then
      call this%date_t1%update(this%nc_timestep)
      call this%date_t2%update(this%nc_timestep)
    end if

    return
  end subroutine update
  !===========================================
  subroutine read_first_timesteps(this, date)
    !---------------------------------------------
    ! Read the first two timesteps of the fieldset
    ! The dates are not updated and simulation time is not compared
    ! to next_read_dt.
    ! read_idx should be incremented twice after this.
    !---------------------------------------------
    class(t_fieldset), intent(inout) :: this
    type(t_datetime), intent(in)     :: date

    call this%update(date, ignore_check=.true., update_dates=.false.)
    call this%update(date, ignore_check=.true., update_dates=.false.)

    return
  end subroutine read_first_timesteps
  !===========================================
  subroutine read_field(this, p_field, field_idx)
    class(t_fieldset), intent(in)             :: this
    class(t_field_dynamic), target, intent(inout) :: p_field
    integer, intent(in)                       :: field_idx
    character                                 :: c_field
    character(len=LEN_CHAR_S)                 :: varname
    integer                                   :: n_dims
    real(rk), dimension(:, :, :), allocatable :: buffer
    integer, allocatable                      :: start(:), count(:)
    integer                                   :: i, j
    integer                                   :: seamask(this%nx, this%ny)
    real(rk)                                  :: missing_value
    logical                                   :: t_nan, b_nan

    t_nan = .false.
    b_nan = .false.
    n_dims = 0
    missing_value = FILLVALUE_BIG

    if (field_idx == this%u_idx) then
      c_field = "u"
    else if (field_idx == this%v_idx) then
      c_field = "v"
    else
      c_field = "x"
    end if

    select type (p_field)
    type is (t_field_dynamic_2d)
      varname = p_field%get_name()
      n_dims = p_field%get_dim()
      missing_value = p_field%get_missing_value()
    type is (t_field_dynamic_3d)
      varname = p_field%get_name()
      n_dims = p_field%get_dim()
      missing_value = p_field%get_missing_value()
    class default
      call throw_error("fieldset :: read_field", "Field "//trim(p_field%get_name())//" is not 2D or 3D")
    end select

#ifdef DEBUG
    DBG, "Reading "//trim(varname)
#endif

    if (n_dims == 2) then
      start = [1, 1, this%read_idx]
      count = [this%nx, this%ny, 1]
    else if (n_dims == 3) then
      start = [1, 1, 1, this%read_idx]
      count = [this%nx, this%ny, this%nz, 1]
    else
      call throw_error("fieldset :: read_field", "Wrong number of dimensions: "//trim(varname))
    end if

    allocate (buffer(count(1), count(2), count(3)))

    debug(start)
    debug(count)

    if (n_dims == 2) then
      call nc_read_real_3d(trim(this%current_path), trim(varname), start, count, buffer)
    else
      call nc_read_real_4d(trim(this%current_path), trim(varname), start, count, buffer)
    end if

    seamask = this%domain%get_seamask()

    select case (c_field)
    case ("v")
      do j = 1, this%ny
        do i = 1, this%nx
          if (seamask(i, j) == DOM_LAND) then
            buffer(i, j, :) = ZERO
          end if
          ! This is needed due to the A grid
          if (this%u_mask(i, j) > 0) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%v_mask(i, j) > 0) then
            buffer(i, j, :) = buffer(i + 1, j, :)
          else if (this%v_mask(i, j) < 0) then
            buffer(i + 1, j, :) = buffer(i, j, :)
          end if
        end do
      end do
    case ("u")
      do j = 1, this%ny
        do i = 1, this%nx
          if (seamask(i, j) == DOM_LAND) then
            buffer(i, j, :) = ZERO
          end if
          if (this%v_mask(i, j) > 0) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%u_mask(i, j) > 0) then
            buffer(i, j, :) = buffer(i, j + 1, :)
          else if (this%u_mask(i, j) < 0) then
            buffer(i, j + 1, :) = buffer(i, j, :)
          end if
        end do
      end do
    case default
      do j = 1, this%ny
        do i = 1, this%nx
          if (seamask(i, j) == DOM_LAND) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
    end select

    where (buffer == missing_value) buffer = ZERO

    ! It only makes sense to check for NaNs in 3D fields if the
    ! vertical coordinates are following the bathymetry (adaptive).
    if (n_dims == 3 .and. this%zax_style /= STATIC_DEPTH_VALUES) then
      t_nan = all(buffer(:, :, this%zax_top_idx) == ZERO)
      b_nan = all(buffer(:, :, this%zax_bot_idx) == ZERO)
    end if

    select type (p_field)
    type is (t_field_dynamic_2d)
      call p_field%set(buffer(:, :, 1))
    type is (t_field_dynamic_3d)
      call p_field%set(buffer, t_nan, b_nan)
    end select

  end subroutine read_field
  !===========================================
  subroutine read_field_subdomains(this, p_field, field_idx)
    class(t_fieldset), intent(in)                 :: this
    class(t_field_dynamic), target, intent(inout) :: p_field
    integer, intent(in)                           :: field_idx
    character                                     :: c_field
    character(len=LEN_CHAR_S)                     :: varname
    character(len=LEN_CHAR_L)                     :: subdom_filename
    integer                                       :: n_dims
    real(rk), dimension(:, :, :), allocatable     :: buffer, buffer_subdom
    integer, allocatable                          :: start(:), count(:)
    integer                                       :: i, j, i_subdom, ioff, joff, istart, jstart
    integer                                       :: seamask(this%nx, this%ny)
    real(rk)                                      :: missing_value
    logical                                       :: t_nan, b_nan

    t_nan = .false.
    b_nan = .false.
    n_dims = 0
    missing_value = FILLVALUE_BIG

    if (field_idx == this%u_idx) then
      c_field = "u"
    else if (field_idx == this%v_idx) then
      c_field = "v"
    else
      c_field = "x"
    end if

    select type (p_field)
    type is (t_field_dynamic_2d)
      varname = p_field%get_name()
      n_dims = p_field%get_dim()
      missing_value = p_field%get_missing_value()
    type is (t_field_dynamic_3d)
      varname = p_field%get_name()
      n_dims = p_field%get_dim()
      missing_value = p_field%get_missing_value()
    class default
      call throw_error("fieldset :: read_field_subdomains", "Field "//trim(p_field%get_name())//" is not 2D or 3D")
    end select

#ifdef DEBUG
    DBG, "Reading "//trim(varname)
    debug(field_idx)
    debug(this%u_idx)
    debug(this%v_idx)
#endif

    if (n_dims == 2) then
      allocate (buffer(this%nx, this%ny, 1))
      allocate (count(3))
      count(3) = 1
      start = [1, 1, this%read_idx]
    else if (n_dims == 3) then
      allocate (buffer(this%nx, this%ny, this%nz))
      allocate (count(4))
      count(3) = this%nz
      count(4) = 1
      start = [1, 1, 1, this%read_idx]
    else
      call throw_error("fieldset :: read_field", "Wrong number of dimensions: "//trim(varname))
    end if

    ! Makes the uninitialised ubound warning go away, but might not be necessary actually
    buffer = ZERO

    do i_subdom = 0, this%nproc - 1
      write (subdom_filename, "(a,i0.4,a)") trim(this%current_path)//"/"//trim(this%file_prefix), &
        i_subdom, trim(this%file_suffix)//".nc"

      ioff = this%pmap(i_subdom + 1, 1); joff = this%pmap(i_subdom + 1, 2)
      istart = 1 + ioff; jstart = 1 + joff
      if (ioff .lt. 0) then
        count(1) = this%nxp - abs(ioff)
        istart = 1
      else
        count(1) = this%nxp
      end if
      if (joff .lt. 0) then
        count(2) = this%nyp - abs(joff)
        jstart = 1
      else
        count(2) = this%nyp
      end if

      allocate (buffer_subdom(count(1), count(2), count(3)))

      if (n_dims == 2) then
        call nc_read_real_3d(trim(subdom_filename), trim(varname), start, count, buffer_subdom)
      else
        call nc_read_real_4d(trim(subdom_filename), trim(varname), start, count, buffer_subdom)
      end if
      buffer(istart:istart + count(1) - 1, jstart:jstart + count(2) - 1, 1:count(3)) = buffer_subdom

      deallocate (buffer_subdom)

    end do

    seamask = this%domain%get_seamask()

    select case (c_field)
    case ("v")
      DBG, "Doing case 'v'"
      do j = 1, this%ny
        do i = 1, this%nx
          if (seamask(i, j) == DOM_LAND) then
            buffer(i, j, :) = ZERO
          end if
          ! This is needed due to the A grid
          if (this%u_mask(i, j) > 0) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
      do j = 1, this%ny
        do i = 1, this%nx
          if (this%v_mask(i, j) > 0) then
            buffer(i, j, :) = buffer(i + 1, j, :)
          else if (this%v_mask(i, j) < 0) then
            buffer(i + 1, j, :) = buffer(i, j, :)
          end if
        end do
      end do
    case ("u")
      DBG, "Doing case 'u'"
      do j = 1, this%ny
        do i = 1, this%nx
#ifdef DEBUG
          ! TODO: buffer(i, j+1) is modified after this pass, so it gets
          ! overwritten. Maybe add a "to be modified" flag to the mask, or some
          ! kind of use of the where(u_mask == ...) statement
          ! ! WAIT... I might be wrong...
          if ((i == 1390) .and. (j == 441)) then ! 1390 and 441 on GOF data!!!
            DBG, i, j
            DBG, seamask(i, j)
            DBG, this%v_mask(i, j)
          end if
#endif
          if (seamask(i, j) == DOM_LAND) then
            buffer(i, j, :) = ZERO
          end if
          if (this%v_mask(i, j) > 0) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
      do j = 1, this%ny
        do i = 1, this%nx
#ifdef DEBUG
          if ((i == 1390) .and. (j == 441)) then
            DBG, i, j
            DBG, seamask(i, j)
            DBG, this%u_mask(i, j)
            DBG, buffer(i, j, 1)
          end if
#endif
          if (this%u_mask(i, j) > 0) then
            buffer(i, j, :) = buffer(i, j + 1, :)
          else if (this%u_mask(i, j) < 0) then
            buffer(i, j + 1, :) = buffer(i, j, :)
          end if
        end do
      end do
    case default
      DBG, "Doing case 'default'"
      do j = 1, this%ny
        do i = 1, this%nx
          if (seamask(i, j) == DOM_LAND) then
            buffer(i, j, :) = ZERO
          end if
        end do
      end do
    end select

    where (buffer == missing_value) buffer = ZERO

    if (n_dims == 3 .and. this%zax_style /= STATIC_DEPTH_VALUES) then
      t_nan = all(buffer(:, :, this%zax_top_idx) == ZERO)
      b_nan = all(buffer(:, :, this%zax_bot_idx) == ZERO)
    end if

    select type (p_field)
    type is (t_field_dynamic_2d)
      call p_field%set(buffer(:, :, 1))
    type is (t_field_dynamic_3d)
      call p_field%set(buffer, t_nan, b_nan)
    end select

  end subroutine read_field_subdomains
  !===========================================
  ! INFO FUNCTIONS
  !===========================================
  subroutine list_fields(this)
    class(t_fieldset), intent(in) :: this

    call this%fields%get_info()

  end subroutine list_fields
  !===========================================
  logical function has_field(this, field_name) result(res)
    class(t_fieldset), intent(in) :: this
    character(*), intent(in)      :: field_name

    res = .false.
    if (this%get_field_index(field_name) > 0) then
      res = .true.
    end if

    return
  end function has_field
  !===========================================
  logical function top_is_nan(this, field_name)
    class(t_fieldset), intent(in) :: this
    character(*), intent(in) :: field_name
    class(t_variable), pointer :: p_field

    top_is_nan = .false.
    call this%fields%get_item(trim(field_name), p_field)
    select type (p_field)
    type is (t_field_static_3d)
      top_is_nan = p_field%top_is_nan()
    type is (t_field_dynamic_3d)
      top_is_nan = p_field%top_is_nan()
    class default
      call throw_error("fieldset :: top_is_nan", "Field is not a 3D field")
    end select

    return
  end function top_is_nan
  !===========================================
  logical function bottom_is_nan(this, field_name)
    class(t_fieldset), intent(in) :: this
    character(*), intent(in) :: field_name
    class(t_variable), pointer :: p_field

    bottom_is_nan = .false.
    call this%fields%get_item(trim(field_name), p_field)
    select type (p_field)
    type is (t_field_static_3d)
      bottom_is_nan = p_field%bottom_is_nan()
    type is (t_field_dynamic_3d)
      bottom_is_nan = p_field%bottom_is_nan()
    class default
      call throw_error("fieldset :: bottom_is_nan", "Field is not a 3D field")
    end select

    return
  end function bottom_is_nan
  !===========================================
  integer function zax_direction(this)
    class(t_fieldset), intent(in) :: this
    zax_direction = this%zax_dir
    return
  end function zax_direction
  !===========================================
  real(rk) function sealevel(this, t, i, j) result(res)
    !---------------------------------------------
    ! Should I keep this here?
    !---------------------------------------------
    class(t_fieldset), intent(in) :: this
    real(rk), intent(in)          :: t, i, j
    class(t_variable), pointer    :: p_field
    real(rk)                      :: zax(this%nz)

    res = ZERO
    if (this%fields%key_exists("ELEV")) then
      call this%fields%get_item("ELEV", p_field)
      select type (p_field)
      type is (t_field_dynamic_2d)
        res = p_field%get(t, i, j)
      class default
        call throw_error("fieldset :: sealevel", "ELEV is not a dynamic 2D field. I don't know what to do.")
      end select
    else
      select case (this%zax_style)
      case (DEPTH_VALUES, LAYER_THICKNESS)
        zax = this%get_zax(t, int(i), int(j))
        res = zax(this%nz)
      case (STATIC_DEPTH_VALUES)
        zax = this%get_zax()
        res = zax(this%nz)
      end select
    end if

    return
  end function sealevel

end module mod_fieldset
