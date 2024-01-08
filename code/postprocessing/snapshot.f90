#include "cppdefs.h"
module mod_snapshot
  !----------------------------------------------------------------
  ! [module description]
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use mod_measure
  use mod_particle
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_basic_snapshot, t_property_snapshot
  !---------------------------------------------
  type, extends(t_measure) :: t_basic_snapshot
    real(rk), allocatable :: data(:, :)
  contains
    procedure :: run => run_basic_snapshot
    procedure :: get => get_basic_snapshot_data
    procedure :: set => set_basic_snapshot_data
    procedure :: reset => reset_basic_snapshot
    procedure :: get_measure_type => get_basic_snapshot_type
    procedure :: get_variable_name => get_basic_snapshot_variable_name
    procedure :: get_measure_description => get_basic_snapshot_description
  end type t_basic_snapshot
  !---------------------------------------------
  interface t_basic_snapshot
    module procedure :: ctor_basic_snapshot
  end interface t_basic_snapshot
  !---------------------------------------------
  type, extends(t_basic_snapshot) :: t_property_snapshot
    character(len=LEN_CHAR_L) :: variable_name
  contains
    procedure :: run => run_property_snapshot
    procedure :: get_variable_name => get_property_snapshot_variable_name
    procedure :: get_measure_description => get_property_snapshot_description
  end type t_property_snapshot
  !---------------------------------------------
  interface t_property_snapshot
    module procedure :: ctor_property_snapshot
  end interface t_property_snapshot
  !===================================================
contains
  !===========================================
  type(t_basic_snapshot) function ctor_basic_snapshot(name, nlon, nlat, unit) result(this)
    !---------------------------------------------
    ! Creating the snapshot object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat

    dbghead(postprocess :: ctor_basic_snapshot)
    FMT1, "======== Init basic snapshot ========"
    FMT2, "name: ", name
    FMT2, "unit: ", unit
    FMT2, "nlon: ", nlon
    FMT2, "nlat: ", nlat

    call this%set_name(name)
    call this%set_units(unit)
    allocate (this%data(nlon, nlat), source=ZERO)

    dbgtail(postprocess :: ctor_basic_snapshot)
  end function ctor_basic_snapshot
  !===========================================
  type(t_property_snapshot) function ctor_property_snapshot(name, nlon, nlat, unit, variable_name) result(this)
    !---------------------------------------------
    ! Creating the snapshot object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat
    character(len=*), intent(in) :: variable_name

    dbghead(postprocess :: ctor_property_snapshot)
    FMT1, "======== Init property snapshot ========"

    this%t_basic_snapshot = t_basic_snapshot(name, nlon, nlat, unit)
    this%variable_name = variable_name

    dbgtail(postprocess :: ctor_property_snapshot)
  end function ctor_property_snapshot
  !===========================================
  subroutine run_basic_snapshot(this, p, i, j)
    class(t_basic_snapshot), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j

    !$omp atomic update
    this%data(i, j) = this%data(i, j) + ONE

  end subroutine run_basic_snapshot
  !===========================================
  subroutine run_property_snapshot(this, p, i, j)
    class(t_property_snapshot), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j
    real(rk) :: value

    value = p%get_variable(this%variable_name)

    !$omp atomic update
    this%data(i, j) = this%data(i, j) + value

  end subroutine run_property_snapshot
  !===========================================
  function get_basic_snapshot_data(this) result(res)
    class(t_basic_snapshot), intent(inout) :: this
    real(rk), dimension(:, :), allocatable :: res

    res = this%data

  end function get_basic_snapshot_data
  !===========================================
  subroutine set_basic_snapshot_data(this, data)
    class(t_basic_snapshot), intent(inout) :: this
    real(rk), dimension(:, :), intent(in) :: data
    character(len=LEN_CHAR_L) :: msg

    if (any(shape(this%data) /= shape(data))) then
      write (msg, "('data shape mismatch: ', 2i0, ' /= ', 2i0)") shape(this%data), shape(data)
      call throw_error("basic_snapshot :: set_basic_snapshot_data", trim(msg))
    end if

    this%data = data

  end subroutine set_basic_snapshot_data
  !===========================================
  subroutine reset_basic_snapshot(this)
    class(t_basic_snapshot), intent(inout) :: this

    this%data = ZERO

  end subroutine reset_basic_snapshot
  !===========================================
  function get_basic_snapshot_type(this) result(res)
    class(t_basic_snapshot), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    res = "snapshot"

  end function get_basic_snapshot_type
  !===========================================
  function get_basic_snapshot_variable_name(this) result(res)
    class(t_basic_snapshot), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    res = ""

  end function get_basic_snapshot_variable_name
  !===========================================
  function get_property_snapshot_variable_name(this) result(res)
    class(t_property_snapshot), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    res = this%variable_name

  end function get_property_snapshot_variable_name
  !===========================================
  function get_basic_snapshot_description(this) result(res)
    class(t_basic_snapshot), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    res = "instantaneous counts in bins"

  end function get_basic_snapshot_description
  !===========================================
  function get_property_snapshot_description(this) result(res)
    class(t_property_snapshot), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    res = "instantaneous "//trim(this%variable_name)//" in bins"

  end function get_property_snapshot_description

end module mod_snapshot
