#include "cppdefs.h"
module mod_accumulator
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
  public :: t_basic_accumulator, t_property_accumulator
  !---------------------------------------------
  type, extends(t_measure) :: t_basic_accumulator
    real(rk), allocatable :: data(:, :)
  contains
    procedure :: run => run_basic_accumulator
    procedure :: get => get_basic_accumulator_data
    procedure :: set => set_basic_accumulator_data
    procedure :: reset => reset_basic_accumulator
    procedure :: get_measure_type => get_basic_accumulator_type
    procedure :: get_variable_name => get_basic_accumulator_variable_name
    procedure :: get_measure_description => get_basic_accumulator_description
  end type t_basic_accumulator
  !---------------------------------------------
  interface t_basic_accumulator
    module procedure :: ctor_basic_accumulator
  end interface t_basic_accumulator
  !---------------------------------------------
  type, extends(t_basic_accumulator) :: t_property_accumulator
    character(len=LEN_CHAR_L) :: variable_name
  contains
    procedure :: run => run_property_accumulator
    procedure :: get_variable_name => get_property_accumulator_variable_name
    procedure :: get_measure_description => get_property_accumulator_description
  end type t_property_accumulator
  !---------------------------------------------
  interface t_property_accumulator
    module procedure :: ctor_property_accumulator
  end interface t_property_accumulator
  !===================================================
contains
  !===========================================
  type(t_basic_accumulator) function ctor_basic_accumulator(name, nlon, nlat, unit) result(this)
    !---------------------------------------------
    ! Creating the accumulator object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat

    dbghead(postprocess :: ctor_basic_accumulator)
    FMT1, "======== Init basic accumulator ========"
    FMT2, "name: ", name
    FMT2, "unit: ", unit
    FMT2, "nlon: ", nlon
    FMT2, "nlat: ", nlat

    call this%set_name(name)
    call this%set_units(unit)
    allocate (this%data(nlon, nlat), source=ZERO)

    dbgtail(postprocess :: ctor_basic_accumulator)
  end function ctor_basic_accumulator
  !===========================================
  type(t_property_accumulator) function ctor_property_accumulator(name, nlon, nlat, unit, variable_name) result(this)
    !---------------------------------------------
    ! Creating the accumulator object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat
    character(len=*), intent(in) :: variable_name

    dbghead(postprocess :: ctor_property_accumulator)
    FMT1, "======== Init property accumulator ========"

    this%t_basic_accumulator = t_basic_accumulator(name, nlon, nlat, unit)
    this%variable_name = variable_name

    dbgtail(postprocess :: ctor_property_accumulator)
  end function ctor_property_accumulator
  !===========================================
  subroutine run_basic_accumulator(this, p, i, j)
    class(t_basic_accumulator), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j

    !$omp atomic update
    this%data(i, j) = this%data(i, j) + ONE

  end subroutine run_basic_accumulator
  !===========================================
  subroutine run_property_accumulator(this, p, i, j)
    class(t_property_accumulator), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j
    real(rk) :: value

    dbghead(postprocessor :: run_property_accumulator)

    debug(i); debug(j)

    value = p%get_variable(this%variable_name)

    debug(this%variable_name)
    debug(value)

    !$omp atomic update
    this%data(i, j) = this%data(i, j) + value

    dbgtail(postprocessor :: run_property_accumulator)
  end subroutine run_property_accumulator
  !===========================================
  function get_basic_accumulator_data(this) result(res)
    class(t_basic_accumulator), intent(inout) :: this
    real(rk), dimension(:, :), allocatable :: res

    res = this%data

  end function get_basic_accumulator_data
  !===========================================
  subroutine set_basic_accumulator_data(this, data)
    class(t_basic_accumulator), intent(inout) :: this
    real(rk), dimension(:, :), intent(in) :: data
    character(len=LEN_CHAR_L) :: msg

    if (any(shape(this%data) /= shape(data))) then
      write (msg, "('data shape mismatch: ', 2i0, ' /= ', 2i0)") shape(this%data), shape(data)
      call throw_error("basic_accumulator :: set_basic_accumulator_data", trim(msg))
    end if
    this%data = data

  end subroutine set_basic_accumulator_data
  !===========================================
  subroutine reset_basic_accumulator(this)
    class(t_basic_accumulator), intent(inout) :: this

    ! The accumulator does not reset

  end subroutine reset_basic_accumulator
  !===========================================
  function get_basic_accumulator_type(this) result(res)
    class(t_basic_accumulator), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    res = "accumulator"

  end function get_basic_accumulator_type
  !===========================================
  function get_basic_accumulator_variable_name(this) result(res)
    class(t_basic_accumulator), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    res = ""

  end function get_basic_accumulator_variable_name
  !===========================================
  function get_property_accumulator_variable_name(this) result(res)
    class(t_property_accumulator), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    res = this%variable_name

  end function get_property_accumulator_variable_name
  !===========================================
  function get_basic_accumulator_description(this) result(res)
    class(t_basic_accumulator), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    res = "accumulated counts in bins"

  end function get_basic_accumulator_description
  !===========================================
  function get_property_accumulator_description(this) result(res)
    class(t_property_accumulator), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    res = "accumulated "//trim(this%variable_name)//" in bins"

  end function get_property_accumulator_description

end module mod_accumulator
