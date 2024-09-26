#include "cppdefs.h"
#include "particle.h"
module mod_accumulator
  !----------------------------------------------------------------
  ! Accumulated properties (i.e. the grid is NOT reset after each time step).
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use mod_measure
  use mod_particle
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_counter_accumulator, t_settled_counter_accumulator, t_property_accumulator, t_settled_property_accumulator
  public :: make_accumulator
  !---------------------------------------------
  type, abstract, extends(t_measure) :: t_basic_accumulator
    integer :: nlon, nlat, ndep
    real(rk), allocatable :: data(:, :, :)
    character(len=LEN_CHAR_L) :: variable_name = ""
    character(len=LEN_CHAR_L) :: measure_type = "accumulator"
    character(len=LEN_CHAR_L) :: measure_description = "accumulated counts in bins"
  contains
    procedure :: get => get_data_basic_accumulator
    procedure :: set => set_data_basic_accumulator
    procedure :: reset => reset_basic_accumulator
    procedure :: get_measure_type => get_type_basic_accumulator
    procedure :: get_variable_name => get_variable_name_basic_accumulator
    procedure :: get_measure_description => get_description_basic_accumulator
  end type t_basic_accumulator
  !---------------------------------------------
  type, extends(t_basic_accumulator) :: t_counter_accumulator
  contains
    procedure :: run => run_counter_accumulator
  end type t_counter_accumulator
  !---------------------------------------------
  interface t_counter_accumulator
    module procedure :: ctor_counter_accumulator
  end interface t_counter_accumulator
  !---------------------------------------------
  type, extends(t_basic_accumulator) :: t_settled_counter_accumulator
  contains
    procedure :: run => run_settled_counter_accumulator
  end type t_settled_counter_accumulator
  !---------------------------------------------
  interface t_settled_counter_accumulator
    module procedure :: ctor_settled_counter_accumulator
  end interface t_settled_counter_accumulator
  !---------------------------------------------
  type, extends(t_basic_accumulator) :: t_property_accumulator
    ! character(len=LEN_CHAR_L) :: variable_name
  contains
    procedure :: run => run_property_accumulator
    ! procedure :: get_variable_name => get_variable_name_property_accumulator
    ! procedure :: get_measure_description => get_description_property_accumulator
  end type t_property_accumulator
  !---------------------------------------------
  interface t_property_accumulator
    module procedure :: ctor_property_accumulator
  end interface t_property_accumulator
  !---------------------------------------------
  type, extends(t_basic_accumulator) :: t_settled_property_accumulator
    ! character(len=LEN_CHAR_L) :: variable_name
  contains
    procedure :: run => run_settled_property_accumulator
    ! procedure :: get_variable_name => get_variable_name_property_accumulator
    ! procedure :: get_measure_description => get_description_property_accumulator
  end type t_settled_property_accumulator
  !---------------------------------------------
  interface t_settled_property_accumulator
    module procedure :: ctor_settled_property_accumulator
  end interface t_settled_property_accumulator
  !---------------------------------------------
  interface make_accumulator
    module procedure :: make_counter_accumulator, make_property_accumulator
  end interface make_accumulator
  !===================================================
contains
  !===========================================
  function make_counter_accumulator(name, nlon, nlat, ndep, unit, bottom) result(res)
    !---------------------------------------------
    ! Factory functions
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat, ndep
    logical, intent(in) :: bottom
    class(t_basic_accumulator), allocatable :: res

    if (bottom) then
      res = t_settled_counter_accumulator(name, nlon, nlat, unit)
    else
      res = t_counter_accumulator(name, nlon, nlat, ndep, unit)
    end if
  end function make_counter_accumulator
  !===========================================
  function make_property_accumulator(name, nlon, nlat, ndep, unit, bottom, variable_name) result(res)
    !---------------------------------------------
    ! Factory functions
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat, ndep
    character(len=*), intent(in) :: variable_name
    logical, intent(in) :: bottom
    class(t_basic_accumulator), allocatable :: res

    if (bottom) then
      res = t_settled_property_accumulator(name, nlon, nlat, unit, variable_name)
    else
      res = t_property_accumulator(name, nlon, nlat, ndep, unit, variable_name)
    end if
  end function make_property_accumulator
  !===========================================
  type(t_counter_accumulator) function ctor_counter_accumulator(name, nlon, nlat, ndep, unit) result(this)
    !---------------------------------------------
    ! Creating the accumulator object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat, ndep

    dbghead(postprocess :: ctor_counter_accumulator)
    FMT1, "======== Init counter accumulator ========"
    FMT2, "name: ", name
    FMT2, "unit: ", unit
    FMT2, "nlon: ", nlon
    FMT2, "nlat: ", nlat
    FMT2, "ndep: ", ndep

    this%nlon = nlon
    this%nlat = nlat
    this%ndep = ndep

    ! this%measure_type = "counter accumulator"
    this%measure_description = "accumulated counts of suspended particles"

    call this%set_name(name)
    call this%set_units(unit)
    call this%set_dim(3)
    allocate (this%data(nlon, nlat, ndep), source=ZERO)

    dbgtail(postprocess :: ctor_counter_accumulator)
  end function ctor_counter_accumulator
  !===========================================
  type(t_settled_counter_accumulator) function ctor_settled_counter_accumulator(name, nlon, nlat, unit) result(this)
    !---------------------------------------------
    ! Creating the accumulator object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat

    dbghead(postprocess :: ctor_settled_counter_accumulator)
    FMT1, "======== Init settled counter accumulator ========"
    FMT2, "name: ", name
    FMT2, "unit: ", unit
    FMT2, "nlon: ", nlon
    FMT2, "nlat: ", nlat
    FMT2, "ndep: 1"

    this%nlon = nlon
    this%nlat = nlat
    this%ndep = 1

    ! this%measure_type = "settled counter accumulator"
    this%measure_description = "accumulated counts of settled particles"

    call this%set_name(name)
    call this%set_units(unit)
    call this%set_dim(2)
    allocate (this%data(nlon, nlat, 1), source=ZERO)

    dbgtail(postprocess :: ctor_settled_counter_accumulator)
  end function ctor_settled_counter_accumulator
  !===========================================
  type(t_property_accumulator) function ctor_property_accumulator(name, nlon, nlat, ndep, unit, variable_name) result(this)
    !---------------------------------------------
    ! Creating the accumulator object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat, ndep
    character(len=*), intent(in) :: variable_name

    dbghead(postprocess :: ctor_property_accumulator)
    FMT1, "======== Init property accumulator ========"
    FMT2, "name: ", name
    FMT2, "unit: ", unit
    FMT2, "nlon: ", nlon
    FMT2, "nlat: ", nlat
    FMT2, "ndep: ", ndep

    this%nlon = nlon
    this%nlat = nlat
    this%ndep = ndep
    this%variable_name = variable_name

    ! this%measure_type = "property accumulator"
    this%measure_description = "accumulated "//trim(this%variable_name)//" of suspended particles"

    call this%set_name(name)
    call this%set_units(unit)
    call this%set_dim(3)
    allocate (this%data(nlon, nlat, ndep), source=ZERO)

    dbgtail(postprocess :: ctor_property_accumulator)
  end function ctor_property_accumulator
  !===========================================
 type(t_settled_property_accumulator) function ctor_settled_property_accumulator(name, nlon, nlat, unit, variable_name) result(this)
    !---------------------------------------------
    ! Creating the accumulator object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat
    character(len=*), intent(in) :: variable_name

    dbghead(postprocess :: ctor_settled_property_accumulator)
    FMT1, "======== Init settled property accumulator ========"
    FMT2, "name: ", name
    FMT2, "unit: ", unit
    FMT2, "nlon: ", nlon
    FMT2, "nlat: ", nlat
    FMT2, "ndep: 1"

    this%nlon = nlon
    this%nlat = nlat
    this%ndep = 1
    this%variable_name = variable_name

    ! this%measure_type = "settled property accumulator"
    this%measure_description = "accumulated "//trim(this%variable_name)//" of settled particles"

    call this%set_name(name)
    call this%set_units(unit)
    call this%set_dim(2)
    allocate (this%data(nlon, nlat, 1), source=ZERO)

    dbgtail(postprocess :: ctor_settled_property_accumulator)
  end function ctor_settled_property_accumulator
  !===========================================
  subroutine run_counter_accumulator(this, p, i, j, k)
    class(t_counter_accumulator), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j, k

    if (p%state /= ST_SUSPENDED) return

    !$omp atomic update
    this%data(i, j, k) = this%data(i, j, k) + ONE

  end subroutine run_counter_accumulator
  !===========================================
  subroutine run_settled_counter_accumulator(this, p, i, j, k)
    class(t_settled_counter_accumulator), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j, k

    if (p%state /= ST_BOTTOM) return

    !$omp atomic update
    this%data(i, j, 1) = this%data(i, j, 1) + ONE

  end subroutine run_settled_counter_accumulator
  !===========================================
  subroutine run_property_accumulator(this, p, i, j, k)
    class(t_property_accumulator), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j, k
    real(rk) :: value

    dbghead(postprocessor :: run_property_accumulator)

    debug(i); debug(j)

    if (p%state /= ST_SUSPENDED) return

    value = p%get_variable(this%variable_name)

    debug(this%variable_name)
    debug(value)

    !$omp atomic update
    this%data(i, j, k) = this%data(i, j, k) + value

    dbgtail(postprocessor :: run_property_accumulator)
  end subroutine run_property_accumulator
  !===========================================
  subroutine run_settled_property_accumulator(this, p, i, j, k)
    class(t_settled_property_accumulator), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j, k
    real(rk) :: value

    dbghead(postprocessor :: run_settled_property_accumulator)

    debug(i); debug(j)

    if (p%state /= ST_BOTTOM) return

    value = p%get_variable(this%variable_name)

    debug(this%variable_name)
    debug(value)

    !$omp atomic update
    this%data(i, j, 1) = this%data(i, j, 1) + value

    dbgtail(postprocessor :: run_settled_property_accumulator)
  end subroutine run_settled_property_accumulator
  !===========================================
  function get_data_basic_accumulator(this) result(res)
    class(t_basic_accumulator), intent(inout) :: this
    real(rk), dimension(:, :, :), allocatable :: res

    res = this%data

  end function get_data_basic_accumulator
  !===========================================
  subroutine set_data_basic_accumulator(this, data)
    class(t_basic_accumulator), intent(inout) :: this
    real(rk), dimension(:, :, :), intent(in) :: data
    character(len=LEN_CHAR_L) :: msg

    if (any(shape(this%data) /= shape(data))) then
      write (msg, "('data shape mismatch: ', 2i0, ' /= ', 2i0)") shape(this%data), shape(data)
      call throw_error("basic_accumulator :: set_basic_accumulator_data", trim(msg))
    end if
    this%data = data

  end subroutine set_data_basic_accumulator
  !===========================================
  subroutine reset_basic_accumulator(this)
    class(t_basic_accumulator), intent(inout) :: this

    ! The accumulator does not reset

  end subroutine reset_basic_accumulator
  !===========================================
  function get_type_basic_accumulator(this) result(res)
    class(t_basic_accumulator), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    ! res = "accumulator"
    res = this%measure_type

  end function get_type_basic_accumulator
  !===========================================
  function get_variable_name_basic_accumulator(this) result(res)
    class(t_basic_accumulator), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    ! res = ""
    res = this%variable_name

  end function get_variable_name_basic_accumulator
  ! !===========================================
  ! function get_variable_name_property_accumulator(this) result(res)
  !   class(t_property_accumulator), intent(inout) :: this
  !   character(len=LEN_CHAR_L) :: res

  !   res = this%variable_name

  ! end function get_variable_name_property_accumulator
  !===========================================
  function get_description_basic_accumulator(this) result(res)
    class(t_basic_accumulator), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    ! res = "accumulated counts in bins"
    res = this%measure_description

  end function get_description_basic_accumulator
  ! !===========================================
  ! function get_description_property_accumulator(this) result(res)
  !   class(t_property_accumulator), intent(inout) :: this
  !   character(len=LEN_CHAR_L) :: res

  !   res = "accumulated "//trim(this%variable_name)//" in bins"

  ! end function get_description_property_accumulator

end module mod_accumulator
