#include "cppdefs.h"
#include "particle.h"
module mod_snapshot
  !----------------------------------------------------------------
  ! Instantaneous measures (snapshots) of the particle properties.
  ! The grid is reset every time step.
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use mod_measure
  use mod_particle
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_counter_snapshot, t_settled_counter_snapshot, t_property_snapshot, t_settled_property_snapshot
  public :: make_snapshot
  !---------------------------------------------
  type, abstract, extends(t_measure) :: t_basic_snapshot
    integer :: nlon, nlat, ndep
    real(rk), allocatable :: data(:, :, :)
    character(len=LEN_CHAR_L) :: variable_name = ""
    character(len=LEN_CHAR_L) :: measure_type = "snapshot"
    character(len=LEN_CHAR_L) :: measure_description = "instantaneous counts in bins"
  contains
    procedure :: reset => reset_basic_snapshot
    procedure :: get => get_data_basic_snapshot
    procedure :: set => set_data_basic_snapshot
    procedure :: get_measure_type => get_type_basic_snapshot
    procedure :: get_variable_name => get_variable_name_basic_snapshot
    procedure :: get_measure_description => get_description_basic_snapshot
  end type t_basic_snapshot
  !---------------------------------------------
  type, extends(t_basic_snapshot) :: t_counter_snapshot
  contains
    procedure :: run => run_counter_snapshot
    ! procedure :: get_measure_type => get_type_counter_snapshot
  end type t_counter_snapshot
  !---------------------------------------------
  interface t_counter_snapshot
    module procedure :: ctor_counter_snapshot
  end interface t_counter_snapshot
  !---------------------------------------------
  type, extends(t_basic_snapshot) :: t_settled_counter_snapshot
  contains
    procedure :: run => run_settled_counter_snapshot
    ! procedure :: get_measure_type => get_type_settled_counter_snapshot
  end type t_settled_counter_snapshot
  !---------------------------------------------
  interface t_settled_counter_snapshot
    module procedure :: ctor_settled_counter_snapshot
  end interface t_settled_counter_snapshot
  !---------------------------------------------
  type, extends(t_basic_snapshot) :: t_property_snapshot
    ! character(len=LEN_CHAR_L) :: variable_name
  contains
    procedure :: run => run_property_snapshot
    ! procedure :: get_measure_type => get_type_property_snapshot
    ! procedure :: get_variable_name => get_variable_name_property_snapshot
    ! procedure :: get_measure_description => get_description_property_snapshot
  end type t_property_snapshot
  !---------------------------------------------
  interface t_property_snapshot
    module procedure :: ctor_property_snapshot
  end interface t_property_snapshot
  !---------------------------------------------
  type, extends(t_basic_snapshot) :: t_settled_property_snapshot
    ! character(len=LEN_CHAR_L) :: variable_name
  contains
    procedure :: run => run_settled_property_snapshot
  end type t_settled_property_snapshot
  !---------------------------------------------
  interface t_settled_property_snapshot
    module procedure :: ctor_settled_property_snapshot
  end interface t_settled_property_snapshot
  !---------------------------------------------
  interface make_snapshot
    module procedure :: make_counter_snapshot, make_property_snapshot
  end interface make_snapshot
  !===================================================
contains
  !===========================================
  function make_counter_snapshot(name, nlon, nlat, ndep, unit, bottom) result(res)
    !---------------------------------------------
    ! Factory functions
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat, ndep
    logical, intent(in) :: bottom
    class(t_basic_snapshot), allocatable :: res

    if (bottom) then
      res = t_settled_counter_snapshot(name, nlon, nlat, unit)
    else
      res = t_counter_snapshot(name, nlon, nlat, ndep, unit)
    end if

  end function make_counter_snapshot
  !===========================================
  function make_property_snapshot(name, nlon, nlat, ndep, unit, bottom, variable_name) result(res)
    !---------------------------------------------
    ! Factory functions
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat, ndep
    logical, intent(in) :: bottom
    character(len=*), intent(in) :: variable_name
    class(t_basic_snapshot), allocatable :: res

    if (bottom) then
      res = t_settled_property_snapshot(name, nlon, nlat, unit, variable_name)
    else
      res = t_property_snapshot(name, nlon, nlat, ndep, unit, variable_name)
    end if

  end function make_property_snapshot
  !===========================================
  type(t_counter_snapshot) function ctor_counter_snapshot(name, nlon, nlat, ndep, unit) result(this)
    !---------------------------------------------
    ! Creating the snapshot object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat, ndep

    dbghead(postprocess :: ctor_counter_snapshot)
    FMT1, "======== Init counter snapshot ========"
    FMT2, "name: ", name
    FMT2, "unit: ", unit
    FMT2, "nlon: ", nlon
    FMT2, "nlat: ", nlat
    FMT2, "ndep: ", ndep

    this%nlon = nlon
    this%nlat = nlat
    this%ndep = ndep

    ! this%measure_type = "counter snapshot"
    this%measure_description = "instantaneous counts of suspended particles"

    call this%set_name(name)
    call this%set_units(unit)
    call this%set_dim(3)
    allocate (this%data(nlon, nlat, ndep), source=ZERO)

    dbgtail(postprocess :: ctor_counter_snapshot)
  end function ctor_counter_snapshot
  !===========================================
  type(t_settled_counter_snapshot) function ctor_settled_counter_snapshot(name, nlon, nlat, unit) result(this)
    !---------------------------------------------
    ! Creating the snapshot object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat

    dbghead(postprocess :: ctor_settled_counter_snapshot)
    FMT1, "======== Init settled counter snapshot ========"
    FMT2, "name: ", name
    FMT2, "unit: ", unit
    FMT2, "nlon: ", nlon
    FMT2, "nlat: ", nlat
    FMT2, "ndep: 1"

    this%nlon = nlon
    this%nlat = nlat
    this%ndep = 1

    ! this%measure_type = "settled counter snapshot"
    this%measure_description = "instantaneous counts of settled particles"

    call this%set_name(name)
    call this%set_units(unit)
    call this%set_dim(2)
    allocate (this%data(nlon, nlat, 1), source=ZERO)

    dbgtail(postprocess :: ctor_settled_counter_snapshot)
  end function ctor_settled_counter_snapshot
  !===========================================
  type(t_property_snapshot) function ctor_property_snapshot(name, nlon, nlat, ndep, unit, variable_name) result(this)
    !---------------------------------------------
    ! Creating the snapshot object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat, ndep
    character(len=*), intent(in) :: variable_name

    dbghead(postprocess :: ctor_property_snapshot)
    FMT1, "======== Init property snapshot ========"
    FMT2, "name: ", name
    FMT2, "unit: ", unit
    FMT2, "nlon: ", nlon
    FMT2, "nlat: ", nlat
    FMT2, "ndep: ", ndep

    this%nlon = nlon
    this%nlat = nlat
    this%ndep = ndep
    this%variable_name = variable_name

    ! this%measure_type = "property snapshot"
    this%measure_description = "instantaneous "//trim(this%variable_name)//" of suspended particles"

    call this%set_name(name)
    call this%set_units(unit)
    call this%set_dim(3)
    allocate (this%data(nlon, nlat, ndep), source=ZERO)

    dbgtail(postprocess :: ctor_property_snapshot)
  end function ctor_property_snapshot
  !===========================================
  type(t_settled_property_snapshot) function ctor_settled_property_snapshot(name, nlon, nlat, unit, variable_name) result(this)
    !---------------------------------------------
    ! Creating the snapshot object
    !---------------------------------------------
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: unit
    integer, intent(in) :: nlon, nlat
    character(len=*), intent(in) :: variable_name

    dbghead(postprocess :: ctor_settled_property_snapshot)
    FMT1, "======== Init settled property snapshot ========"
    FMT2, "name: ", name
    FMT2, "unit: ", unit
    FMT2, "nlon: ", nlon
    FMT2, "nlat: ", nlat
    FMT2, "ndep: 1"

    this%nlon = nlon
    this%nlat = nlat
    this%ndep = 1
    this%variable_name = variable_name

    ! this%measure_type = "settled property snapshot"
    this%measure_description = "instantaneous "//trim(this%variable_name)//" of settled particles"

    call this%set_name(name)
    call this%set_units(unit)
    call this%set_dim(2)
    allocate (this%data(nlon, nlat, 1), source=ZERO)

    dbgtail(postprocess :: ctor_settled_property_snapshot)
  end function ctor_settled_property_snapshot
  !===========================================
  subroutine run_counter_snapshot(this, p, i, j, k)
    class(t_counter_snapshot), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j, k

    if (p%state /= ST_SUSPENDED) return

    !$omp atomic update
    this%data(i, j, k) = this%data(i, j, k) + ONE

  end subroutine run_counter_snapshot
  !===========================================
  subroutine run_settled_counter_snapshot(this, p, i, j, k)
    class(t_settled_counter_snapshot), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j, k

    if (p%state /= ST_BOTTOM) return

    !$omp atomic update
    this%data(i, j, 1) = this%data(i, j, 1) + ONE

  end subroutine run_settled_counter_snapshot
  !===========================================
  subroutine run_property_snapshot(this, p, i, j, k)
    class(t_property_snapshot), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j, k
    real(rk) :: value

    if (p%state /= ST_SUSPENDED) return

    value = p%get_variable(this%variable_name)

    !$omp atomic update
    this%data(i, j, k) = this%data(i, j, k) + value

  end subroutine run_property_snapshot
  !===========================================
  subroutine run_settled_property_snapshot(this, p, i, j, k)
    class(t_settled_property_snapshot), intent(inout) :: this
    type(t_particle), intent(in) :: p
    integer, intent(in) :: i, j, k
    real(rk) :: value

    if (p%state /= ST_BOTTOM) return

    value = p%get_variable(this%variable_name)

    !$omp atomic update
    this%data(i, j, 1) = this%data(i, j, 1) + value

  end subroutine run_settled_property_snapshot
  !===========================================
  function get_data_basic_snapshot(this) result(res)
    class(t_basic_snapshot), intent(inout) :: this
    real(rk), dimension(:, :, :), allocatable :: res

    res = this%data

  end function get_data_basic_snapshot
  !===========================================
  subroutine set_data_basic_snapshot(this, data)
    class(t_basic_snapshot), intent(inout) :: this
    real(rk), dimension(:, :, :), intent(in) :: data
    character(len=LEN_CHAR_L) :: msg

    if (any(shape(this%data) /= shape(data))) then
      write (msg, "('data shape mismatch: ', 3i0, ' /= ', 3i0)") shape(this%data), shape(data)
      call throw_error("basic_snapshot :: set_data_basic_snapshot", trim(msg))
    end if

    this%data = data

  end subroutine set_data_basic_snapshot
  !===========================================
  subroutine reset_basic_snapshot(this)
    class(t_basic_snapshot), intent(inout) :: this

    this%data = ZERO

  end subroutine reset_basic_snapshot
  !===========================================
  function get_type_basic_snapshot(this) result(res)
    class(t_basic_snapshot), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    ! res = "snapshot"
    res = this%measure_type

  end function get_type_basic_snapshot
  ! !===========================================
  ! function get_type_counter_snapshot(this) result(res)
  !   class(t_counter_snapshot), intent(inout) :: this
  !   character(len=LEN_CHAR_L) :: res

  !   res = "counter snapshot"

  ! end function get_type_counter_snapshot
  ! !===========================================
  ! function get_type_settled_counter_snapshot(this) result(res)
  !   class(t_settled_counter_snapshot), intent(inout) :: this
  !   character(len=LEN_CHAR_L) :: res

  !   res = "settled counter snapshot"

  ! end function get_type_settled_counter_snapshot
  ! !===========================================
  ! function get_type_property_snapshot(this) result(res)
  !   class(t_property_snapshot), intent(inout) :: this
  !   character(len=LEN_CHAR_L) :: res

  !   res = "property snapshot"

  ! end function get_type_property_snapshot
  !===========================================
  function get_variable_name_basic_snapshot(this) result(res)
    class(t_basic_snapshot), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    ! res = ""
    res = this%variable_name

  end function get_variable_name_basic_snapshot
  !===========================================
  ! function get_variable_name_property_snapshot(this) result(res)
  !   class(t_property_snapshot), intent(inout) :: this
  !   character(len=LEN_CHAR_L) :: res

  !   res = this%variable_name

  ! end function get_variable_name_property_snapshot
  !===========================================
  function get_description_basic_snapshot(this) result(res)
    class(t_basic_snapshot), intent(inout) :: this
    character(len=LEN_CHAR_L) :: res

    ! res = "instantaneous counts in bins"
    res = this%measure_description

  end function get_description_basic_snapshot
  ! !===========================================
  ! function get_description_property_snapshot(this) result(res)
  !   class(t_property_snapshot), intent(inout) :: this
  !   character(len=LEN_CHAR_L) :: res

  !   res = "instantaneous "//trim(this%variable_name)//" in bins"

  ! end function get_description_property_snapshot

end module mod_snapshot
