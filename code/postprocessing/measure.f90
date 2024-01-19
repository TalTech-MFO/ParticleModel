#include "cppdefs.h"
module mod_measure
  !----------------------------------------------------------------
  ! [module description]
  !----------------------------------------------------------------
  use mod_errors
  use mod_precdefs
  use mod_variable
  use mod_particle
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_measure
  !---------------------------------------------
  type, abstract, extends(t_variable) :: t_measure
  contains
    procedure(run_func), deferred              :: run
    procedure(get_func), deferred              :: get
    procedure(set_func), deferred              :: set
    procedure(reset_func), deferred            :: reset
    procedure(get_measure_info_func), deferred :: get_measure_type
    procedure(get_measure_info_func), deferred :: get_variable_name
    procedure(get_measure_info_func), deferred :: get_measure_description
  end type t_measure
  !---------------------------------------------
  abstract interface
    subroutine run_func(this, p, i, j, k)
      import :: t_measure, t_particle
      class(t_measure), intent(inout) :: this
      type(t_particle), intent(in) :: p
      integer, intent(in) :: i, j, k
    end subroutine run_func
    function get_func(this) result(res)
      import :: t_measure, rk
      class(t_measure), intent(inout) :: this
      real(rk), dimension(:, :, :), allocatable :: res
    end function get_func
    subroutine reset_func(this)
      import :: t_measure
      class(t_measure), intent(inout) :: this
    end subroutine reset_func
    function get_measure_info_func(this) result(res)
      import :: t_measure, LEN_CHAR_L
      class(t_measure), intent(inout) :: this
      character(len=LEN_CHAR_L) :: res
    end function get_measure_info_func
    subroutine set_func(this, data)
      import :: t_measure, rk
      class(t_measure), intent(inout) :: this
      real(rk), intent(in) :: data(:, :, :)
    end subroutine set_func
  end interface
end module mod_measure
