#include "cppdefs.h"
module utils
  use mod_precdefs
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_timer, insert_before_extension
  !---------------------------------------------
  type t_timer
    private
    real(rk) :: m_start = ZERO, m_stop = ZERO
    real(rk) :: duration = ZERO
    real(rk) :: m_max = ZERO
    real(rk) :: m_min = ZERO
    real(rk) :: n = ZERO
    logical :: accumulate_duration = .false.
    character(len=LEN_CHAR_S) :: name
  contains
    procedure :: start
    procedure :: stop
    procedure :: show
    procedure :: set_accumulate
    procedure :: reset
  end type t_timer
  !---------------------------------------------
  !===================================================
contains
  !===========================================
  subroutine set_accumulate(this)
    class(t_timer), intent(inout) :: this

    this%accumulate_duration = .true.

  end subroutine set_accumulate
  !===========================================
  subroutine reset(this)
    class(t_timer), intent(inout) :: this

    this%n = ZERO
    this%duration = ZERO
    this%m_start = ZERO
    this%m_stop = ZERO
    this%name = ""

  end subroutine reset
  !===========================================
  subroutine start(this, name)
    class(t_timer), intent(inout) :: this
    character(len=*), optional, intent(in) :: name
    real(rk) :: current_time

    if (present(name)) this%name = name

    call cpu_time(current_time)
    this%m_start = current_time

  end subroutine start
  !===========================================
  subroutine stop(this)
    class(t_timer), intent(inout) :: this
    real(rk) :: current_time
    real(rk) :: current_duration

    call cpu_time(current_time)
    this%m_stop = current_time
    if (this%accumulate_duration) then
      current_duration = (this%m_stop - this%m_start)
      this%duration = this%duration + current_duration
      this%n = this%n + ONE
      if (this%m_max < current_duration) then
        this%m_max = current_duration
      end if
      if (this%m_min == ZERO) this%m_min = current_duration
      if (this%m_min > current_duration) then
        this%m_min = current_duration
      end if
    end if

  end subroutine stop
  !===========================================
  subroutine show(this, average)
    class(t_timer), intent(in) :: this
    logical, intent(in), optional :: average
    logical :: show_avg

    show_avg = .false.

    if (present(average)) show_avg = average

    if (show_avg) then
 FMT2, "TIMER: ", trim(this%name), ": Average: ", this%duration / this%n, " [sec], max: ", this%m_max, " [sec] min: ", this%m_min, " [sec]"
    else
      FMT2, "TIMER: ", trim(this%name), ":", this%m_stop - this%m_start, " [sec] "
    end if

  end subroutine show
  !===========================================
  function insert_before_extension(filename, substr) result(new_filename)
    character(len=*), intent(in) :: filename, substr
    character(len=LEN_CHAR_L) :: new_filename
    integer :: pos

    ! Find the position of the last dot in the filename
    pos = find_last_dot(filename)

    ! Create the new filename
    if (pos > 0) then
      new_filename = filename(:pos - 1)//substr//filename(pos:)
    else
      new_filename = filename//substr
    end if
  end function insert_before_extension
  !===========================================
  function find_last_dot(str) result(pos)
    character(len=*), intent(in) :: str
    integer :: pos, i

    pos = -1
    do i = len_trim(str), 1, -1
      if (str(i:i) == '.') then
        pos = i
        exit
      end if
    end do
  end function find_last_dot
end module utils
