#include "cppdefs.h"
module mod_statevector
  !----------------------------------------------------------------
  ! State vector module.
  ! The state vector holds the particle variables to integrate.
  ! The idea for the state vector was inspired by/stolen from the MOHID Lagrangian model.
  !----------------------------------------------------------------
  use mod_common
  use mod_particle
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_statevector
  !---------------------------------------------
  type t_tracer_pointer
    class(t_particle), pointer :: ptr
  end type t_tracer_pointer
  !---------------------------------------------
  type t_statevector
    type(t_tracer_pointer) :: trc        !< Pointer to the particle
    real(rk), allocatable  :: current(:) !< State array
    logical                :: active
    integer                :: state      !< Particle state
    integer                :: nvars
  contains
    procedure :: to_tracer
    procedure :: copy
    procedure :: clean
  end type t_statevector
  !===================================================
contains
  !===========================================
  subroutine to_tracer(this)
    class(t_statevector), intent(in) :: this
    class(t_particle), pointer :: p_trc

    if (associated(this%trc%ptr)) then
      p_trc => this%trc%ptr
      call p_trc%set_state_array(this%current)
    else
      call throw_error("statevector :: to_tracer", "Particle pointer is not associated.")
    end if

    return
  end subroutine to_tracer
  !===========================================
  subroutine copy(this, source)
    !------------------------------------------------
    ! Copy the state vector. (source -> this)
    !------------------------------------------------
    class(t_statevector), intent(inout) :: this
    class(t_statevector), intent(in) :: source

    this%nvars = source%nvars
    this%state = source%state
    this%active = source%active

    if (associated(this%trc%ptr)) then
      nullify (this%trc%ptr)
    end if
    this%trc%ptr => source%trc%ptr

    if (allocated(this%current)) then
      deallocate (this%current)
    end if

    allocate (this%current(this%nvars))
    this%current = source%current

    return
  end subroutine copy
  !===========================================
  subroutine clean(this)
    !------------------------------------------------
    ! Clean up the state vector.
    !------------------------------------------------
    class(t_statevector) :: this

    if (associated(this%trc%ptr)) then
      nullify (this%trc%ptr)
    end if
    if (allocated(this%current)) then
      deallocate (this%current)
    end if
    this%state = 0
    this%active = .false.

    return
  end subroutine clean

end module mod_statevector
