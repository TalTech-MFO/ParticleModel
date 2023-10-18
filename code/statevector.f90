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
  public :: t_statevector, t_statearray
  !---------------------------------------------
  ! Container to hold a pointer to a particle.
  type t_tracer_pointer
    class(t_particle), pointer :: ptr
  end type t_tracer_pointer
  !---------------------------------------------
  ! Basic container to hold the current state of the particle. This is passed to the kernel.
  type t_statearray
    type(t_tracer_pointer) :: trc        !< Pointer to the particle
    real(rk), allocatable  :: current(:) !< State array
    logical                :: active
    integer                :: state      !< Particle state
    integer                :: nvars
  contains
    procedure :: to_tracer => to_tracer_statearray
    procedure :: copy => copy_statearray
    procedure :: clean => clean_statearray
  end type t_statearray
  !---------------------------------------------
  ! State vector to hold the state of all particles.
  type t_statevector
    type(t_statearray), allocatable :: sa(:)
  contains
    procedure :: get_size => get_size_statevector
    procedure :: to_tracer => to_tracer_statevector
    procedure :: copy => copy_statevector
    procedure :: clean => clean_statevector
  end type t_statevector
  !===================================================
contains
  !===========================================
  pure function get_size_statevector(this)
    !-------------------------------------------
    ! Get the size of the state vector.
    !-------------------------------------------
    class(t_statevector), intent(in) :: this
    integer :: get_size_statevector

    get_size_statevector = size(this%sa)

    return
  end function get_size_statevector
  !===========================================
  subroutine to_tracer_statearray(this)
    !-------------------------------------------
    ! Set the state array to the tracer.
    !-------------------------------------------
    class(t_statearray), intent(in) :: this
    class(t_particle), pointer :: p_trc

    if (associated(this%trc%ptr)) then
      p_trc => this%trc%ptr
      call p_trc%set_state_array(this%current)
    else
      call throw_error("statearray :: to_tracer", "Particle pointer is not associated.")
    end if

    return
  end subroutine to_tracer_statearray
  !===========================================
  subroutine to_tracer_statevector(this)
    !-------------------------------------------
    ! Set all state arrays to the tracer.
    !-------------------------------------------
    class(t_statevector), intent(in) :: this
    integer :: i

    do i = 1, this%get_size()
      call this%sa(i)%to_tracer()
    end do

    return
  end subroutine to_tracer_statevector
  !===========================================
  subroutine copy_statearray(this, source)
    !------------------------------------------------
    ! Copy the state vector. (source -> this)
    !------------------------------------------------
    class(t_statearray), intent(inout) :: this
    class(t_statearray), intent(in) :: source

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
  end subroutine copy_statearray
  !===========================================
  subroutine copy_statevector(this, source)
    !------------------------------------------------
    ! Copy the state vector. (source -> this)
    !------------------------------------------------
    class(t_statevector), intent(inout) :: this
    class(t_statevector), intent(in) :: source
    integer :: i

    if (allocated(this%sa)) then
      deallocate (this%sa)
    end if

    allocate (this%sa(source%get_size()))
    do i = 1, source%get_size()
      call this%sa(i)%copy(source%sa(i))
    end do

    return
  end subroutine copy_statevector
  !===========================================
  subroutine clean_statearray(this)
    !------------------------------------------------
    ! Clean up the state vector.
    !------------------------------------------------
    class(t_statearray) :: this

    if (associated(this%trc%ptr)) then
      nullify (this%trc%ptr)
    end if
    if (allocated(this%current)) then
      deallocate (this%current)
    end if
    this%state = 0
    this%active = .false.

    return
  end subroutine clean_statearray
  !===========================================
  subroutine clean_statevector(this)
    !------------------------------------------------
    ! Clean up the state vector.
    !------------------------------------------------
    class(t_statevector) :: this
    integer :: i

    if (allocated(this%sa)) then
      do i = 1, this%get_size()
        call this%sa(i)%clean()
      end do
      deallocate (this%sa)
    end if

    return
  end subroutine clean_statevector

end module mod_statevector
