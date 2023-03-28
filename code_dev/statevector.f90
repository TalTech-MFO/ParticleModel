#include "cppdefs.h"
module mod_statevector
  !----------------------------------------------------------------
  ! State vector module.
  ! The state vector holds the particle variables to integrate.
  ! The idea for the state vector was inspired by/stolen from the MOHID Lagrangian model.
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
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
    type(t_tracer_pointer), allocatable :: trc(:) !< Pointer to the particle
    real(rk), allocatable               :: state(:, :) !< State vector
    integer                             :: nparticles !< Number of particles at the current time step
    integer                             :: nvars
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
    integer :: i

    do i = 1, this%nparticles
      if (associated(this%trc(i)%ptr)) then
        p_trc => this%trc(i)%ptr
        call p_trc%set_state_array(this%state(i, :))
      else
        call throw_error("statevector :: to_tracer", "Particle pointer is not associated.")
      end if
    end do

    return
  end subroutine to_tracer
  !===========================================
  subroutine copy(this, source)
    !------------------------------------------------
    ! Copy the state vector. (source -> this)
    !------------------------------------------------
    class(t_statevector), intent(inout) :: this
    class(t_statevector), intent(in) :: source

    this%nparticles = source%nparticles
    this%nvars = source%nvars

    if (allocated(this%state)) then
      deallocate (this%state)
    end if

    allocate (this%state(this%nvars, this%nparticles))
    this%state = source%state

    if (allocated(this%trc)) then
      deallocate (this%trc)
    end if

    allocate (this%trc(this%nparticles))
    this%trc = source%trc

    return
  end subroutine copy
  !===========================================
  subroutine clean(this)
    !------------------------------------------------
    ! Clean up the state vector.
    !------------------------------------------------
    class(t_statevector) :: this

    if (allocated(this%state)) then
      deallocate (this%state)
    end if
    if (allocated(this%trc)) then
      deallocate (this%trc)
    end if

    this%nparticles = 0

    return
  end subroutine clean

end module mod_statevector
