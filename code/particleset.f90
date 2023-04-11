#include "cppdefs.h"
#include "ompdefs.h"
#include "particle.h"
module mod_particleset
  !----------------------------------------------------------------
  ! Particle set module.
  ! This module defines the particleset class which holds all the
  ! particle objects (similarly to the fieldset class).
  ! The particle set class handles releasing new particles and gathering
  ! information to and from the particles.
  !----------------------------------------------------------------
  use mod_common
  use mod_particle
  use mod_statevector
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_particleset
  !---------------------------------------------
  type t_particleset
    private
    type(t_particle), pointer :: particles(:)
    type(t_statevector)       :: sv
  contains
    private
    procedure :: init
    procedure :: release
    procedure, public :: tracer_to_statevector
    procedure, public :: statevector_to_tracer
  end type t_particleset
  !---------------------------------------------
  integer :: ierr
  !===================================================
contains
  !===========================================
  subroutine init(this)
    class(t_particleset), intent(in) :: this
  end subroutine init
  !===========================================
  subroutine release(this)
    class(t_particleset), intent(in) :: this
  end subroutine release
  !===========================================
  subroutine tracer_to_statevector(this)
    class(t_particleset), intent(inout) :: this
    integer :: ipart, nvars

    ! ! Allocate statevector array
    allocate (this%sv%sa(size(this%particles)), stat=ierr)
    if (ierr /= 0) then
      call throw_error("particleset :: tracer_to_statevector", "Could not allocate statevector array")
    end if

    ! TODO: OpenMP
    do ipart = 1, size(this%particles)
      nvars = this%particles(ipart)%get_num_vars()
      allocate (this%sv%sa(ipart)%current(nvars), stat=ierr)
      if (ierr /= 0) then
        call throw_error("particleset :: tracer_to_statevector", "Could not allocate statevector current array")
      end if
      this%sv%sa(ipart)%trc%ptr => this%particles(ipart)
      this%sv%sa(ipart)%current(:) = this%particles(ipart)%get_state_array()
      this%sv%sa(ipart)%active = this%particles(ipart)%st%is_active
      this%sv%sa(ipart)%state = this%particles(ipart)%st%state
      this%sv%sa(ipart)%nvars = nvars
    end do

    return
  end subroutine tracer_to_statevector
  !===========================================
  subroutine statevector_to_tracer(this)
    class(t_particleset), intent(in) :: this

    call this%sv%to_tracer()

    return
  end subroutine statevector_to_tracer

end module mod_particleset
