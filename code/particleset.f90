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
    type(t_particle), pointer        :: particles(:)
    type(t_statevector), allocatable :: sv(:)
  contains
    private
    procedure :: init
    procedure :: release
    procedure :: to_statevector
    procedure :: from_statevector
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
  subroutine to_statevector(this, sv)
    class(t_particleset), intent(inout) :: this
    class(t_statevector), intent(inout) :: sv
    integer :: ipart, nvars

    allocate (this%sv(size(this%particles)), stat=ierr)
    if (ierr /= 0) then
      call throw_error("particleset :: to_statevector", "Could not allocate statevector array")
    end if

    ! TODO: OpenMP parallel loop
    do ipart = 1, size(this%particles)
      nvars = this%particles(ipart)%get_nvars()
      allocate (this%sv(ipart)%current(nvars), stat=ierr)
      if (ierr /= 0) then
        call throw_error("particleset :: to_statevector", "Could not allocate statevector current array")
      end if
      this%sv(ipart)%trc%ptr => this%particles(ipart)
      this%sv(ipart)%current(:) = this%particles(ipart)%get_state_array()
      this%sv(ipart)%active = this%particles(ipart)%st%is_active
      this%sv(ipart)%state = this%particles(ipart)%st%state
      this%sv(ipart)%nvars = nvars
    end do

    return
  end subroutine to_statevector
  !===========================================
  subroutine from_statevector(this)
    class(t_particleset), intent(in) :: this
    integer :: ipart

    do ipart = 1, size(this%sv)
      call this%sv(ipart)%to_tracer()
    end do

    return
  end subroutine from_statevector

end module mod_particleset
