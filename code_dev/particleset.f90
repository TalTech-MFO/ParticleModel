#include "cppdefs.h"
#include "ompdefs.h"
module mod_particleset
  !----------------------------------------------------------------
  ! Particle set module.
  ! This module defines the particleset class which holds all the
  ! particle objects (similarly to the fieldset class).
  ! The particle set class handles releasing new particles and gathering
  ! information to and from the particles.
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
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
    type(t_statevector) :: sv
    integer :: nvars
  contains
    private
    procedure :: init
    procedure :: release
    procedure :: to_statevector
    procedure :: from_statevector
  end type t_particleset
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
    integer :: ipart, sv_index
    integer :: active_particles

    ! Allocate sv
    active_particles = 0
    START_OMP_DO reduction(+:active_particles)
    do ipart = 1, size(this%particles)
      if (this%particles(ipart)%is_active) active_particles = active_particles + 1
    end do
    END_OMP_DO

    allocate (this%sv%trc(active_particles), this%sv%state(active_particles, this%nvars))
    this%sv%nparticles = active_particles
    this%sv%nvars = this%nvars

    ! TODO: OpenMP parallel loop
    sv_index = 0
    do ipart = 1, size(this%particles)
      if (.not. this%particles(ipart)%is_active) cycle
      sv_index = sv_index + 1
      this%sv%trc(sv_index)%ptr => this%particles(ipart)
      this%sv%state(sv_index, :) = this%particles(ipart)%get_state_array()
    end do

    return
  end subroutine to_statevector
  !===========================================
  subroutine from_statevector(this)
    class(t_particleset), intent(in) :: this

    call this%sv%to_tracer()

    return
  end subroutine from_statevector

end module mod_particleset
