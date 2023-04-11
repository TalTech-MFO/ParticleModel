#include "cppdefs.h"
module mod_simulation
  !----------------------------------------------------------------
  ! Simulation module.
  ! This is the main module of the program.
  ! It contains the main (time) loop and handles the setup of the
  ! simulation and the output.
  !----------------------------------------------------------------
  use mod_common
  use mod_datetime
  use mod_particleset
  use mod_fieldset
  use mod_solver
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_simulation
  !---------------------------------------------
  type t_simulation
    type(t_datetime)    :: start_time !< Start time of the simulation
    type(t_datetime)    :: end_time !< End time of the simulation
    type(t_datetime)    :: current_time !< Current time of the simulation
    type(t_particleset) :: particleset !< Particleset to hold the particles
    type(t_fieldset)    :: fieldset !< Fieldset to hold the fields
    type(t_solver)      :: solver !< Solver to move the particles
    real(rk)            :: dt !< Time step
  contains
    procedure :: init
    procedure :: run
    procedure :: finalize
  end type t_simulation
  !===================================================
contains
  !===========================================
  subroutine init(this)
    class(t_simulation), intent(inout) :: this
    !-------------------------------------------
    ! Initialize the simulation.
    !-------------------------------------------
    call this%solver%init()
    call this%solver%kernel%add_process(t_advection)

  end subroutine init
  !===========================================
  subroutine run(this)
    !-------------------------------------------
    ! Run the simulation.
    !-------------------------------------------
    class(t_simulation), intent(inout) :: this
    integer :: itime = 0

    do while (this%current_time < this%end_time)
      ! Update fields
      ! Release particles
      call particleset%tracer_to_statevector()
      call this%solver%run(particleset%sv)
      call particleset%statevector_to_tracer()
      ! Update time
      call this%current_time%update(this%dt)
      ! Write output
      ! Write restart
      itime = itime + 1
    end do

  end subroutine run
  !===========================================
  subroutine finalize(this)
    class(t_simulation), intent(inout) :: this
    !-------------------------------------------
    ! Finalize the simulation.
    !-------------------------------------------

  end subroutine finalize

end module mod_simulation
