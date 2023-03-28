#include "cppdefs.h"
module mod_solver
  !----------------------------------------------------------------
  ! Solver module.
  ! Defines a solver class which contains the integration methods
  ! (euler, rk2, rk4, etc.).
  ! The solver's run method runs the particle loop and calls the kernel
  ! which computes the derivatives of the particle variables.
  !----------------------------------------------------------------
  use mod_precdefs
  use mod_errors
  use mod_kernel
  use mod_particle
  use mod_statevector
  use mod_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_solver
  !---------------------------------------------
  type t_solver
    integer :: solver_type = 1
    type(t_kernel) :: kernel
  contains
    procedure :: init
    procedure :: run
    procedure :: run_euler, run_rk2, run_rk4
  end type t_solver
  !===================================================
contains
  !===========================================
  subroutine init(this, solver_type)
    class(t_solver), intent(inout) :: this
    integer, intent(in) :: solver_type

    this%solver_type = solver_type
    call this%kernel%init()

    return
  end subroutine init
  !===========================================
  subroutine run(this, sv, fieldset, time, dt)
    class(t_solver), intent(in) :: this
    class(t_statevector), intent(inout) :: sv
    class(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk), intent(in) :: dt

    select case (this%solver_type)
    case (1)
      call this%run_euler(sv, fieldset, time, dt)
    case (2)
      call this%run_rk2(sv, fieldset, time, dt)
    case (3)
      call this%run_rk4(sv, fieldset, time, dt)
    case default
      call this%run_euler(sv, fieldset, time, dt)
    end select

    return
  end subroutine run
  !===========================================
  subroutine run_euler(this, sv, fieldset, time, dt)
    class(t_solver), intent(in) :: this
    type(t_statevector), intent(inout) :: sv
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk), intent(in) :: dt
    integer :: ipart

    ! TODO: OpenMP parallel loop
    do ipart = 1, sv%nparticles
      sv%state(ipart, :) = sv%state(ipart, :) + this%kernel%run(sv%state(ipart, :), fieldset, time, dt) * dt
    end do

    return
  end subroutine run_euler
  !===========================================
  subroutine run_rk2(this, sv, fieldset, time, dt)
    class(t_solver), intent(in) :: this
    type(t_statevector), intent(inout) :: sv
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk), intent(in) :: dt
    integer :: ipart
    type(t_statevector) :: k1

    call k1%copy(sv)

    ! TODO: OpenMP parallel loop
    do ipart = 1, sv%nparticles
      k1%state(ipart, :) = sv%state(ipart, :) + this%kernel%run(sv%state(ipart, :), fieldset, time, dt) * dt
      sv%state(ipart, :) = sv%state(ipart, :) + &
                           0.5_rk * (k1%state(ipart, :) + this%kernel%run(k1%state(ipart, :), fieldset, time, dt) * dt)
    end do

    call k1%clean()

    return
  end subroutine run_rk2
  !===========================================
  subroutine run_rk4(this, sv, fieldset, time, dt)
    class(t_solver), intent(in) :: this
    type(t_statevector), intent(inout) :: sv
    type(t_fieldset), intent(in) :: fieldset
    real(rk), intent(in) :: time
    real(rk), intent(in) :: dt

    call throw_error("solver :: run_rk4", "Not implemented yet.")
  end subroutine run_rk4

end module mod_solver
