#include "cppdefs.h"
module mod_process
  !----------------------------------------------------------------
  ! Base class for a process (advection, diffusion, biofouling etc.)
  ! This class has to be inherited by the specific process and has
  ! two public functions: init and run.
  ! The init function is called once by the kernel at the beginning
  ! of the simulation. The init function must initialize all the variables
  ! needed by the run function. Can have its own namelist.
  ! The run function is called by the kernel at each time step and
  ! returns a derivative of the particle's state vector. Currently,
  ! input and output must be 1D arrays (i.e., one particle at a time).
  !----------------------------------------------------------------
  use mod_common
  use mod_fieldset
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: t_process
  !---------------------------------------------
  type, abstract :: t_process
  contains
    procedure(i_init), deferred :: init
    procedure(i_run), deferred :: run
  end type t_process
  !===================================================
  abstract interface
    subroutine i_init(this)
      import t_process
      class(t_process), intent(inout) :: this
    end subroutine i_init
    function i_run(this, sv, fieldset, time, dt) result(res)
      import t_process, t_fieldset, rk
      class(t_process), intent(in) :: this
      type(t_fieldset), intent(in) :: fieldset
      real(rk), intent(in) :: time, dt
      real(rk), intent(in) :: sv(:)
      real(rk) :: res(size(sv))
    end function i_run
  end interface
  !===================================================
contains
  !===========================================
  ! [subroutine/function definition]

end module mod_process
