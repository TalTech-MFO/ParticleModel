#include "cppdefs.h"
#include "particle.h"
#include "ompdefs.h"
module mod_loop
  !----------------------------------------------------------------
  ! Main loop
  !----------------------------------------------------------------
#ifdef USE_OMP
  use omp_lib
#endif
  use mod_precdefs
  use mod_errors
  use mod_params, only: do_velocity, do_diffusion, do_biofouling, run_3d, advection_method
  use mod_advection
  use mod_vertical_motion
  use mod_diffusion
  use mod_biofouling, only: biofouling
  use field_vars, only: fieldset
  use mod_particle, only: t_particle
  use mod_domain_vars, only: domain
  use mod_particle_vars, only: particles, inputstep, &
                               max_age, runparts, kill_beached, kill_boundary, release_particles
  use time_vars, only: theDate, run_start_dt, run_end_dt, dt
  use mod_output, only: outputstep, restartstep, write_data, &
#ifdef OUT_SUPPORT_ACTIVE
                        write_data_only_active, &
#endif
                        write_restart, write_all_particles, write_active_particles
  use postprocess_vars, only: postprocessor, enable_postprocessing
  implicit none
  private
  !===================================================
  !---------------------------------------------
  public :: loop
  !---------------------------------------------
  !===================================================
contains
  !===========================================
  subroutine loop
    integer           :: itime = 0
    integer           :: ipart
    real(rk)          :: time
    character(len=8)  :: d
    character(len=10) :: t
#ifndef SAY_LESS
    integer           :: active_particles, inactive_particles
    character(len=LEN_CHAR_L) :: info
#endif

    FMT1, "======== Starting time loop ========"

#ifdef USE_OMP
    FMT2, "Using OpenMP with ", omp_get_max_threads(), " threads"
#endif

    ! Read appropriate fields:
    call fieldset%read_first_timesteps(run_start_dt)

    ! Start integration loop
    do while (theDate < run_end_dt)

      !   - update fields
      call fieldset%update(theDate)
      time = fieldset%get_time(theDate)

      !   - reset postprocessor fields (if needed)
      if (enable_postprocessing) then 
        call postprocessor%reset_measures()
      endif

      !   - release particles
      call release_particles(itime, theDate, fieldset, time)

#ifndef SAY_LESS
      if (mod(itime, PROGRESSINFO) .eq. 0) then
        call date_and_time(date=d, time=t)
  write (info, "(a,i22,a)") "| "//theDate%nice_format()//" | itime = ", itime, " | Time: "//t(1:2)//":"//t(3:4)//":"//t(5:10)//"  |"
        FMT2, LINE, LINE, LINE
        FMT2, trim(info)
        FMT2, LINE, LINE, LINE
        if (itime .ne. 0) then
          active_particles = 0
          inactive_particles = 0

          START_OMP_DO reduction(+:active_particles, inactive_particles)
          do ipart = 1, runparts
            if (particles(ipart)%is_active) then
              active_particles = active_particles + 1
            else
              inactive_particles = inactive_particles + 1
            end if
          end do
          END_OMP_DO
          FMT2, "| ", active_particles, " active particles"
          FMT2, "| ", inactive_particles, " inactive particles"
        end if
      end if
#endif

      !---------------------------------------------
      ! Start particle loop
      START_OMP_DO shared(particles, fieldset, domain) schedule(guided)
      do ipart = 1, runparts
        DBG, LINE
        DBG, "Particle nr.", ipart
        DBG, LINE
        !---------------------------------------------
        ! Skip inactive particles
        if (.not. particles(ipart)%is_active) cycle

        !---------------------------------------------
        ! Advect only if the particle is alive (is_active=.true.) and active (state=3/4)
        if (particles(ipart)%state >= ST_SUSPENDED) then
          ! - do advection
          call advect(particles(ipart), fieldset, time, advection_method, run_3d)
          ! - do biofouling
          if (do_biofouling) call biofouling(particles(ipart), fieldset, time)
          ! - do vertical velocity
          if (do_velocity .and. run_3d) call vertical_velocity(particles(ipart), fieldset, time)
          ! - do diffusion
          if (do_diffusion) call diffuse(particles(ipart), fieldset, time, run_3d)
        end if

        !---------------------------------------------
        ! Update particles
#ifdef DEBUG
        DBG, "Particle before updating"
        call particles(ipart)%print_info()
#endif
        call particles(ipart)%update(fieldset, time, dt)

#ifndef USE_VECTORIZED_POSTPROCESSING
        !---------------------------------------------
        ! Postprocessing
        if (enable_postprocessing) then 
          call postprocessor%after_timestep(particles(ipart))
        end if
#endif

      end do
      END_OMP_DO

#ifdef USE_VECTORIZED_POSTPROCESSING
      !---------------------------------------------
      ! Postprocessing
      if (enable_postprocessing) then
        call postprocessor%after_timestep(particles(:runparts))
      end if
#endif

      !---------------------------------------------
      ! Update time
      call theDate%update(dt)

      !---------------------------------------------
      ! Write output
      if ((mod(itime, outputstep) .eq. 0) .and. (runparts .gt. 0)) then
        if (write_all_particles) call write_data(runparts)
#ifdef OUT_SUPPORT_ACTIVE
        if (write_active_particles) call write_data_only_active(runparts)
#endif
      end if

      if (enable_postprocessing) then 
        call postprocessor%save(itime, theDate)
      end if

      !---------------------------------------------
      ! Write restart (save after updating the date so it could be used as initial state later)
      if ((restartstep > 0) .and. (mod(itime, restartstep) == 0)) then
        call write_restart(runparts)
      end if

      !---------------------------------------------
      itime = itime + 1

    end do

    !---------------------------------------------
    ! Write restart at end of simulation
    if (restartstep == 0) then
      call write_restart(runparts)
    end if
    call postprocessor%write_restart(theDate)

    call date_and_time(date=d, time=t)
    FMT2, LINE
    FMT2, t(1:2), ":", t(3:4), ":", t(5:10)
    FMT2, "Finished all time steps (", itime, ")"
    FMT2, LINE

    return
  end subroutine loop
end module mod_loop
