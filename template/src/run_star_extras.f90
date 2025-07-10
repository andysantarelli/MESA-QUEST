! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   This file is part of MESA.
!
!   MESA is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   MESA is distributed in the hope that it will be useful, 
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************
!
! This module implements black hole accretion physics in stellar evolution models.
! It includes Bondi accretion with cavity formation and sound speed smoothing
! for numerical stability.
!
! ***********************************************************************
 
      module run_star_extras
      
      use star_lib
      use star_def
      use const_def
      use math_lib
      use utils_lib, only: mesa_error
      use auto_diff
      
      implicit none
      
      ! Smoothing parameters for sound speed stabilization
      real(dp), parameter :: temporal_smooth_factor = 0.0d0   ! Weight for temporal averaging
      integer, parameter :: n_zones_smooth = 50               ! Number of zones for spatial averaging
      
      ! Storage for smoothing algorithm
      real(dp), save :: c_s_smooth_prev = 0d0
      logical, save :: first_call = .true.
      
      contains

      include "other_cgrav.inc"
      
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! Set procedure pointers for custom routines
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items
         s% other_cgrav => my_other_cgrav

      end subroutine extras_controls
      
      
      ! Calculate spatially smoothed sound speed for numerical stability
      real(dp) function get_smoothed_sound_speed(s) result(c_s_smooth)
         type (star_info), pointer :: s
         real(dp) :: c_s_spatial, sum_cs
         integer :: k, k_start, k_end, nz
         
         nz = s% nz
         
         ! Spatial averaging over innermost zones
         k_start = max(1, nz - n_zones_smooth + 1)
         k_end = nz
         
         sum_cs = 0d0
         do k = k_start, k_end
            sum_cs = sum_cs + s% csound(k) 
         end do
         
         c_s_spatial = sum_cs / n_zones_smooth
         
         ! Apply temporal averaging if enabled
         if (first_call) then
            c_s_smooth = c_s_spatial
            first_call = .false.
         else
            c_s_smooth = (1d0 - temporal_smooth_factor) * c_s_spatial + &
                         temporal_smooth_factor * c_s_smooth_prev
         end if
         
         ! Store for next timestep
         c_s_smooth_prev = c_s_smooth
         
      end function get_smoothed_sound_speed


      ! Rate-limited Bondi radius update
      real(dp) function get_bondi_radius(s, R_B_new) result(R_B_limited)
         type (star_info), pointer :: s
         real(dp), intent(in) :: R_B_new
         real(dp) :: max_change, max_fractional_change
         real(dp) :: dt, c_s

         dt  = s% dt  ! time step (s)
         c_s = get_smoothed_sound_speed(s)  ! sound speed (cm/s)
         
         if (first_bondi_call) then
            R_B_limited = R_B_new
            first_bondi_call = .false.
         else
            max_fractional_change = dt * c_s / (2 * R_B_new)
            max_change = max_fractional_change * R_B_prev
            
            ! Limit the change
            if (R_B_new > R_B_prev + max_change) then
               R_B_limited = R_B_prev + max_change
            else if (R_B_new < R_B_prev - max_change) then
               R_B_limited = R_B_prev - max_change
            else
               R_B_limited = R_B_new
            end if
         end if
         
         R_B_prev = R_B_limited
         
      end function get_bondi_radius
      
      
      ! Adjust radial coordinates to match new center radius
      subroutine do1_relax_R_center(s, new_Rcenter, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: new_Rcenter ! cm
         integer, intent(out) :: ierr
         real(dp) :: dm, rho, dr3, rp13
         integer :: k
         
         ierr = 0
         s% R_center = new_Rcenter
         
         ! Adjust lnR's to maintain cell densities
         rp13 = s% R_center**3
         do k = s% nz, 1, -1
            dm = s% dm(k)
            rho = s% rho(k)
            dr3 = dm / (rho * four_thirds_pi) ! Cell volume
            s% xh(s% i_lnR, k) = log(rp13 + dr3) * one_third
            rp13 = rp13 + dr3
         end do
         
      end subroutine do1_relax_R_center
      
      
      ! Main black hole accretion physics routine
      subroutine black_hole_accretion(id, s, startup, ierr)
          integer, intent(in) :: id
          logical, intent(in) :: startup
          type (star_info), pointer :: s
          integer, intent(out) :: ierr
          
          ! Physical variables
          real(dp) :: G, c2, c_s, c_s_smooth, rho, gamma1, opacity, dt
          real(dp) :: nabla_ad, P_rad, P_gas
          real(dp) :: L_BH, M_BH, M_BH_new, M_cav, M_rat
          real(dp) :: R_B, R_B_smooth, R_B_raw, R
          real(dp) :: M_dot, M_dot_BH, dm, M_dot_Bondi
          real(dp) :: rad_eff, con_eff, timestep_factor
          real(dp) :: core_avg_rho, core_avg_eps, new_core_mass
          
          ! Control parameters
          rad_eff = s% x_ctrl(1)        ! Radiative efficiency (epsilon)
          con_eff = s% x_ctrl(2)        ! Convective efficiency (eta)
          timestep_factor = s% x_ctrl(3) ! Timestep scaling factor
          
          ! Extract stellar properties at center
          dt      = s% dt
          G       = s% cgrav(s% nz)
          c_s     = s% csound(s% nz)
          rho     = s% rho(s% nz)
          opacity = s% opacity(s% nz)
          P_rad   = s% prad(s% nz)
          P_gas   = s% pgas(s% nz)
          gamma1  = s% gamma1(s% nz)
          R       = 10**(s% log_surface_radius) * Rsun
          nabla_ad = 1d0 - 1d0 / gamma1
          c2 = clight**2

          M_BH = s% xtra(1) ! Black hole mass in grams
          M_rat = M_BH / s% mstar
          
          ! Get smoothed sound speed for stability
          c_s_smooth = get_smoothed_sound_speed(s)
          
          ! Calculate Bondi accretion rate with smoothed sound speed
          M_dot_Bondi = 16d0 * pi * rho * (G * M_BH)**2 / c_s_smooth / c2 * &
                        con_eff / rad_eff * (1d0 - rad_eff) / gamma1
          L_BH = rad_eff / (1d0 - rad_eff) * M_dot_Bondi * c2
          
          M_dot = L_BH / (rad_eff * c2)
          dm = (1d0 - rad_eff) * M_dot * dt
          M_BH_new = M_BH + dm
          
          ! Calculate Bondi radius
          R_B_smooth = 2d0 * G * M_BH_new / c_s_smooth**2
          R_B_raw = 2d0 * G * M_BH_new / c_s**2
          R_B = get_bondi_radius(s, R_B_smooth)
          
          ! Calculate cavity mass and core properties
          M_cav = 8d0 * pi / 3d0 * rho * R_B**3
          new_core_mass = (M_BH_new + M_cav) / Msun
          core_avg_eps = L_BH / (new_core_mass * Msun)
          core_avg_rho = (new_core_mass * Msun) / (4d0 / 3d0 * pi * R_B**3)

          ! Set adaptive timestep
          s% max_timestep = timestep_factor * M_BH / ((1d0 - rad_eff) * M_dot)
          
          ! Store results in xtra array for output
          s% xtra(1)  = M_BH_new
          s% xtra(2)  = L_BH
          s% xtra(3)  = R_B
          s% xtra(4)  = R_B_raw
          s% xtra(5)  = M_dot
          s% xtra(6)  = safe_log10(dm) - safe_log10(dt)
          s% xtra(7)  = rad_eff
          s% xtra(8)  = opacity
          s% xtra(9)  = M_cav
          s% xtra(10)  = P_rad
          s% xtra(11) = P_gas
          s% xtra(12) = nabla_ad
          s% xtra(13) = M_dot_Bondi
          s% xtra(14) = L_BH
          s% xtra(15) = c_s_smooth
          s% xtra(16) = c_s

          ! Apply core modifications
          if (startup) then
              call star_relax_core( &
                  id, new_core_mass, s% job% dlg_core_mass_per_step, &
                  s% job% relax_core_years_for_dt, &
                  core_avg_rho, core_avg_eps, ierr)
          else
              s% M_center = new_core_mass * Msun
              s% mstar = s% mstar - rad_eff * M_dot * dt
              s% xmstar = s% mstar - s% M_center
              s% L_center = L_BH
              call do1_relax_R_center(s, R_B, ierr)
          end if

          
          ! Output diagnostics
          write(*, '(a)') '--- Black Hole Properties ---'
          write(*, '(a, f12.6)') 'M/M_sun: ',    s% mstar / Msun
          write(*, '(a, f12.6)') 'M_BH/M_sun: ', M_BH_new / Msun
          write(*, '(a, f12.6)') 'L_BH/L_sun: ', L_BH / Lsun
          write(*, '(a, f12.6)') 'R_B/R_star: ', R_B / R
          write(*, '(a, f12.6)') 'Radiative efficiency: ', rad_eff
          write(*, '(a, es12.4)') 'Sound speed (raw): ', c_s
          write(*, '(a, es12.4)') 'Sound speed (smooth): ', c_s_smooth
          write(*, '(a, f12.6)') 'R_B (smooth): ', R_B / Rsun
          write(*, '(a)') '-----------------------------'
          
      end subroutine black_hole_accretion
          

      subroutine extras_startup(id, restart, ierr)
          integer, intent(in) :: id
          logical, intent(in) :: restart
          integer, intent(out) :: ierr
          type (star_info), pointer :: s

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          ! Initialize black hole mass and perform startup calculations
          if (s% x_logical_ctrl(1)) then
              s% xtra(1) = s% job% new_core_mass * Msun
              call black_hole_accretion(id, s, .true., ierr) 
          end if

      end subroutine extras_startup
      
      
      integer function extras_start_step(id)
          integer, intent(in) :: id
          integer :: ierr
          type (star_info), pointer :: s

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return
          extras_start_step = 0

          ! Update black hole properties at start of each step
          if (s% x_logical_ctrl(1)) then
              call black_hole_accretion(id, s, .false., ierr) 
          end if
          
      end function extras_start_step
      
      
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
         
         ! Terminate if black hole mass exceeds stellar mass
         if (s% xtra(1) >= s% mstar) then
            extras_check_model = terminate
            write(*, *) 'Terminating: M_BH >= M_star'
            termination_code_str(t_xtra1) = 'black hole mass exceeds stellar mass'
         end if

         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 16
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! Initialize arrays
         names(1:16) = 'empty'
         vals(1:16) = -1d99
         
         ! Set black hole history columns
         if (s% x_logical_ctrl(1)) then 
             names(1) = "M_BH"
             vals(1) = s% xtra(1) / Msun
             names(2) = "L_BH"
             vals(2) = s% xtra(2) / Lsun
             names(3) = "R_B"
             vals(3) = s% xtra(3) / Rsun
             names(4) = "R_B_raw"
             vals(4) = s% xtra(4) / Rsun
             names(5) = "M_dot"
             vals(5) = s% xtra(4) / Msun
             names(6) = "log10_dm_dt"
             vals(6) = s% xtra(5)
             names(7) = "rad_eff"
             vals(7) = s% xtra(6)
             names(8) = "kap_center"
             vals(8) = s% xtra(7)
             names(9) = "M_cav"
             vals(9) = s% xtra(8) / Msun
             names(10) = "prad_center"
             vals(10) = s% xtra(9)
             names(11) = "pgas_center"
             vals(11) = s% xtra(10)
             names(12) = "nabla_ad_center"
             vals(12) = s% xtra(11)
             names(13) = "M_dot_Bondi" 
             vals(13) = s% xtra(12)
             names(14) = "L_BH"
             vals(14) = s% xtra(13) / Lsun
             names(15) = "cs_smooth"
             vals(15) = s% xtra(14)
             names(16) = "cs_center"
             vals(16) = s% xtra(15)
         end if

      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 1
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'M_BH'
         vals(1) = s% xtra(1) / Msun

      end subroutine data_for_extra_profile_header_items


      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         real(dp) :: current_mass, current_interval
         integer, save :: last_mass_checkpoint = 0

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Save profiles at 100 Msun intervals of black hole mass
         current_mass = s% xtra(1) / Msun
         current_interval = int(current_mass / 100.0d0)
   
         if (current_interval > last_mass_checkpoint) then
            s% need_to_save_profiles_now = .true.
            last_mass_checkpoint = current_interval
            write(*, '(a, f8.2, a)') 'Saving profile at BH mass: ', current_mass, ' Msun'
         end if

         extras_finish_step = keep_going
         
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
         
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve

      end module run_star_extras
