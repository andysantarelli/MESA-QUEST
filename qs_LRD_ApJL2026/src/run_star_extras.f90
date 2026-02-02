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
      use utils_lib, only: mesa_error, mkdir
      use auto_diff
      
      implicit none
      private
      public :: extras_controls  ! Add this line
      include "test_suite_extras_def.inc"
      

      ! defs

      real(dp) :: M_BH_new, L_BH, R_B, R_B_raw
      real(dp) :: M_dot, log_dm_dt, rad_eff, opacity
      real(dp) :: M_cav, P_rad, P_gas, nabla_ad
      real(dp) :: M_dot_Bondi, L_Bondi, c_s_smooth, c_s, R_sc
           
      contains

      include "other_cgrav.inc"
      include "quasi_star.inc"

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
         s% other_photo_write => photo_write
         s% other_photo_read => photo_read

      end subroutine extras_controls


      real(dp) function get_R_sc(s) result(R_i)
         type (star_info), pointer :: s
         real(dp) :: K_i, M_rat, R, t, m_lo, m_hi, xi_lo, xi_hi, K_lo, K_hi
         integer :: i
         integer, parameter :: nrows = 465
         real(dp), dimension(3, nrows) :: x

         ! Table layout
         ! x(1,i) = K_i
         ! x(2,i) = M_ratio = M_BH / M_star
         ! x(3,i) = xi_r    = R / R_i  (so R_i = R / xi_r)

         R     = s% r(1)                ! total radius (cm)
         M_rat = s% xtra(1) / s% mstar       ! M_BH / M_star, not mbh new?

         open (99, file = "Ki_tab.dat", status = "old")
         do i = 1, nrows
            read (99,*) x(1,i), x(2,i), x(3,i)
         end do
         close (99)

         ! use first row below table range
         if (M_rat <= x(2,1)) then
            K_i = x(1,1)
            R_i = R / x(3,1)
            return
         end if

         ! above table range
         if (M_rat >= x(2,nrows)) then
            K_i = x(1,nrows)
            R_i = R / x(3,nrows)
            return
         end if

         ! dind bracketing interval [i-1, i] with x(2,i-1) <= M_rat <= x(2,i)
         do i = 2, nrows
            if (M_rat <= x(2,i)) then
               m_lo = x(2,i-1); m_hi = x(2,i)
               t    = (M_rat - m_lo) / (m_hi - m_lo)

               ! Linear interpolation for K_i and xi_r
               K_lo = x(1,i-1); K_hi = x(1,i)
               xi_lo = x(3,i-1); xi_hi = x(3,i)

               K_i = (1.0_dp - t)*K_lo + t*K_hi
               R_i = R / ((1.0_dp - t)*xi_lo + t*xi_hi)
               return
            end if
         end do

         ! fallback to last row
         K_i = x(1,nrows)
         R_i = R / x(3,nrows)
      end function get_R_sc
      
      ! Adjust radial coordinates to match new center radius
      subroutine do1_relax_R_center(s, new_Rcenter, ierr)
         type (star_info), pointer :: s
         real(dp), intent(in) :: new_Rcenter ! cm
         integer, intent(out) :: ierr
         real(dp) :: dm, rho, dr3, rp13, r300, r3p1, vol
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
            s%lnR(k) = s%xh(s%i_lnR, k)
            s%r(k) = exp(s%lnR(k))
            rp13 = rp13 + dr3
         end do

!         ! adjust densities
!         do k = 1, s%nz
!            r300 = pow3(s%r(k))
!            if (k < s%nz) then
!               r3p1 = pow3(s%r(k + 1))
!            else
!               r3p1 = pow3(s%r_center)
!            end if
!            vol = (4d0*pi/3d0)*(r300 - r3p1)
!            s%rho(k) = s%dm(k)/vol
!            s%lnd(k) = log(s%rho(k))
!            s%xh(s%i_lnd, k) = s%lnd(k)
!            if (is_bad(s%lnd(k))) then
!               write (*, 2) 'bad lnd vol dm r300 r3p1', k, s%lnd(k), vol, s%dm(k), r300, r3p1
!               call mesa_error(__FILE__, __LINE__, 'remesh for quasi star')
!            end if
!         end do


         
      end subroutine do1_relax_R_center
      
      
      ! Main black hole accretion physics routine
      subroutine black_hole_accretion(id, s, startup, restart, ierr)
          integer, intent(in) :: id
          logical, intent(in) :: startup, restart
          type (star_info), pointer :: s
          integer, intent(out) :: ierr
          
          integer :: k
          ! Physical variables
          real(dp) :: G, c2, rho, dt, quasi_star_mass
          real(dp) :: L_Edd, gamma1
          real(dp) :: M_BH, M_rat
          real(dp) :: R_B, R
          real(dp) :: M_dot, M_dot_BH, dm
          real(dp) :: rad_eff, con_eff, timestep_factor, alpha
          real(dp) :: core_avg_rho, core_avg_eps, new_core_mass, L_center_prev, dL
          
          ! Control parameters
          rad_eff = s% x_ctrl(1)         ! Radiative efficiency (epsilon)
          con_eff = s% x_ctrl(2)         ! Convective efficiency (eta)
          timestep_factor = s% x_ctrl(3) ! Timestep scaling factor
          alpha = s% x_ctrl(4)           ! Accretion factor (alpha)
          quasi_star_mass = s% m(1)/Msun ! Quasi-star mass (Msun)

          !write (*,*) 'quasi_star_mass', quasi_star_mass
          ! Stellar properties at center

          dt      = s% dt
          G       = s% cgrav(s% nz)!0.5d0*(s% cgrav(s% nz)+s% cgrav(s% nz-1))
          c_s     = s% csound(s% nz)!s% csound_face(s% nz)
          rho     = s%rho(s% nz)!s% rho_face(s% nz)
          opacity = s% opacity(s% nz)!get_kap_face(s,s% nz)!s%opacity(s% nz)!get_kap_face(s,s% nz) ! 0.2d0*(1d0 + s% X(s%nz))
          gamma1  = s%gamma1(s% nz)!get_gamma1_face(s,s% nz)
          R       = s% r(s% nz) !0.5d0*(s% r(s% nz)+s% r(s%nz-1))
          c2 = clight**2

          M_BH = M_BH_new ! s% xtra(1) ! Black hole mass (g)
          M_rat = M_BH / s% mstar
          
          ! Calculate Bondi accretion rate with smoothed sound speed
          M_dot_Bondi = 16d0 * pi * rho * pow2(G * M_BH) / c_s / c2 * &
                        con_eff / rad_eff * (1d0 - rad_eff) / gamma1
          L_Bondi = rad_eff / (1d0 - rad_eff) * M_dot_Bondi * c2
          L_Edd = 4*pi * clight * G * quasi_star_mass * Msun / opacity  ! Eddington rate for entire object, erg/s

          L_BH = alpha * L_Bondi
          if (s% x_logical_ctrl(3)) L_BH = alpha * L_Edd
          
          M_dot = L_BH / (rad_eff * c2)
          dm = (1d0 - rad_eff) * M_dot * dt
          M_BH_new = M_BH + dm
          
          ! Calculate inner boundaries
          R_B = 2d0 * G * M_BH_new / pow2(c_s)
          R_sc = get_R_sc(s) 
          
          ! Calculate cavity mass and core properties
          M_cav = 8d0 * pi / 3d0 * rho * pow3(R_B) ! rho prop to r**-3/2
          if (s% x_logical_ctrl(2)) M_cav = 8d0 * pi / 5d0 * rho * pow3(R_sc) ! rho prop to r**-1/2
          new_core_mass = (M_BH_new + M_cav) / Msun
          
          core_avg_eps = L_BH / (new_core_mass * Msun)
          core_avg_rho = (new_core_mass * Msun) / (4d0 / 3d0 * pi * pow3(R_B))
          if (s% x_logical_ctrl(2)) core_avg_rho = (new_core_mass * Msun) / (4d0 / 3d0 * pi * pow3(R_sc))

          ! Set adaptive timestep
          s% max_timestep = timestep_factor * M_BH / ((1d0 - rad_eff) * M_dot)
          
          ! Store results in xtra array for output
          s% xtra(1)  = M_BH_new
          s% xtra(2)  = L_BH
          s% xtra(3)  = R_B
          s% xtra(4)  = M_dot
          log_dm_dt  = safe_log10(dm) - safe_log10(dt)
          s% xtra(5)  = log_dm_dt
          s% xtra(6)  = rad_eff
          s% xtra(7)  = opacity
          s% xtra(8)  = M_cav
          s% xtra(9) = M_dot_Bondi
          s% xtra(10) = L_Bondi
          s% xtra(11) = c_s
          s% xtra(12) = R_sc

          ! Apply core modifications
          if (startup .and. .not. restart) then
              call star_write_model(id, "initial_model.mod", ierr)
!              call star_relax_core( &
!                  id, new_core_mass, s% job% dlg_core_mass_per_step, &
!                  s% job% relax_core_years_for_dt, &
!                  core_avg_rho, core_avg_eps, ierr)

             s% M_center = new_core_mass !M_BH
             s% L_center = L_BH

            if (s% x_logical_ctrl(2)) then
                call do1_relax_R_center(s, R_sc, ierr)
                !call star_relax_R_center(id, R_sc, 1d-2, 1d-5, ierr)
            else
                call do1_relax_R_center(s, R_B, ierr)  ! Use smoothed R_B
                !call star_relax_R_center(id, R_B, 1d-2, 1d-5, ierr)
            end if

          else

            
              s% M_center = new_core_mass * Msun
              s% mstar = s% mstar - rad_eff * M_dot * dt
              s% xmstar = s% mstar - s% M_center
              s% L_center = L_BH

!            !nz = s% nz
!            L_center_prev = s% L_center
!            s% L_center = L_BH
!            dL = L_BH - L_center_prev
!            do k=1,s%nz
!               s% xh(s%i_lum,k) = s% xh(s%i_lum,k) + dL
!            end do


            if (s% x_logical_ctrl(2)) then
                call do1_relax_R_center(s, R_sc, ierr)
                !s% R_center = R_sc
                !call star_relax_R_center(id, R_sc, 1d-2, 1d-5, ierr)
            else
            !                  s% R_center = R_B
                call do1_relax_R_center(s, R_B, ierr)  ! Use smoothed R_B
                !call star_relax_R_center(id, R_B, 1d-2, 1d-5, ierr)
            end if

          end if
          
          ! Output diagnostics
          write(*, '(a)') '--- Black Hole Properties ---'
          write(*,*) 'M/M_sun: ',    s% mstar / Msun
          write(*,*) 'M_BH/M_sun: ', M_BH_new / Msun
          write(*,*) 'Core_mass/M_sun: ', new_core_mass*Msun / s% mstar
          write(*,*) 'M_BH/M_star: ', M_BH_new / s% mstar
          write(*, '(a, f12.6)') 'log(L_BH/L_sun): ', max(-99,log10(L_BH/Lsun))! / Lsun
          write(*, '(a, f12.6)') 'R_B/R_star: ', R_B / s% r(1)!R
          write(*, '(a, f12.6)') 'R_sc/R_star: ', R_sc / s% r(1)
          write(*, '(a, f12.6)') 'Radiative efficiency: ', rad_eff
          write(*, '(a, es12.4)') 'Sound speed (raw): ', c_s
          write(*, '(a)') '-----------------------------'
          
      end subroutine black_hole_accretion

      subroutine superedd_winds(id, s, ierr)
          integer, intent(in) :: id
          type (star_info), pointer :: s
          integer, intent(out) :: ierr
          real(dp) :: M_dot_superedd, dm_superedd

          M_dot_superedd = 1.4d-4 * (s% mstar / Msun)**0.96 * (M_BH_new / Msun)**0.17 ! Msun/yr
          dm_superedd = M_dot_superedd * Msun/secyer * s% dt

          s% mstar = s% mstar - dm_superedd   ! remove wind mass
          s% xmstar = s% mstar - s% M_center  ! correct xmstar

      end subroutine superedd_winds 

      subroutine extras_startup(id, restart, ierr)
          integer, intent(in) :: id
          logical, intent(in) :: restart
          integer, intent(out) :: ierr
          type (star_info), pointer :: s

          ierr = 0
          call star_ptr(id, s, ierr)
          if (ierr /= 0) return

          if (.not. restart) then
             M_BH_new = 0d0
             write(*,*) 'M_BH_new', M_BH_new
             L_BH = 0d0
             R_B = 0d0
             R_B_raw = 0d0
             M_dot= 0d0
             log_dm_dt= -99d0
             rad_eff= 0d0
             opacity = 0d0
             M_cav= 0d0
             P_rad= 0d0
             P_gas= 0d0
             nabla_ad = 0d0
             M_dot_Bondi= 0d0
             L_Bondi= 0d0
             c_s_smooth= 0d0
             c_s= 0d0
             R_sc = 0d0
          end if

          ! Initialize black hole mass and perform startup calculations
          if (.not. restart .and. s% x_logical_ctrl(1)) then
              M_BH_new = s% job% new_core_mass * Msun
               write(*,*) 's% xtra(1)', s% xtra(1)
               write(*,*) 'M_BH_new', M_BH_new

              call black_hole_accretion(id, s, .true., restart, ierr)
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


         ! helps set up initial model with the right entropy, and prevents model from leaving opacity table grids.
         if (s%T(s%nz)> 1d6 ) then
!            s% inject_uniform_extra_heat = 1d8
!            s% min_q_for_uniform_extra_heat = 0
!            s% max_q_for_uniform_extra_heat = 1
         else
            s% inject_uniform_extra_heat = 0d0
            s% use_other_cgrav = .true.
         end if

          ! Update black hole properties at start of each step
          if (s% x_logical_ctrl(1)) then
              call black_hole_accretion(id, s, .false.,.false., ierr)
              if (s% x_logical_ctrl(4)) then
                  call superedd_winds(id, s, ierr)
              end if
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

         if (M_BH_new >= s% mstar) then
            extras_check_model = terminate
            write(*, *) 'Terminating: M_BH >= M_star'
            termination_code_str(t_xtra1) = 'black hole mass exceeds stellar mass'
         end if


         if (s% x_logical_ctrl(2)) then
            if (R_sc - s%r(1) > 1d-2) then
               extras_check_model = terminate
               write(*, *) 'Terminating R_sc > R_surface'
            end if
         else
            if (R_B - s%r(1) > -1d-2) then
               extras_check_model = terminate
               write(*, *) 'Terminating R_B > R_surface'
            end if
         end if

         if (M_BH_new / Msun >= s% x_ctrl(6) .and. s% x_ctrl(6) > 0) then
            extras_check_model = terminate
            write(*, *) 'Terminating: M_BH >= mass limit'
         end if

         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
         
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr, n
         type (star_info), pointer :: s
         character(len=100), allocatable :: strings(:)
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 12
         write(*,*) 'Number of extra history columns:', how_many_extra_history_columns
      end function how_many_extra_history_columns

      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: i
         logical :: make_sed
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! Set black hole history columns
         if (s% x_logical_ctrl(1)) then 
             names(1)  = 'M_BH'
             names(2)  = 'L_BH'
             names(3)  = 'R_B'
             names(4)  = 'M_dot'
             names(5)  = 'log_dm_dt'
             names(6)  = 'rad_eff'
             names(7)  = 'opacity'
             names(8)  = 'M_cav'
             names(9)  = 'M_dot_Bondi'
             names(10) = 'L_Bondi'
             names(11) = 'c_s'
             names(12) = 'R_sc'
             vals(1)  = s% xtra(1) / Msun
             vals(2)  = s% xtra(2) / Lsun
             vals(3)  = s% xtra(3) / Rsun
             vals(4)  = s% xtra(4)
             vals(5)  = s% xtra(5)
             vals(6)  = s% xtra(6)
             vals(7)  = s% xtra(7)
             vals(8)  = s% xtra(8) / Msun
             vals(9)  = s% xtra(9)
             vals(10) = s% xtra(10) / Lsun
             vals(11) = s% xtra(11)
             vals(12) = s% xtra(12) / Rsun
         end if
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer :: k
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         names(1) = 'cgrav_ratio'
         do k=1,s% nz
            vals(k,1) = s% cgrav(k)/standard_cgrav
         end do
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
         vals(1) = M_BH_new / Msun

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
         current_mass = M_BH_new / Msun
         current_interval = int(current_mass / s% x_ctrl(5))  ! x_ctrl(5) is mass interval in Msun
   
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

      subroutine photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         call quasi_star_photo_write(id, iounit)
      end subroutine photo_write

      subroutine photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         call quasi_star_photo_read(id, iounit, ierr)
      end subroutine photo_read

! additional functions for getting face values

      function get_kap_face(s,k) result(kap_face)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: kap_face
         real(dp) :: alfa, beta
         if (k == 1) then
            kap_face = s% opacity(k)
            return
         end if
         call get_face_weights(s, k, alfa, beta)
         kap_face = alfa*s% opacity(k) + beta*s%opacity(k-1)
      end function get_kap_face

      function get_gamma1_face(s,k) result(gamma1_face)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp) :: gamma1_face
         real(dp) :: alfa, beta
         if (k == 1) then
            gamma1_face = s% gamma1(k)
            return
         end if
         call get_face_weights(s, k, alfa, beta)
         gamma1_face = alfa*s% gamma1(k) + beta*s%gamma1(k-1)
      end function get_gamma1_face

      subroutine get_face_weights(s, k, alfa, beta)
         type (star_info), pointer :: s
         integer, intent(in) :: k
         real(dp), intent(out) :: alfa, beta
         ! face_value(k) = alfa*cell_value(k) + beta*cell_value(k-1)
         if (k == 1) call mesa_error(__FILE__,__LINE__,'bad k==1 for get_face_weights')
         alfa = s% dq(k-1)/(s% dq(k-1) + s% dq(k))
         beta = 1d0 - alfa
      end subroutine get_face_weights


      end module run_star_extras
