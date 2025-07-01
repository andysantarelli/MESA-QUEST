# quasi-star-model
This is a MESA model simulating a quasi-star -- a massive, star-like object containing a black hole in its core that acts as a heavy SMBH seeding mechanism. 
This is an early version, created to replicate the work of [Ball 2012](https://ui.adsabs.harvard.edu/abs/2012PhDT.........1B/abstract). 
Additional features, processes, data, and plots will be added as development and research continues. 

The primary changes have been made within the [MESA](https://mesa.sourceforge.net/) [`run_star_extras.f90`](https://github.com/andysantarelli/quasi-star-model/template/src/run_star_extras.f90) file. These were inspired by [earlier work](https://github.com/earlbellinger/black-hole-sun) done with [Earl Bellinger](https://earlbellinger.com). The core changes are found in the subroutine below. 

```fortran
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
          R_B = R_B_smooth
          
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
```

We also implement a smoothing function for the sound speed (and thus the Bondi radius) in order to decrease noise very near the black hole that can be detrimental to the structure of the quasi-star. We include the option to add weights for spatial and temporal smoothing, however the temporal smoothing is off by default. 

```fortran
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
```
