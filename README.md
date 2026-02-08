# MESA Quasi-star Evolutionary Simulation Toolkit (MESA-QUEST)
This is a MESA model simulating a quasi-star -- a massive, star-like object containing a black hole in its core that acts as a heavy SMBH seeding mechanism. 

We include each version here:
 - An early [version](https://github.com/andysantarelli/MESA-QUEST/blob/main/qs_Bondi_RNAAS2025/), created to replicate the work of [Ball 2012](https://ui.adsabs.harvard.edu/abs/2012PhDT.........1B/abstract). Corresponding publication found [here](https://iopscience.iop.org/article/10.3847/2515-5172/adef33/meta). MESA version r24.03.1.
 - A version containing a new implementation based on [Coughlin 2024](https://iopscience.iop.org/article/10.3847/1538-4357/ad5723). Publication in press. MESA version r24.08.1.
 - A [version](https://github.com/andysantarelli/MESA-QUEST/blob/main/qs_LRD_ApJL2026/) using the previous implementation along with the MESA Colors module. Corresponding publication found [here](https://iopscience.iop.org/article/10.3847/2041-8213/ae3713). MESA version r25.12.1.

Additional features, processes, data, and plots will be added upon publication of their respective papers.

If you use MESA-QUEST in a publication, please cite the corresponding publications above, and follow the [MESA best practices](https://docs.mesastar.org/en/latest/using_mesa/best_practices.html) guidelines. 

The primary changes have been made within the [MESA](https://mesa.sourceforge.net/) [`run_star_extras.f90`](https://github.com/andysantarelli/MESA-QUEST/blob/main/qs_LRD_ApJL2026/src/run_star_extras.f90) file. The core changes are found in the subroutine below. 

```fortran
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
          G       = s% cgrav(s% nz)
          c_s     = s% csound(s% nz)
          rho     = s%rho(s% nz)
          opacity = s% opacity(s% nz)
          gamma1  = s%gamma1(s% nz)
          R       = s% r(s% nz) 
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
          log_dm_dt   = safe_log10(dm) - safe_log10(dt)
          s% xtra(5)  = log_dm_dt
          s% xtra(6)  = rad_eff
          s% xtra(7)  = opacity
          s% xtra(8)  = M_cav
          s% xtra(9)  = M_dot_Bondi
          s% xtra(10) = L_Bondi
          s% xtra(11) = c_s
          s% xtra(12) = R_sc

          ! Apply core modifications
          if (startup .and. .not. restart) then
             call star_write_model(id, "initial_model.mod", ierr)

             s% M_center = new_core_mass !M_BH
             s% L_center = L_BH

            if (s% x_logical_ctrl(2)) then
                call do1_relax_R_center(s, R_sc, ierr)
            else
                call do1_relax_R_center(s, R_B, ierr)  
            end if

          else
           
              s% M_center = new_core_mass * Msun
              s% mstar = s% mstar - rad_eff * M_dot * dt
              s% xmstar = s% mstar - s% M_center
              s% L_center = L_BH

            if (s% x_logical_ctrl(2)) then
                call do1_relax_R_center(s, R_sc, ierr)
            else
                call do1_relax_R_center(s, R_B, ierr)  ! Use smoothed R_B
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

```

We include relevant adjustable parameters at the bottom of the [`inlist_project`](https://github.com/andysantarelli/MESA-QUEST/blob/main/qs_LRD_ApJL2026/inlist_project) file. These will be expanded to include additional options and scenarios as we develop these models further. 

```fortran
    x_logical_ctrl(1) = .true.    ! accrete onto black hole
    x_logical_ctrl(2) = .true.    ! use saturated-convection inner boundary
    x_logical_ctrl(3) = .true.    ! use M_star Eddington rate ! .false. means use Ball et al. 2012
    x_logical_ctrl(4) = .false.   ! use super-Eddington winds
    
    x_ctrl(1) = 0.1               ! radiative efficiency, epsilon
    x_ctrl(2) = 0.1               ! convective efficiency, eta (only used for Ball accretion)
    x_ctrl(3) = 1.0               ! BH mass growth timestep factor
    x_ctrl(4) = 1.0               ! accretion factor, alpha
    x_ctrl(5) = 1d5               ! profile interval in Msun units
    x_ctrl(6) = -1                ! black hole mass limit in Msun units (negative value means no limit)
```


![qwest? qwest!](https://github.com/andysantarelli/MESA-QUEST/blob/main/qwest%EF%BC%9F/qwest!.png)
 - Art/character credit: [Matthew J. Wills](https://swordscomic.com/)
