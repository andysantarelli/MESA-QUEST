
! photo utilities for saving variables
      subroutine quasi_star_photo_write(id, iounit)
         integer, intent(in) :: id, iounit
         write(iounit) M_BH_new, L_BH, &
                      R_B, R_B_raw, &
                      M_dot, log_dm_dt, rad_eff, &
                      opacity, M_cav, P_rad, &
                      P_gas, nabla_ad, M_dot_Bondi, &
                      L_Bondi, c_s_smooth, &
                      c_s, R_sc
      end subroutine quasi_star_photo_write


      subroutine quasi_star_photo_read(id, iounit, ierr)
         integer, intent(in) :: id, iounit
         integer, intent(out) :: ierr
         ierr = 0
         read(iounit, iostat=ierr) M_BH_new, L_BH, &
                      R_B, R_B_raw, &
                      M_dot, log_dm_dt, rad_eff, &
                      opacity, M_cav, P_rad, &
                      P_gas, nabla_ad, M_dot_Bondi, &
                      L_Bondi, c_s_smooth, c_s, R_sc
      end subroutine quasi_star_photo_read



