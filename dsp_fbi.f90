MODULE dsp_fbi

  use dsp_header

  implicit none

  real(kind=DP),    save :: k_perp, e0, dfs, alp, sigma_p
  real(kind=DP),    save :: theta, mp
  complex(kind=DP), save :: phi_jpara


! initial values !
    data phi_jpara / ( 0.0_DP, 0._DP ) /
    data k_perp  / 1.0_DP    /
    data e0      / 0.001_DP  /
    data dfs     / 2.d-5     /
    data alp     / 7.d-4     /
    data sigma_p / 5._DP     /

    data theta   / 0.0_DP    /
    data mp      / 0.5_DP    /

  private

  public fbi_init, fbi_func, fbi_set_k_perp, fbi_set_phi_jpara
  public fbi_set_theta


CONTAINS


  SUBROUTINE fbi_init


      write(olog,fmt="(a,1p,1e15.7)")
      write(olog,fmt="(a,1p,2e15.7)") "# phi_jpara = ", phi_jpara
      write(olog,fmt="(a,1p,1e15.7)") "# k_perp = ", k_perp
      write(olog,fmt="(a,1p,1e15.7)") "# e0     = ", e0
      write(olog,fmt="(a,1p,1e15.7)") "# dfs    = ", dfs
      write(olog,fmt="(a,1p,1e15.7)") "# alp    = ", alp
      write(olog,fmt="(a,1p,1e15.7)") "# sigma_p= ", sigma_p
      write(olog,fmt="(a,1p,1e15.7)") "# theta  = ", theta
      write(olog,fmt="(a,1p,1e15.7)") "# mp     = ", mp
      write(olog,fmt="(a,1p,1e15.7)")


  END SUBROUTINE fbi_init


  SUBROUTINE fbi_set_k_perp( k_perp_new )

    real(kind=DP) :: k_perp_new

      k_perp = k_perp_new

  END SUBROUTINE fbi_set_k_perp


  SUBROUTINE fbi_set_theta( theta_new )

    real(kind=DP) :: theta_new

      theta = theta_new

  END SUBROUTINE fbi_set_theta


  SUBROUTINE fbi_set_phi_jpara( phi_jpara_new )

    complex(kind=DP) :: phi_jpara_new

      phi_jpara = phi_jpara_new

  END SUBROUTINE fbi_set_phi_jpara


  SUBROUTINE fbi_func( omega )

    complex(kind=DP), intent(out) :: omega

    complex(kind=DP) :: ui = ( 0._DP, 1._DP )

!!!      omega = ( k_perp * e0 - ui * k_perp**2 * dfs ) &

!!!    write(6,*) "theta = ", theta

      omega = ( k_perp * e0 * ( mp * sin(theta) - cos(theta) ) &
              - ui * k_perp**2 * dfs ) &
            / ( 1._DP + phi_jpara * sigma_p )        &
            - ui * alp * 2._DP


  END SUBROUTINE fbi_func


END MODULE dsp_fbi
