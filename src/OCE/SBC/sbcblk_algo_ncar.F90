MODULE sbcblk_algo_ncar
   !!======================================================================
   !!                   ***  MODULE  sbcblk_algo_ncar  ***
   !! Computes:
   !!   * bulk transfer coefficients C_D, C_E and C_H
   !!   * air temp. and spec. hum. adjusted from zt (2m) to zu (10m) if needed
   !!   * the effective bulk wind speed at 10m U_blk
   !!   => all these are used in bulk formulas in sbcblk.F90
   !!
   !!    Using the bulk formulation/param. of Large & Yeager 2008
   !!
   !!       Routine turb_ncar maintained and developed in AeroBulk
   !!                     (https://github.com/brodeau/aerobulk/)
   !!
   !!                         L. Brodeau, 2015
   !!=====================================================================
   !! History :  3.6  !  2016-02  (L.Brodeau) successor of old turb_ncar of former sbcblk_core.F90
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   turb_ncar  : computes the bulk turbulent transfer coefficients
   !!                   adjusts t_air and q_air from zt to zu m
   !!                   returns the effective bulk wind speed at 10m
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE sbcwave, ONLY   :  cdn_wave ! wave module
#if defined key_si3 || defined key_cice
   USE sbc_ice         ! Surface boundary condition: ice fields
#endif
   !
   USE iom             ! I/O manager library
   USE lib_mpp         ! distribued memory computing library
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE lib_fortran     ! to use key_nosignedzero

   USE sbcblk_phy      ! all thermodynamics functions, rho_air, q_sat, etc... !LB

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: TURB_NCAR   ! called by sbcblk.F90

   INTEGER , PARAMETER ::   nb_itt = 5        ! number of itterations
   !! * Substitutions
#  include "do_loop_substitute.h90"

   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE turb_ncar( zt, zu, sst, t_zt, ssq, q_zt, U_zu, &
      &                  Cd, Ch, Ce, t_zu, q_zu, U_blk,      &
      &                  Cdn, Chn, Cen                       )
      !!----------------------------------------------------------------------------------
      !!                      ***  ROUTINE  turb_ncar  ***
      !!
      !! ** Purpose :   Computes turbulent transfert coefficients of surface
      !!                fluxes according to Large & Yeager (2004) and Large & Yeager (2008)
      !!                If relevant (zt /= zu), adjust temperature and humidity from height zt to zu
      !!                Returns the effective bulk wind speed at 10m to be used in the bulk formulas
      !!
      !!
      !! INPUT :
      !! -------
      !!    *  zt   : height for temperature and spec. hum. of air            [m]
      !!    *  zu   : height for wind speed (usually 10m)                     [m]
      !!    *  sst  : bulk SST                                                [K]
      !!    *  t_zt : potential air temperature at zt                         [K]
      !!    *  ssq  : specific humidity at saturation at SST                  [kg/kg]
      !!    *  q_zt : specific humidity of air at zt                          [kg/kg]
      !!    *  U_zu : scalar wind speed at zu                                 [m/s]
      !!
      !!
      !! OUTPUT :
      !! --------
      !!    *  Cd     : drag coefficient
      !!    *  Ch     : sensible heat coefficient
      !!    *  Ce     : evaporation coefficient
      !!    *  t_zu   : pot. air temperature adjusted at wind height zu       [K]
      !!    *  q_zu   : specific humidity of air        //                    [kg/kg]
      !!    *  U_blk  : bulk wind speed at zu                                 [m/s]
      !!
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in   )                     ::   zt       ! height for t_zt and q_zt                    [m]
      REAL(wp), INTENT(in   )                     ::   zu       ! height for U_zu                             [m]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   sst      ! sea surface temperature                [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   t_zt     ! potential air temperature              [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   ssq      ! sea surface specific humidity           [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   q_zt     ! specific air humidity at zt             [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(jpi,jpj) ::   U_zu     ! relative wind module at zu                [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cd       ! transfer coefficient for momentum         (tau)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ch       ! transfer coefficient for sensible heat (Q_sens)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Ce       ! transfert coefficient for evaporation   (Q_lat)
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   t_zu     ! pot. air temp. adjusted at zu               [K]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   q_zu     ! spec. humidity adjusted at zu           [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   U_blk    ! bulk wind speed at zu                     [m/s]
      REAL(wp), INTENT(  out), DIMENSION(jpi,jpj) ::   Cdn, Chn, Cen ! neutral transfer coefficients
      !
      INTEGER :: j_itt
      LOGICAL :: l_zt_equal_zu = .FALSE.      ! if q and t are given at same height as U
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   Cx_n10        ! 10m neutral latent/sensible coefficient
      REAL(wp), DIMENSION(jpi,jpj) ::   sqrt_Cd_n10   ! root square of Cd_n10
      REAL(wp), DIMENSION(jpi,jpj) ::   zeta_u        ! stability parameter at height zu
      REAL(wp), DIMENSION(jpi,jpj) ::   zpsi_h_u
      REAL(wp), DIMENSION(jpi,jpj) ::   ztmp0, ztmp1, ztmp2
      REAL(wp), DIMENSION(jpi,jpj) ::   stab          ! stability test integer
      !!----------------------------------------------------------------------------------
      l_zt_equal_zu = ( ABS(zu - zt) < 0.01_wp ) ! testing "zu == zt" is risky with double precision

      U_blk = MAX( 0.5_wp , U_zu )   !  relative wind speed at zu (normally 10m), we don't want to fall under 0.5 m/s

      !! First guess of stability:
      ztmp0 = virt_temp(t_zt, q_zt) - virt_temp(sst, ssq) ! air-sea difference of virtual pot. temp. at zt
      stab  = 0.5_wp + sign(0.5_wp,ztmp0)                           ! stab = 1 if dTv > 0  => STABLE, 0 if unstable

      !! Neutral coefficients at 10m:
      IF( ln_cdgw ) THEN      ! wave drag case
         cdn_wave(:,:) = cdn_wave(:,:) + rsmall * ( 1._wp - tmask(:,:,1) )
         ztmp0   (:,:) = cdn_wave(:,:)
      ELSE
      ztmp0 = cd_neutral_10m( U_blk )
      ENDIF

      sqrt_Cd_n10 = SQRT( ztmp0 )

      !! Initializing transf. coeff. with their first guess neutral equivalents :
      Cd = ztmp0
      Ce = 1.e-3_wp*( 34.6_wp * sqrt_Cd_n10 )
      Ch = 1.e-3_wp*sqrt_Cd_n10*(18._wp*stab + 32.7_wp*(1._wp - stab))
      stab = sqrt_Cd_n10   ! Temporaty array !!! stab == SQRT(Cd)
 
      IF( ln_cdgw ) THEN
	Cen = Ce
	Chn = Ch
      ENDIF

      !! First guess of temperature and humidity at height zu:
      t_zu = MAX( t_zt ,  180._wp )   ! who knows what's given on masked-continental regions...
      q_zu = MAX( q_zt , 1.e-6_wp )   !               "

      !! ITERATION BLOCK
      DO j_itt = 1, nb_itt
         !
         ztmp1 = t_zu - sst   ! Updating air/sea differences
         ztmp2 = q_zu - ssq

         ! Updating turbulent scales :   (L&Y 2004 eq. (7))
         ztmp0 = stab*U_blk       ! u*       (stab == SQRT(Cd))
         ztmp1 = Ch/stab*ztmp1    ! theta*   (stab == SQRT(Cd))
         ztmp2 = Ce/stab*ztmp2    ! q*       (stab == SQRT(Cd))

         ! Estimate the inverse of Monin-Obukov length (1/L) at height zu:
         ztmp0 = One_on_L( t_zu, q_zu, ztmp0, ztmp1, ztmp2 )
         
         !! Stability parameters :
         zeta_u   = zu*ztmp0
         zeta_u = sign( min(abs(zeta_u),10._wp), zeta_u )
         zpsi_h_u = psi_h( zeta_u )

         !! Shifting temperature and humidity at zu (L&Y 2004 eq. (9b-9c))
         IF( .NOT. l_zt_equal_zu ) THEN
            !! Array 'stab' is free for the moment so using it to store 'zeta_t'
            stab = zt*ztmp0
            stab = SIGN( MIN(ABS(stab),10._wp), stab )  ! Temporaty array stab == zeta_t !!!
            stab = LOG(zt/zu) + zpsi_h_u - psi_h(stab)                   ! stab just used as temp array again!
            t_zu = t_zt - ztmp1/vkarmn*stab    ! ztmp1 is still theta*  L&Y 2004 eq.(9b)
            q_zu = q_zt - ztmp2/vkarmn*stab    ! ztmp2 is still q*      L&Y 2004 eq.(9c)
            q_zu = max(0._wp, q_zu)
         ENDIF

         ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 eq. 9a)...
         !   In very rare low-wind conditions, the old way of estimating the
         !   neutral wind speed at 10m leads to a negative value that causes the code
         !   to crash. To prevent this a threshold of 0.25m/s is imposed.
         ztmp2 = psi_m(zeta_u)
         IF( ln_cdgw ) THEN      ! surface wave case
            stab = vkarmn / ( vkarmn / sqrt_Cd_n10 - ztmp2 )  ! (stab == SQRT(Cd))
            Cd   = stab * stab
            ztmp0 = (LOG(zu/10._wp) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
            ztmp2 = stab / sqrt_Cd_n10   ! (stab == SQRT(Cd))
            ztmp1 = 1._wp + Chn * ztmp0     
            Ch    = Chn * ztmp2 / ztmp1  ! L&Y 2004 eq. (10b)
            ztmp1 = 1._wp + Cen * ztmp0
            Ce    = Cen * ztmp2 / ztmp1  ! L&Y 2004 eq. (10c)

         ELSE
         ! Update neutral wind speed at 10m and neutral Cd at 10m (L&Y 2004 eq. 9a)...
         !   In very rare low-wind conditions, the old way of estimating the
         !   neutral wind speed at 10m leads to a negative value that causes the code
         !   to crash. To prevent this a threshold of 0.25m/s is imposed.
         ztmp0 = MAX( 0.25_wp , U_blk/(1._wp + sqrt_Cd_n10/vkarmn*(LOG(zu/10._wp) - ztmp2)) ) ! U_n10 (ztmp2 == psi_m(zeta_u))
         ztmp0 = cd_neutral_10m(ztmp0)                                               ! Cd_n10
         Cdn(:,:) = ztmp0
         sqrt_Cd_n10 = sqrt(ztmp0)

         stab    = 0.5_wp + sign(0.5_wp,zeta_u)                        ! update stability
         Cx_n10  = 1.e-3_wp*sqrt_Cd_n10*(18._wp*stab + 32.7_wp*(1._wp - stab))  ! L&Y 2004 eq. (6c-6d)    (Cx_n10 == Ch_n10)
         Chn(:,:) = Cx_n10

         !! Update of transfer coefficients:
         ztmp1 = 1._wp + sqrt_Cd_n10/vkarmn*(LOG(zu/10._wp) - ztmp2)   ! L&Y 2004 eq. (10a) (ztmp2 == psi_m(zeta_u))
         Cd      = ztmp0 / ( ztmp1*ztmp1 )
         stab = SQRT( Cd ) ! Temporary array !!! (stab == SQRT(Cd))

         ztmp0 = (LOG(zu/10._wp) - zpsi_h_u) / vkarmn / sqrt_Cd_n10
         ztmp2 = stab / sqrt_Cd_n10   ! (stab == SQRT(Cd))
         ztmp1 = 1._wp + Cx_n10*ztmp0    ! (Cx_n10 == Ch_n10)
         Ch  = Cx_n10*ztmp2 / ztmp1   ! L&Y 2004 eq. (10b)

         Cx_n10  = 1.e-3_wp * (34.6_wp * sqrt_Cd_n10)  ! L&Y 2004 eq. (6b)    ! Cx_n10 == Ce_n10
         Cen(:,:) = Cx_n10
         ztmp1 = 1._wp + Cx_n10*ztmp0
         Ce  = Cx_n10*ztmp2 / ztmp1  ! L&Y 2004 eq. (10c)
         ENDIF

      END DO !DO j_itt = 1, nb_itt

   END SUBROUTINE turb_ncar


   FUNCTION cd_neutral_10m( pw10 )
      !!----------------------------------------------------------------------------------      
      !! Estimate of the neutral drag coefficient at 10m as a function
      !! of neutral wind  speed at 10m
      !!
      !! Origin: Large & Yeager 2008 eq.(11a) and eq.(11b)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pw10           ! scalar wind speed at 10m (m/s)
      REAL(wp), DIMENSION(jpi,jpj)             :: cd_neutral_10m
      !
      INTEGER  ::     ji, jj     ! dummy loop indices
      REAL(wp) :: zgt33, zw, zw6 ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO_2D_11_11
         !
         zw  = pw10(ji,jj)
         zw6 = zw*zw*zw
         zw6 = zw6*zw6
         !
         ! When wind speed > 33 m/s => Cyclone conditions => special treatment
         zgt33 = 0.5_wp + SIGN( 0.5_wp, (zw - 33._wp) )   ! If pw10 < 33. => 0, else => 1
         !
         cd_neutral_10m(ji,jj) = 1.e-3_wp * ( &
            &       (1._wp - zgt33)*( 2.7_wp/zw + 0.142_wp + zw/13.09_wp - 3.14807E-10_wp*zw6) & ! wind <  33 m/s
            &      +    zgt33   *      2.34_wp )                                                 ! wind >= 33 m/s
         !
         cd_neutral_10m(ji,jj) = MAX(cd_neutral_10m(ji,jj), 1.E-6_wp)
         !
      END_2D
      !
   END FUNCTION cd_neutral_10m


   FUNCTION psi_m( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for momentum
      !!    !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_m
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) :: zx2, zx, zstab   ! local scalars
      !!----------------------------------------------------------------------------------
      DO_2D_11_11
         zx2 = SQRT( ABS( 1._wp - 16._wp*pzeta(ji,jj) ) )
         zx2 = MAX( zx2 , 1._wp )
         zx  = SQRT( zx2 )
         zstab = 0.5_wp + SIGN( 0.5_wp , pzeta(ji,jj) )
         !
         psi_m(ji,jj) =        zstab  * (-5._wp*pzeta(ji,jj))       &          ! Stable
            &          + (1._wp - zstab) * (2._wp*LOG((1._wp + zx)*0.5_wp)   &          ! Unstable
            &               + LOG((1._wp + zx2)*0.5_wp) - 2._wp*ATAN(zx) + rpi*0.5_wp)  !    "
         !
      END_2D
   END FUNCTION psi_m


   FUNCTION psi_h( pzeta )
      !!----------------------------------------------------------------------------------
      !! Universal profile stability function for temperature and humidity
      !!    !! Psis, L&Y 2004 eq. (8c), (8d), (8e)
      !!
      !! pzeta : stability paramenter, z/L where z is altitude measurement
      !!         and L is M-O length
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: psi_h
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pzeta
      !
      INTEGER  ::   ji, jj     ! dummy loop indices
      REAL(wp) :: zx2, zstab  ! local scalars
      !!----------------------------------------------------------------------------------
      !
      DO_2D_11_11
         zx2 = SQRT( ABS( 1._wp - 16._wp*pzeta(ji,jj) ) )
         zx2 = MAX( zx2 , 1._wp )
         zstab = 0.5_wp + SIGN( 0.5_wp , pzeta(ji,jj) )
         !
         psi_h(ji,jj) =         zstab  * (-5._wp*pzeta(ji,jj))        &  ! Stable
            &           + (1._wp - zstab) * (2._wp*LOG( (1._wp + zx2)*0.5_wp ))   ! Unstable
         !
      END_2D
   END FUNCTION psi_h

   !!======================================================================
END MODULE sbcblk_algo_ncar
