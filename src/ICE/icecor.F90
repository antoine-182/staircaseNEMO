MODULE icecor
   !!======================================================================
   !!                     ***  MODULE  icecor  ***
   !!   sea-ice: Corrections on sea-ice variables at the end of the time step
   !!======================================================================
   !! History :  3.0  !  2006-04  (M. Vancoppenolle) Original code
   !!            3.5  !  2014-06  (C. Rousset)       Complete rewriting/cleaning
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!    ice_cor      : corrections on sea-ice variables
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean domain
   USE phycst         ! physical constants
   USE ice            ! sea-ice: variable
   USE ice1D          ! sea-ice: thermodynamic variables
   USE iceitd         ! sea-ice: rebining
   USE icevar         ! sea-ice: operations
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! fortran utilities (glob_sum + no signed zero)
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_cor   ! called by icestp.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 4.0 , NEMO Consortium (2018)
   !! $Id: icecor.F90 12489 2020-02-28 15:55:11Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_cor( kt, kn )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE ice_cor  ***
      !!               
      !! ** Purpose :   Computes corrections on sea-ice global variables at 
      !!              the end of the dynamics (kn=1) and thermodynamics (kn=2)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! number of iteration
      INTEGER, INTENT(in) ::   kn    ! 1 = after dyn ; 2 = after thermo
      !
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   zsal, zzc
      REAL(wp), DIMENSION(jpi,jpj) ::   zafx   ! concentration trends diag
      !!----------------------------------------------------------------------
      ! controls
      IF( ln_timing    )   CALL timing_start('icecor')                                                             ! timing
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icecor', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (0, 'icecor',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation
      !
      IF( kt == nit000 .AND. lwp .AND. kn == 2 ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_cor:  correct sea ice variables if out of bounds ' 
         WRITE(numout,*) '~~~~~~~'
      ENDIF
      !                             !-----------------------------------------------------
      !                             !  ice thickness must exceed himin (for temp. diff.) !
      !                             !-----------------------------------------------------
      WHERE( a_i(:,:,:) >= epsi20 )   ;   h_i(:,:,:) = v_i(:,:,:) / a_i(:,:,:)
      ELSEWHERE                       ;   h_i(:,:,:) = 0._wp
      END WHERE
      WHERE( h_i(:,:,:) < rn_himin )      a_i(:,:,:) = a_i(:,:,:) * h_i(:,:,:) / rn_himin
      !
      !                             !-----------------------------------------------------
      !                             !  ice concentration should not exceed amax          !
      !                             !-----------------------------------------------------
      at_i(:,:) = SUM( a_i(:,:,:), dim=3 )
      DO jl = 1, jpl
         WHERE( at_i(:,:) > rn_amax_2d(:,:) )   a_i(:,:,jl) = a_i(:,:,jl) * rn_amax_2d(:,:) / at_i(:,:)
      END DO
    
      !                             !-----------------------------------------------------
      IF ( nn_icesal == 2 ) THEN    !  salinity must stay in bounds [Simin,Simax]        !
         !                          !-----------------------------------------------------
         zzc = rhoi * r1_Dt_ice
         DO jl = 1, jpl
            DO_2D_11_11
               zsal = sv_i(ji,jj,jl)
               sv_i(ji,jj,jl) = MIN(  MAX( rn_simin*v_i(ji,jj,jl) , sv_i(ji,jj,jl) ) , rn_simax*v_i(ji,jj,jl)  )
               sfx_res(ji,jj) = sfx_res(ji,jj) - ( sv_i(ji,jj,jl) - zsal ) * zzc   ! associated salt flux
            END_2D
         END DO
      ENDIF
      !                             !-----------------------------------------------------
      !                             !  Rebin categories with thickness out of bounds     !
      !                             !-----------------------------------------------------
      IF ( jpl > 1 )   CALL ice_itd_reb( kt )

      !                             !-----------------------------------------------------
      CALL ice_var_zapsmall         !  Zap small values                                  !
      !                             !-----------------------------------------------------

      !                             !-----------------------------------------------------
      IF( kn == 2 ) THEN            !  Ice drift case: Corrections to avoid wrong values !
         DO_2D_00_00
            IF ( at_i(ji,jj) == 0._wp ) THEN    ! what to do if there is no ice
               IF ( at_i(ji+1,jj) == 0._wp )   u_ice(ji  ,jj) = 0._wp   ! right side
               IF ( at_i(ji-1,jj) == 0._wp )   u_ice(ji-1,jj) = 0._wp   ! left side
               IF ( at_i(ji,jj+1) == 0._wp )   v_ice(ji,jj  ) = 0._wp   ! upper side
               IF ( at_i(ji,jj-1) == 0._wp )   v_ice(ji,jj-1) = 0._wp   ! bottom side
            ENDIF
         END_2D
         CALL lbc_lnk_multi( 'icecor', u_ice, 'U', -1., v_ice, 'V', -1. )
      ENDIF

      !                             !-----------------------------------------------------
      SELECT CASE( kn )             !  Diagnostics                                       !
      !                             !-----------------------------------------------------
      CASE( 1 )                        !--- dyn trend diagnostics
         !
         IF( ln_icediachk .OR. iom_use('hfxdhc') ) THEN
            diag_heat(:,:) = - SUM(SUM( e_i (:,:,1:nlay_i,:) - e_i_b (:,:,1:nlay_i,:), dim=4 ), dim=3 ) * r1_Dt_ice &      ! W.m-2
               &             - SUM(SUM( e_s (:,:,1:nlay_s,:) - e_s_b (:,:,1:nlay_s,:), dim=4 ), dim=3 ) * r1_Dt_ice
            diag_sice(:,:) =   SUM(     sv_i(:,:,:)          - sv_i_b(:,:,:)                  , dim=3 ) * r1_Dt_ice * rhoi
            diag_vice(:,:) =   SUM(     v_i (:,:,:)          - v_i_b (:,:,:)                  , dim=3 ) * r1_Dt_ice * rhoi
            diag_vsnw(:,:) =   SUM(     v_s (:,:,:)          - v_s_b (:,:,:)                  , dim=3 ) * r1_Dt_ice * rhos
         ENDIF
         !                       ! concentration tendency (dynamics)
         IF( iom_use('afxdyn') .OR. iom_use('afxthd') .OR. iom_use('afxtot') ) THEN 
            zafx(:,:) = SUM( a_i(:,:,:) - a_i_b(:,:,:), dim=3 ) * r1_Dt_ice 
            CALL iom_put( 'afxdyn' , zafx )
         ENDIF
         !
      CASE( 2 )                        !--- thermo trend diagnostics & ice aging
         !
         oa_i(:,:,:) = oa_i(:,:,:) + a_i(:,:,:) * rDt_ice   ! ice natural aging incrementation
         !
         IF( ln_icediachk .OR. iom_use('hfxdhc') ) THEN
            diag_heat(:,:) = diag_heat(:,:) &
               &             - SUM(SUM( e_i (:,:,1:nlay_i,:) - e_i_b (:,:,1:nlay_i,:), dim=4 ), dim=3 ) * r1_Dt_ice &
               &             - SUM(SUM( e_s (:,:,1:nlay_s,:) - e_s_b (:,:,1:nlay_s,:), dim=4 ), dim=3 ) * r1_Dt_ice
            diag_sice(:,:) = diag_sice(:,:) &
               &             + SUM(     sv_i(:,:,:)          - sv_i_b(:,:,:)                  , dim=3 ) * r1_Dt_ice * rhoi
            diag_vice(:,:) = diag_vice(:,:) &
               &             + SUM(     v_i (:,:,:)          - v_i_b (:,:,:)                  , dim=3 ) * r1_Dt_ice * rhoi
            diag_vsnw(:,:) = diag_vsnw(:,:) &
               &             + SUM(     v_s (:,:,:)          - v_s_b (:,:,:)                  , dim=3 ) * r1_Dt_ice * rhos
            CALL iom_put ( 'hfxdhc' , diag_heat ) 
         ENDIF
         !                       ! concentration tendency (total + thermo)
         IF( iom_use('afxdyn') .OR. iom_use('afxthd') .OR. iom_use('afxtot') ) THEN 
            zafx(:,:) = zafx(:,:) + SUM( a_i(:,:,:) - a_i_b(:,:,:), dim=3 ) * r1_Dt_ice
            CALL iom_put( 'afxthd' , SUM( a_i(:,:,:) - a_i_b(:,:,:), dim=3 ) * r1_Dt_ice )
            CALL iom_put( 'afxtot' , zafx )
         ENDIF
         !
      END SELECT
      !
      ! controls
      IF( sn_cfctl%l_prtctl ) &
         &                 CALL ice_prt3D   ('icecor')                                                             ! prints
      IF( ln_icectl .AND. kn == 2 ) &
         &                 CALL ice_prt     ( kt, iiceprt, jiceprt, 2, ' - Final state - ' )                       ! prints
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icecor', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'icecor',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation
      IF( ln_timing    )   CALL timing_stop ('icecor')                                                             ! timing
      !
   END SUBROUTINE ice_cor

#else
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icecor
