MODULE dynadv_cen2
   !!======================================================================
   !!                       ***  MODULE  dynadv  ***
   !! Ocean dynamics: Update the momentum trend with the flux form advection
   !!                 using a 2nd order centred scheme
   !!======================================================================
   !! History :  2.0  ! 2006-08  (G. Madec, S. Theetten)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_adv_cen2  : flux form momentum advection (ln_dynadv_cen2=T) using a 2nd order centred scheme
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control

   ! USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   ! USE prtctl         ! Print control
   ! USE in_out_manager ! I/O manager
   ! USE lib_mpp        ! MPP library
   ! USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_adv_cen2   ! routine called by step.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynadv_cen2.F90 13151 2020-06-24 12:38:26Z gm $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_adv_cen2( kt, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_cen2  ***
      !!
      !! ** Purpose :   Compute the now momentum advection trend in flux form
      !!              and the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!
      !! ** Action  :   (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) updated with the now vorticity term trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt           ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs    ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv     ! ocean velocities and RHS of momentum equation
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zfu_t, zfu_f, zfu_uw, zfu, zfuf
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zfv_t, zfv_f, zfv_vw, zfv, zfvf, zfw
      REAL(wp) :: zetu, zetv, zefu, zefv
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  zuu, zvv
      REAL(wp), DIMENSION(jpi,jpj) :: z2d
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_adv_cen2 : 2nd order flux form momentum advection'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      IF( l_trddyn ) THEN           ! trends: store the input trends
         zfu_uw(:,:,:) = puu(:,:,:,Krhs)
         zfv_vw(:,:,:) = pvv(:,:,:,Krhs)
      ENDIF
      !
      !                             !==  Horizontal advection  ==!
      !
#if defined key_fluxcornermirror
      IF( kt == nit000 ) THEN
         WRITE(numout,*) 'dyn_adv_cen2 : mirror flow !'
      ENDIF
      !
      zuu(:,:,:) = puu(:,:,:,Kmm) ; zvv = pvv(:,:,:,Kmm)
      !
      z2d(:,:) = 1._wp
      DO jk = 1, jpkm1
        DO_2D_00_00
              IF ( (     tmask(ji  ,jj+1,1) + tmask(ji+1,jj+1,1)  &
                   &   + tmask(ji  ,jj  ,1) + tmask(ji+1,jj  ,1)  ) == 3._wp) THEN ! salient angle
                  !
                  IF      ( tmask(ji  ,jj  ,1)==0._wp ) THEN   ! South
                   zuu(ji,jj,jk) = - pvv(ji+1,jj  ,jk,Kmm)
                   zvv(ji,jj,jk) = - puu(ji  ,jj+1,jk,Kmm)
                   z2d(ji,jj) = 2._wp
                  ELSE IF ( tmask(ji+1,jj  ,1)==0._wp ) THEN   ! East
                    zuu(ji  ,jj,jk) = + pvv(ji,jj  ,jk,Kmm)
                    zvv(ji+1,jj,jk) = + puu(ji,jj+1,jk,Kmm)
                  ELSE IF ( tmask(ji+1,jj+1,1)==0._wp ) THEN   ! North
                    zuu(ji  ,jj+1,jk) = - pvv(ji,jj,jk,Kmm)
                    zvv(ji+1,jj  ,jk) = - puu(ji,jj,jk,Kmm)
                  ELSE IF ( tmask(ji  ,jj+1,1)==0._wp ) THEN   ! West
                    zuu(ji,jj+1,jk) = + pvv(ji+1,jj  ,jk,Kmm)
                    zvv(ji,jj  ,jk) = + puu(ji  ,jj,jk,Kmm)
                  END IF
                  !
              END IF
            !
        END_2D
      END DO
      !
      DO jk = 1, jpkm1                    ! horizontal transport
         zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u(:,:,jk,Kmm) * puu(:,:,jk,Kmm) ! puu advective (zuu advected)
         zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v(:,:,jk,Kmm) * pvv(:,:,jk,Kmm)
         !
         zfuf(:,:,jk) = 0.25_wp * e2u(:,:) * e3u(:,:,jk,Kmm) * zuu(:,:,jk) ! mirrored flux
         zfvf(:,:,jk) = 0.25_wp * e1v(:,:) * e3v(:,:,jk,Kmm) * zvv(:,:,jk)
         !
         DO_2D_10_10
            zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji+1,jj  ,jk,Kmm) )
            zfv_f(ji  ,jj  ,jk) = ( zfv(ji,jj,jk) + zfv(ji+1,jj,jk) ) * ( zuu(ji,jj,jk) + zuu(ji  ,jj+1,jk) ) * z2d(ji,jj)
            zfu_f(ji  ,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji,jj+1,jk) ) * ( zvv(ji,jj,jk) + zvv(ji+1,jj  ,jk) ) * z2d(ji,jj)
            zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji  ,jj+1,jk,Kmm) )
            !
            ! zfu_t(ji+1,jj  ,jk) = ( zfuf(ji,jj,jk) + zfuf(ji+1,jj,jk) ) * ( zuu(ji,jj,jk) + zuu(ji+1,jj  ,jk) )
            ! zfv_f(ji  ,jj  ,jk) = ( zfvf(ji,jj,jk) + zfvf(ji+1,jj,jk) ) * ( zuu(ji,jj,jk) + zuu(ji  ,jj+1,jk) )
            ! zfu_f(ji  ,jj  ,jk) = ( zfuf(ji,jj,jk) + zfuf(ji,jj+1,jk) ) * ( zvv(ji,jj,jk) + zvv(ji+1,jj  ,jk) )
            ! zfv_t(ji  ,jj+1,jk) = ( zfvf(ji,jj,jk) + zfvf(ji,jj+1,jk) ) * ( zvv(ji,jj,jk) + zvv(ji  ,jj+1,jk) )
            !
         END_2D
         !

         !
         DO_2D_00_00
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)       &
               &                                    +    zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  )    &
               &                                    * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)      &
               &                                    +    zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  )   &
               &                                    * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)
         END_2D
      END DO
#elif defined key_fluxcornerdot
      IF( kt == nit000 ) THEN
         WRITE(numout,*) 'dyn_adv_cen2 : mirror flow (dot)!'
      ENDIF
      !
      zuu(:,:,:) = puu(:,:,:,Kmm) ; zvv = pvv(:,:,:,Kmm)
      !
      DO jk = 1, jpkm1
        DO_2D_00_00
              IF ( (     tmask(ji  ,jj+1,1) + tmask(ji+1,jj+1,1)  &
                   &   + tmask(ji  ,jj  ,1) + tmask(ji+1,jj  ,1)  ) == 3._wp) THEN ! salient angle
                  !
                  IF      ( tmask(ji  ,jj  ,1)==0._wp ) THEN   ! South
                   zuu(ji,jj,jk) = puu(ji,jj+1,jk,Kmm)
                   zvv(ji,jj,jk) = pvv(ji+1,jj,jk,Kmm)
                  END IF
                  !
              END IF
            !
        END_2D
      END DO
      !
      DO jk = 1, jpkm1                    ! horizontal transport
          zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u(:,:,jk,Kmm) * puu(:,:,jk,Kmm) ! puu advective (zuu advected)
          zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v(:,:,jk,Kmm) * pvv(:,:,jk,Kmm)
          !
          zfuf(:,:,jk) = 0.25_wp * e2u(:,:) * e3u(:,:,jk,Kmm) * zuu(:,:,jk) ! mirrored flux
          zfvf(:,:,jk) = 0.25_wp * e1v(:,:) * e3v(:,:,jk,Kmm) * zvv(:,:,jk)
         DO_2D_10_10
            zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji+1,jj  ,jk,Kmm) )
            zfv_f(ji  ,jj  ,jk) = ( zfvf(ji,jj,jk) + zfvf(ji+1,jj,jk) ) * ( zuu(ji,jj,jk) + zuu(ji  ,jj+1,jk) )
            zfu_f(ji  ,jj  ,jk) = ( zfuf(ji,jj,jk) + zfuf(ji,jj+1,jk) ) * ( zvv(ji,jj,jk) + zvv(ji+1,jj  ,jk) )
            zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji  ,jj+1,jk,Kmm) )
            !
         END_2D
         !

         !
         DO_2D_00_00
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)       &
               &                                    +    zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  )    &
               &                                    * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)      &
               &                                    +    zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  )   &
               &                                    * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)
         END_2D
      END DO
#elif defined key_fluxcorner45
      IF( kt == nit000 ) THEN
         WRITE(numout,*) 'dyn_adv_cen2 : fmasked !'
      ENDIF
      !
      DO jk = 1, jpkm1                    ! horizontal transport
         zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u(:,:,jk,Kmm) * puu(:,:,jk,Kmm)
         zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v(:,:,jk,Kmm) * pvv(:,:,jk,Kmm)
         DO_2D_10_10
            !
            zfu_t(ji+1,jj  ,jk) =                   ( zfu(ji,jj,jk) + zfu(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji+1,jj  ,jk,Kmm) )
            zfv_f(ji  ,jj  ,jk) = fmask(ji,jj,jk) * ( zfv(ji,jj,jk) + zfv(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) )
            zfu_f(ji  ,jj  ,jk) = fmask(ji,jj,jk) * ( zfu(ji,jj,jk) + zfu(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) )
            zfv_t(ji  ,jj+1,jk) =                   ( zfv(ji,jj,jk) + zfv(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji  ,jj+1,jk,Kmm) )
         END_2D
         !
         DO_2D_00_00
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)       &
               &                                    +    zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  )    &
               &                                    * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)      &
               &                                    +    zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  )   &
               &                                    * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)
         END_2D
      END DO
#elif defined key_flux4corner
      IF( kt == nit000 ) THEN
         WRITE(numout,*) 'dyn_adv_cen2 : x4 in corners (corner) !'
      ENDIF
      z2d(:,:) = 1._wp
      DO jk = 1, jpkm1
        DO_2D_00_00
              IF ( (     tmask(ji  ,jj+1,1) + tmask(ji+1,jj+1,1)  &
                   &   + tmask(ji  ,jj  ,1) + tmask(ji+1,jj  ,1)  ) == 3._wp) THEN ! salient angle
                  !
                  z2d(ji,jj)=0._wp
                  !
              END IF
        END_2D
      END DO
      !
      DO jk = 1, jpkm1                    ! horizontal transport
         zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u(:,:,jk,Kmm) * puu(:,:,jk,Kmm)
         zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v(:,:,jk,Kmm) * pvv(:,:,jk,Kmm)
         DO_2D_10_10
            !
            zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji+1,jj  ,jk,Kmm) )
            zfv_f(ji  ,jj  ,jk) = ( zfv(ji,jj,jk) + zfv(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) ) * z2d(ji,jj)
            zfu_f(ji  ,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) ) * z2d(ji,jj)
            zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji  ,jj+1,jk,Kmm) )
            !
         END_2D
         !
         !
         DO_2D_00_00
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)       &
               &                                    +    zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  )    &
               &                                    * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)      &
               &                                    +    zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  )   &
               &                                    * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)
         END_2D
      END DO
#elif defined key_flux4
      IF( kt == nit000 ) THEN
         WRITE(numout,*) 'dyn_adv_cen2 : x4 in corners 45° !'
      ENDIF
      !
      DO jk = 1, jpkm1                    ! horizontal transport
         zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u(:,:,jk,Kmm) * puu(:,:,jk,Kmm)
         zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v(:,:,jk,Kmm) * pvv(:,:,jk,Kmm)
         DO_2D_10_10
            !
            zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji+1,jj  ,jk,Kmm) )
            zfv_f(ji  ,jj  ,jk) = ( zfv(ji,jj,jk) + zfv(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) ) * ( 4._wp - 3._wp*fmask(ji,jj,jk) )
            zfu_f(ji  ,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) ) * ( 4._wp - 3._wp*fmask(ji,jj,jk) )
            zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji  ,jj+1,jk,Kmm) )
            !
         END_2D
         !
         !
         DO_2D_00_00
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)       &
               &                                    +    zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  )    &
               &                                    * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)      &
               &                                    +    zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  )   &
               &                                    * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)
         END_2D
      END DO
#elif defined key_flux2
      IF( kt == nit000 ) THEN
        ! mean for speeds (not mass transport) in the corners is doubled.
         WRITE(numout,*) 'dyn_adv_cen2 : x2 in corners !'
      ENDIF
      !
      DO jk = 1, jpkm1                    ! horizontal transport
         zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u(:,:,jk,Kmm) * puu(:,:,jk,Kmm)
         zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v(:,:,jk,Kmm) * pvv(:,:,jk,Kmm)
         DO_2D_10_10
            !
            zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji+1,jj  ,jk,Kmm) )
            zfv_f(ji  ,jj  ,jk) = ( zfv(ji,jj,jk) + zfv(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) ) * ( 2._wp - fmask(ji,jj,jk) )
            zfu_f(ji  ,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) ) * ( 2._wp - fmask(ji,jj,jk) )
            zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji  ,jj+1,jk,Kmm) )
            !
         END_2D
         !
         !
         DO_2D_00_00
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)       &
               &                                    +    zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  )    &
               &                                    * r1_e1e2u(ji,jj) / e3u(ji,jj,jk,Kmm)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)      &
               &                                    +    zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  )   &
               &                                    * r1_e1e2v(ji,jj) / e3v(ji,jj,jk,Kmm)
         END_2D
      END DO
#else
      DO jk = 1, jpkm1                    ! horizontal transport
         zfu(:,:,jk) = 0.25_wp * e2u(:,:) * e3u(:,:,jk,Kmm) * puu(:,:,jk,Kmm)
         zfv(:,:,jk) = 0.25_wp * e1v(:,:) * e3v(:,:,jk,Kmm) * pvv(:,:,jk,Kmm)
         DO_2D_10_10
            zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji+1,jj  ,jk,Kmm) )
            zfv_f(ji  ,jj  ,jk) = ( zfv(ji,jj,jk) + zfv(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) )
            zfu_f(ji  ,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) )
            zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji  ,jj+1,jk,Kmm) )
         END_2D
         DO_2D_00_00
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)    &
               &                           + zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  ) * r1_e1e2u(ji,jj)   &
               &                           / e3u(ji,jj,jk,Kmm)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)    &
               &                           + zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  ) * r1_e1e2v(ji,jj)   &
               &                           / e3v(ji,jj,jk,Kmm)
         END_2D
      END DO
#endif
      IF( l_trddyn ) THEN           ! trends: send trend to trddyn for diagnostic
         zfu_uw(:,:,:) = puu(:,:,:,Krhs) - zfu_uw(:,:,:)
         zfv_vw(:,:,:) = pvv(:,:,:,Krhs) - zfv_vw(:,:,:)
         CALL trd_dyn( zfu_uw, zfv_vw, jpdyn_keg, kt, Kmm )
         zfu_t(:,:,:) = puu(:,:,:,Krhs)
         zfv_t(:,:,:) = pvv(:,:,:,Krhs)
      ENDIF
      !
      !                             !==  Vertical advection  ==!
      !
      DO_2D_00_00
         zfu_uw(ji,jj,jpk) = 0._wp   ;   zfv_vw(ji,jj,jpk) = 0._wp
         zfu_uw(ji,jj, 1 ) = 0._wp   ;   zfv_vw(ji,jj, 1 ) = 0._wp
      END_2D
      IF( ln_linssh ) THEN                ! linear free surface: advection through the surface
         DO_2D_00_00
            zfu_uw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * ww(ji,jj,1) + e1e2t(ji+1,jj) * ww(ji+1,jj,1) ) * puu(ji,jj,1,Kmm)
            zfv_vw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * ww(ji,jj,1) + e1e2t(ji,jj+1) * ww(ji,jj+1,1) ) * pvv(ji,jj,1,Kmm)
         END_2D
      ENDIF
      DO jk = 2, jpkm1                    ! interior advective fluxes
         DO_2D_01_01
            zfw(ji,jj,jk) = 0.25_wp * e1e2t(ji,jj) * ww(ji,jj,jk)
         END_2D
         DO_2D_00_00
            zfu_uw(ji,jj,jk) = ( zfw(ji,jj,jk) + zfw(ji+1,jj  ,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji,jj,jk-1,Kmm) )
            zfv_vw(ji,jj,jk) = ( zfw(ji,jj,jk) + zfw(ji  ,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji,jj,jk-1,Kmm) )
         END_2D
      END DO
      DO_3D_00_00( 1, jpkm1 )
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - ( zfu_uw(ji,jj,jk) - zfu_uw(ji,jj,jk+1) ) * r1_e1e2u(ji,jj)   &
            &                                      / e3u(ji,jj,jk,Kmm)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - ( zfv_vw(ji,jj,jk) - zfv_vw(ji,jj,jk+1) ) * r1_e1e2v(ji,jj)   &
            &                                      / e3v(ji,jj,jk,Kmm)
      END_3D
      !
      IF( l_trddyn ) THEN                 ! trends: send trend to trddyn for diagnostic
         zfu_t(:,:,:) = puu(:,:,:,Krhs) - zfu_t(:,:,:)
         zfv_t(:,:,:) = pvv(:,:,:,Krhs) - zfv_t(:,:,:)
         CALL trd_dyn( zfu_t, zfv_t, jpdyn_zad, kt, Kmm )
      ENDIF
      !                                   ! Control print
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' cen2 adv - Ua: ', mask1=umask,   &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=           ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
   END SUBROUTINE dyn_adv_cen2

   !!==============================================================================
END MODULE dynadv_cen2
