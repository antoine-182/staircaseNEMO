MODULE dynldf_lap_blp
   !!======================================================================
   !!                   ***  MODULE  dynldf_lap_blp  ***
   !! Ocean dynamics:  lateral viscosity trend (laplacian and bilaplacian)
   !!======================================================================
   !! History : 3.7  ! 2014-01  (G. Madec, S. Masson)  Original code, re-entrant laplacian
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf_lap   : update the momentum trend with the lateral viscosity using an iso-level   laplacian operator
   !!   dyn_ldf_blp   : update the momentum trend with the lateral viscosity using an iso-level bilaplacian operator
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE ldfslp         ! iso-neutral slopes
   USE zdf_oce        ! ocean vertical physics
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_ldf_lap  ! called by dynldf.F90
   PUBLIC dyn_ldf_blp  ! called by dynldf.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynldf_lap_blp.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf_lap( kt, Kbb, Kmm, pu, pv, pu_rhs, pv_rhs, kpass )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_lap  ***
      !!
      !! ** Purpose :   Compute the before horizontal momentum diffusive
      !!      trend and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The Laplacian operator apply on horizontal velocity is
      !!      writen as :   grad_h( ahmt div_h(U )) - curl_h( ahmf curl_z(U) )
      !!
      !! ** Action : - pu_rhs, pv_rhs increased by the harmonic operator applied on pu, pv.
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kbb, Kmm   ! ocean time level indices
      INTEGER                         , INTENT(in   ) ::   kpass      ! =1/2 first or second passage
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pu, pv     ! before velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs   ! velocity trend   [m/s2]
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zsign        ! local scalars
      REAL(wp) ::   zua, zva     ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zcur, zdiv
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ldf : iso-level harmonic (laplacian) operator, pass=', kpass
         WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      IF( kpass == 1 ) THEN   ;   zsign =  1._wp      ! bilaplacian operator require a minus sign
      ELSE                    ;   zsign = -1._wp      !  (eddy viscosity coef. >0)
      ENDIF
      !
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO_2D_01_01
            !                                      ! ahm * e3 * curl  (computed from 1 to jpim1/jpjm1)
!!gm open question here : e3f  at before or now ?    probably now...
!!gm note that ahmf has already been multiplied by fmask
            zcur(ji-1,jj-1) =  &
               &      ahmf(ji-1,jj-1,jk) * e3f(ji-1,jj-1,jk) * r1_e1e2f(ji-1,jj-1)      &
               &  * (  e2v(ji  ,jj-1) * pv(ji  ,jj-1,jk) - e2v(ji-1,jj-1) * pv(ji-1,jj-1,jk)  &
               &     - e1u(ji-1,jj  ) * pu(ji-1,jj  ,jk) + e1u(ji-1,jj-1) * pu(ji-1,jj-1,jk)  )
            !                                      ! ahm * div        (computed from 2 to jpi/jpj)
!!gm note that ahmt has already been multiplied by tmask
            zdiv(ji,jj)     =   &
               &   ahmt(ji,jj,jk) * r1_e1e2t(ji,jj) / e3t(ji,jj,jk,Kbb)      &
               &     * (  e2u(ji,jj)*e3u(ji,jj,jk,Kbb) * pu(ji,jj,jk)        &
               &        - e2u(ji-1,jj)*e3u(ji-1,jj,jk,Kbb) * pu(ji-1,jj,jk)  &
               &        + e1v(ji,jj)*e3v(ji,jj,jk,Kbb) * pv(ji,jj,jk)        &
               &        - e1v(ji,jj-1)*e3v(ji,jj-1,jk,Kbb) * pv(ji,jj-1,jk)  )
         END_2D
         !
         DO_2D_00_00
            pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign * (                             &
               &              - ( zcur(ji  ,jj) - zcur(ji,jj-1) ) * r1_e2u(ji,jj)   &
               &              / e3u(ji,jj,jk,Kmm)   &
               &              + ( zdiv(ji+1,jj) - zdiv(ji,jj  ) ) * r1_e1u(ji,jj)        )
               !
            pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign * (                              &
               &                ( zcur(ji,jj  ) - zcur(ji-1,jj) ) * r1_e1v(ji,jj)   &
               &              / e3v(ji,jj,jk,Kmm)   &
               &              + ( zdiv(ji,jj+1) - zdiv(ji  ,jj) ) * r1_e2v(ji,jj)       )
         END_2D
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      !
   END SUBROUTINE dyn_ldf_lap


   SUBROUTINE dyn_ldf_blp( kt, Kbb, Kmm, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_ldf_blp  ***
      !!
      !! ** Purpose :   Compute the before lateral momentum viscous trend
      !!              and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The lateral viscous trends is provided by a bilaplacian
      !!      operator applied to before field (forward in time).
      !!      It is computed by two successive calls to dyn_ldf_lap routine
      !!
      !! ** Action :   pt(:,:,:,:,Krhs)   updated with the before rotated bilaplacian diffusion
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kbb, Kmm   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pu, pv     ! before velocity fields
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs   ! momentum trend
      !
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zulap, zvlap   ! laplacian at u- and v-point
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf_blp : bilaplacian operator momentum '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      zulap(:,:,:) = 0._wp
      zvlap(:,:,:) = 0._wp
      !
      CALL dyn_ldf_lap( kt, Kbb, Kmm, pu, pv, zulap, zvlap, 1 )   ! rotated laplacian applied to pt (output in zlap,Kbb)
      !
      CALL lbc_lnk_multi( 'dynldf_lap_blp', zulap, 'U', -1., zvlap, 'V', -1. )             ! Lateral boundary conditions
      !
      CALL dyn_ldf_lap( kt, Kbb, Kmm, zulap, zvlap, pu_rhs, pv_rhs, 2 )   ! rotated laplacian applied to zlap (output in pt(:,:,:,:,Krhs))
      !
   END SUBROUTINE dyn_ldf_blp

   !!======================================================================
END MODULE dynldf_lap_blp
