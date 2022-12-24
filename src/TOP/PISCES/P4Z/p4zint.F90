MODULE p4zint
   !!======================================================================
   !!                         ***  MODULE p4zint  ***
   !! TOP :   PISCES interpolation and computation of various accessory fields
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
   !!   p4z_int        :  interpolation and computation of various accessory fields
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_int  
   REAL(wp) ::   xksilim = 16.5e-6_wp   ! Half-saturation constant for the Si half-saturation constant computation

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zint.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_int( kt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_int  ***
      !!
      !! ** Purpose :   interpolation and computation of various accessory fields
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kbb, Kmm ! time level indices
      !
      INTEGER  :: ji, jj                 ! dummy loop indices
      REAL(wp) :: zvar                   ! local variable
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_int')
      !
      ! Computation of phyto and zoo metabolic rate
      ! -------------------------------------------
      tgfunc (:,:,:) = EXP( 0.063913 * ts(:,:,:,jp_tem,Kmm) )
      tgfunc2(:,:,:) = EXP( 0.07608  * ts(:,:,:,jp_tem,Kmm) )

      ! Computation of the silicon dependant half saturation  constant for silica uptake
      ! ---------------------------------------------------
      DO ji = 1, jpi
         DO jj = 1, jpj
            zvar = tr(ji,jj,1,jpsil,Kbb) * tr(ji,jj,1,jpsil,Kbb)
            xksimax(ji,jj) = MAX( xksimax(ji,jj), ( 1.+ 7.* zvar / ( xksilim * xksilim + zvar ) ) * 1e-6 )
         END DO
      END DO
      !
      IF( nday_year == nyear_len(1) ) THEN
         xksi   (:,:) = xksimax(:,:)
         xksimax(:,:) = 0._wp
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_int')
      !
   END SUBROUTINE p4z_int

   !!======================================================================
END MODULE p4zint
