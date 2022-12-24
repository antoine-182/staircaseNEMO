MODULE diamlr
   !!======================================================================
   !!                       ***  MODULE  diamlr  ***
   !! Management of the IOM context for multiple-linear-regression analysis
   !!======================================================================
   !! History :       !  2019  (S. Mueller)
   !!----------------------------------------------------------------------

   USE par_oce        , ONLY :   wp, jpi, jpj
   USE phycst         , ONLY :   rpi
   USE in_out_manager , ONLY :   lwp, numout, ln_timing
   USE iom            , ONLY :   iom_put, iom_use, iom_update_file_name
   USE dom_oce        , ONLY :   adatrj
   USE timing         , ONLY :   timing_start, timing_stop
#if defined key_iomput
   USE xios
#endif
   USE tide_mod

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PUBLIC ::   lk_diamlr = .FALSE.

   PUBLIC ::   dia_mlr_init, dia_mlr_iom_init, dia_mlr
!!an to make it work on jeanzay
   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2019)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE dia_mlr_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_mlr_init  ***
      !!
      !! ** Purpose : initialisation of IOM context management for 
      !!              multiple-linear-regression analysis
      !!
      !!----------------------------------------------------------------------

      lk_diamlr = .TRUE.

      IF(lwp) THEN
         WRITE(numout, *)
         WRITE(numout, *) 'dia_mlr_init : initialisation of IOM context management for'
         WRITE(numout, *) '~~~~~~~~~~~~   multiple-linear-regression analysis'
      END IF

   END SUBROUTINE dia_mlr_init

   SUBROUTINE dia_mlr_iom_init
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE dia_mlr_iom_init  ***
      !!
      !! ** Purpose : IOM context setup for multiple-linear-regression
      !!              analysis
      !!
      !!----------------------------------------------------------------------
#if defined key_iom
#else
      IF( .FALSE. ) write(numout,*) 'dia_mlr_iom_init: should not see this'    ! useless statement to avoid compiler warnings
#endif

   END SUBROUTINE dia_mlr_iom_init

   SUBROUTINE dia_mlr
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dia_mlr  ***
      !!
      !! ** Purpose : update time used in multiple-linear-regression analysis
      !!
      !!----------------------------------------------------------------------

      REAL(wp), DIMENSION(jpi,jpj) ::   zadatrj2d

      IF( ln_timing )   CALL timing_start('dia_mlr')

      ! Update time to the continuous time since the start of the model run
      ! (value of adatrj converted to time in units of seconds)
      !
      ! A 2-dimensional field of constant value is sent, and subsequently used
      ! directly or transformed to a scalar or a constant 3-dimensional field as
      ! required.
      zadatrj2d(:,:) = adatrj*86400.0_wp
      IF ( iom_use('diamlr_time') ) CALL iom_put('diamlr_time', zadatrj2d)
      
      IF( ln_timing )   CALL timing_stop('dia_mlr')

   END SUBROUTINE dia_mlr

END MODULE diamlr
