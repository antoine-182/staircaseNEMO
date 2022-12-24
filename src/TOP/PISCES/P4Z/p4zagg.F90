MODULE p4zagg
   !!======================================================================
   !!                         ***  MODULE p4zagg  ***
   !! TOP :  PISCES  aggregation of particles
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   p4z_agg       :  Compute aggregation of particles
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_agg         ! called in p4zbio.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zagg.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_agg ( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_agg  ***
      !!
      !! ** Purpose :   Compute aggregation of particles
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      INTEGER, INTENT(in) ::   Kbb, Krhs ! time level indices
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zagg, zagg1, zagg2, zagg3, zagg4
      REAL(wp) ::   zaggpoc1, zaggpoc2, zaggpoc3, zaggpoc4
      REAL(wp) ::   zaggpoc , zaggfe, zaggdoc, zaggdoc2, zaggdoc3
      REAL(wp) ::   zaggpon , zaggdon, zaggdon2, zaggdon3
      REAL(wp) ::   zaggpop, zaggdop, zaggdop2, zaggdop3
      REAL(wp) ::   zaggtmp, zfact, zmax
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_agg')
      !
      !  Exchange between organic matter compartments due to coagulation/disaggregation
      !  ---------------------------------------------------
      IF( ln_p4z ) THEN
         !
         DO_3D_11_11( 1, jpkm1 )
            !
            zfact = xstep * xdiss(ji,jj,jk)
            !  Part I : Coagulation dependent on turbulence
            zagg1 = 25.9  * zfact * tr(ji,jj,jk,jppoc,Kbb) * tr(ji,jj,jk,jppoc,Kbb)
            zagg2 = 4452. * zfact * tr(ji,jj,jk,jppoc,Kbb) * tr(ji,jj,jk,jpgoc,Kbb)

            ! Part II : Differential settling

            !  Aggregation of small into large particles
            zagg3 =  47.1 * xstep * tr(ji,jj,jk,jppoc,Kbb) * tr(ji,jj,jk,jpgoc,Kbb)
            zagg4 =  3.3  * xstep * tr(ji,jj,jk,jppoc,Kbb) * tr(ji,jj,jk,jppoc,Kbb)

            zagg   = zagg1 + zagg2 + zagg3 + zagg4
            zaggfe = zagg * tr(ji,jj,jk,jpsfe,Kbb) / ( tr(ji,jj,jk,jppoc,Kbb) + rtrn )

            ! Aggregation of DOC to POC : 
            ! 1st term is shear aggregation of DOC-DOC
            ! 2nd term is shear aggregation of DOC-POC
            ! 3rd term is differential settling of DOC-POC
            zaggdoc  = ( ( 0.369 * 0.3 * tr(ji,jj,jk,jpdoc,Kbb) + 102.4 * tr(ji,jj,jk,jppoc,Kbb) ) * zfact       &
            &            + 2.4 * xstep * tr(ji,jj,jk,jppoc,Kbb) ) * 0.3 * tr(ji,jj,jk,jpdoc,Kbb)
            ! transfer of DOC to GOC : 
            ! 1st term is shear aggregation
            ! 2nd term is differential settling 
            zaggdoc2 = ( 3.53E3 * zfact + 0.1 * xstep ) * tr(ji,jj,jk,jpgoc,Kbb) * 0.3 * tr(ji,jj,jk,jpdoc,Kbb)
            ! tranfer of DOC to POC due to brownian motion
            zaggdoc3 =  114. * 0.3 * tr(ji,jj,jk,jpdoc,Kbb) *xstep * 0.3 * tr(ji,jj,jk,jpdoc,Kbb)

            !  Update the trends
            tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) - zagg + zaggdoc + zaggdoc3
            tr(ji,jj,jk,jpgoc,Krhs) = tr(ji,jj,jk,jpgoc,Krhs) + zagg + zaggdoc2
            tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) - zaggfe
            tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zaggfe
            tr(ji,jj,jk,jpdoc,Krhs) = tr(ji,jj,jk,jpdoc,Krhs) - zaggdoc - zaggdoc2 - zaggdoc3
            !
            conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zagg + zaggdoc + zaggdoc3
            prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zagg + zaggdoc2
            !
         END_3D
      ELSE    ! ln_p5z
        !
         DO_3D_11_11( 1, jpkm1 )
            !
            zfact = xstep * xdiss(ji,jj,jk)
            !  Part I : Coagulation dependent on turbulence
            zaggtmp = 25.9  * zfact * tr(ji,jj,jk,jppoc,Kbb)
            zaggpoc1 = zaggtmp * tr(ji,jj,jk,jppoc,Kbb)
            zaggtmp = 4452. * zfact * tr(ji,jj,jk,jpgoc,Kbb)
            zaggpoc2 = zaggtmp * tr(ji,jj,jk,jppoc,Kbb)

            ! Part II : Differential settling

            !  Aggregation of small into large particles
            zaggtmp =  47.1 * xstep * tr(ji,jj,jk,jpgoc,Kbb)
            zaggpoc3 = zaggtmp * tr(ji,jj,jk,jppoc,Kbb)
            zaggtmp =  3.3  * xstep * tr(ji,jj,jk,jppoc,Kbb)
            zaggpoc4 = zaggtmp * tr(ji,jj,jk,jppoc,Kbb)

            zaggpoc   = zaggpoc1 + zaggpoc2 + zaggpoc3 + zaggpoc4
            zaggpon = zaggpoc * tr(ji,jj,jk,jppon,Kbb) / ( tr(ji,jj,jk,jppoc,Kbb) + rtrn)
            zaggpop = zaggpoc * tr(ji,jj,jk,jppop,Kbb) / ( tr(ji,jj,jk,jppoc,Kbb) + rtrn)
            zaggfe = zaggpoc * tr(ji,jj,jk,jpsfe,Kbb) / ( tr(ji,jj,jk,jppoc,Kbb)  + rtrn )

            ! Aggregation of DOC to POC : 
            ! 1st term is shear aggregation of DOC-DOC
            ! 2nd term is shear aggregation of DOC-POC
            ! 3rd term is differential settling of DOC-POC
            zaggtmp = ( ( 0.369 * 0.3 * tr(ji,jj,jk,jpdoc,Kbb) + 102.4 * tr(ji,jj,jk,jppoc,Kbb) ) * zfact       &
            &            + 2.4 * xstep * tr(ji,jj,jk,jppoc,Kbb) )
            zaggdoc  = zaggtmp * 0.3 * tr(ji,jj,jk,jpdoc,Kbb)
            zaggdon  = zaggtmp * 0.3 * tr(ji,jj,jk,jpdon,Kbb)
            zaggdop  = zaggtmp * 0.3 * tr(ji,jj,jk,jpdop,Kbb)

            ! transfer of DOC to GOC : 
            ! 1st term is shear aggregation
            ! 2nd term is differential settling 
            zaggtmp = ( 3.53E3 * zfact + 0.1 * xstep ) * tr(ji,jj,jk,jpgoc,Kbb)
            zaggdoc2 = zaggtmp * 0.3 * tr(ji,jj,jk,jpdoc,Kbb)
            zaggdon2 = zaggtmp * 0.3 * tr(ji,jj,jk,jpdon,Kbb)
            zaggdop2 = zaggtmp * 0.3 * tr(ji,jj,jk,jpdop,Kbb)

            ! tranfer of DOC to POC due to brownian motion
            zaggtmp = ( 114. * 0.3 * tr(ji,jj,jk,jpdoc,Kbb) ) * xstep
            zaggdoc3 =  zaggtmp * 0.3 * tr(ji,jj,jk,jpdoc,Kbb)
            zaggdon3 =  zaggtmp * 0.3 * tr(ji,jj,jk,jpdon,Kbb)
            zaggdop3 =  zaggtmp * 0.3 * tr(ji,jj,jk,jpdop,Kbb)

            !  Update the trends
            tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) - zaggpoc + zaggdoc + zaggdoc3
            tr(ji,jj,jk,jppon,Krhs) = tr(ji,jj,jk,jppon,Krhs) - zaggpon + zaggdon + zaggdon3
            tr(ji,jj,jk,jppop,Krhs) = tr(ji,jj,jk,jppop,Krhs) - zaggpop + zaggdop + zaggdop3
            tr(ji,jj,jk,jpgoc,Krhs) = tr(ji,jj,jk,jpgoc,Krhs) + zaggpoc + zaggdoc2
            tr(ji,jj,jk,jpgon,Krhs) = tr(ji,jj,jk,jpgon,Krhs) + zaggpon + zaggdon2
            tr(ji,jj,jk,jpgop,Krhs) = tr(ji,jj,jk,jpgop,Krhs) + zaggpop + zaggdop2
            tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) - zaggfe
            tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zaggfe
            tr(ji,jj,jk,jpdoc,Krhs) = tr(ji,jj,jk,jpdoc,Krhs) - zaggdoc - zaggdoc2 - zaggdoc3
            tr(ji,jj,jk,jpdon,Krhs) = tr(ji,jj,jk,jpdon,Krhs) - zaggdon - zaggdon2 - zaggdon3
            tr(ji,jj,jk,jpdop,Krhs) = tr(ji,jj,jk,jpdop,Krhs) - zaggdop - zaggdop2 - zaggdop3
            !
            conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zaggpoc + zaggdoc + zaggdoc3
            prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zaggpoc + zaggdoc2
            !
         END_3D
         !
      ENDIF
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('agg')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_agg')
      !
   END SUBROUTINE p4z_agg

   !!======================================================================
END MODULE p4zagg
