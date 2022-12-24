MODULE p4zbio
   !!======================================================================
   !!                         ***  MODULE p4zbio  ***
   !! TOP :   PISCES bio-model
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
   !!   p4z_bio        :   computes the interactions between the different
   !!                      compartments of PISCES
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zmort         !  Mortality terms for phytoplankton
   USE p4zmicro        !  Sources and sinks of microzooplankton
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p5zlim          !  Co-limitations of differents nutrients
   USE p5zprod         !  Growth rate of the 2 phyto groups
   USE p5zmort         !  Mortality terms for phytoplankton
   USE p5zmicro        !  Sources and sinks of microzooplankton
   USE p5zmeso         !  Sources and sinks of mesozooplankton
   USE p4zrem          !  Remineralisation of organic matter
   USE p4zpoc          !  Remineralization of organic particles
   USE p4zagg          !  Aggregation of particles
   USE p4zfechem
   USE p4zligand       !  Prognostic ligand model
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
  
   IMPLICIT NONE
   PRIVATE

   PUBLIC  p4z_bio    

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zbio.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_bio ( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_bio  ***
      !!
      !! ** Purpose :   Ecosystem model in the whole ocean: computes the
      !!              different interactions between the different compartments
      !!              of PISCES
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER, INTENT(in) :: Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER             :: ji, jj, jk, jn
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_bio')
      !
      !     ASSIGN THE SHEAR RATE THAT IS USED FOR AGGREGATION
      !     OF PHYTOPLANKTON AND DETRITUS

      xdiss(:,:,:) = 1.
!!gm the use of nmld should be better here?
      DO_3D_11_11( 2, jpkm1 )
!!gm  :  use nmln  and test on jk ...  less memory acces
         IF( gdepw(ji,jj,jk+1,Kmm) > hmld(ji,jj) )   xdiss(ji,jj,jk) = 0.01
      END_3D

      CALL p4z_opt     ( kt, knt, Kbb, Kmm       )     ! Optic: PAR in the water column
      CALL p4z_sink    ( kt, knt, Kbb, Kmm, Krhs )     ! vertical flux of particulate organic matter
      CALL p4z_fechem  ( kt, knt, Kbb, Kmm, Krhs )     ! Iron chemistry/scavenging
      !
      IF( ln_p4z ) THEN
         CALL p4z_lim  ( kt, knt, Kbb, Kmm       )     ! co-limitations by the various nutrients
         CALL p4z_prod ( kt, knt, Kbb, Kmm, Krhs )     ! phytoplankton growth rate over the global ocean. 
         !                                          ! (for each element : C, Si, Fe, Chl )
         CALL p4z_mort ( kt,      Kbb,      Krhs )     ! phytoplankton mortality
         !                                          ! zooplankton sources/sinks routines 
         CALL p4z_micro( kt, knt, Kbb,      Krhs )     ! microzooplankton
         CALL p4z_meso ( kt, knt, Kbb,      Krhs )     ! mesozooplankton
      ELSE
         CALL p5z_lim  ( kt, knt, Kbb, Kmm       )     ! co-limitations by the various nutrients
         CALL p5z_prod ( kt, knt, Kbb, Kmm, Krhs )     ! phytoplankton growth rate over the global ocean. 
         !                                          ! (for each element : C, Si, Fe, Chl )
         CALL p5z_mort ( kt,      Kbb,      Krhs      )     ! phytoplankton mortality
         !                                          ! zooplankton sources/sinks routines 
         CALL p5z_micro( kt, knt, Kbb,      Krhs )           ! microzooplankton
         CALL p5z_meso ( kt, knt, Kbb,      Krhs )           ! mesozooplankton
      ENDIF
      !
      CALL p4z_agg     ( kt, knt, Kbb,      Krhs )     ! Aggregation of particles
      CALL p4z_rem     ( kt, knt, Kbb, Kmm, Krhs )     ! remineralization terms of organic matter+scavenging of Fe
      CALL p4z_poc     ( kt, knt, Kbb, Kmm, Krhs )     ! Remineralization of organic particles
      !
      IF( ln_ligand )  &
      & CALL p4z_ligand( kt, knt, Kbb,      Krhs )
      !                                                             !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bio ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_bio')
      !
   END SUBROUTINE p4z_bio

   !!======================================================================
END MODULE p4zbio
