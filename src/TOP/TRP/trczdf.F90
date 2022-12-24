MODULE trczdf
   !!==============================================================================
   !!                 ***  MODULE  trczdf  ***
   !! Ocean Passive tracers : vertical diffusive trends 
   !!=====================================================================
   !! History :  9.0  ! 2005-11  (G. Madec)  Original code
   !!       NEMO 3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA 
   !!            4.0  ! 2017-04  (G. Madec)  remove the explicit case
   !!----------------------------------------------------------------------
#if defined key_top
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_zdf      : update the tracer trend with the vertical diffusion
   !!----------------------------------------------------------------------
   USE trc           ! ocean passive tracers variables
   USE oce_trc       ! ocean dynamics and active tracers
   USE trd_oce       ! trends: ocean variables
   USE trazdf        ! tracer: vertical diffusion
!!gm do we really need this ?
   USE trcldf        ! passive tracers: lateral diffusion
!!gm
   USE trdtra        ! trends manager: tracers 
   USE prtctl_trc    ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_zdf         ! called by step.F90 
   
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trczdf.F90 12489 2020-02-28 15:55:11Z davestorkey $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_zdf( kt, Kbb, Kmm, Krhs, ptr, Kaa )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_zdf  ***
      !!
      !! ** Purpose :   compute the vertical ocean tracer physics using
      !!              an implicit time-stepping scheme.
      !!---------------------------------------------------------------------
      INTEGER                                   , INTENT(in   ) ::   kt                   ! ocean time-step index
      INTEGER                                   , INTENT(in   ) ::   Kbb, Kmm, Krhs, Kaa  ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra,jpt), INTENT(inout) ::   ptr                  ! passive tracers and RHS of tracer equation
      !
      INTEGER               ::  jk, jn
      CHARACTER (len=22)    :: charout
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra) ::   ztrtrd   ! 4D workspace
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_zdf')
      !
      IF( l_trdtrc )   ztrtrd(:,:,:,:)  = ptr(:,:,:,:,Krhs)
      !
      CALL tra_zdf_imp( kt, nittrc000, 'TRC', rDt_trc, Kbb, Kmm, Krhs, ptr, Kaa, jptra )    !   implicit scheme          
      !
      IF( l_trdtrc )   THEN                      ! save the vertical diffusive trends for further diagnostics
         DO jn = 1, jptra
            DO jk = 1, jpkm1
               ztrtrd(:,:,jk,jn) = ( ( ptr(:,:,jk,jn,Kaa) - ptr(:,:,jk,jn,Kbb) ) / rDt_trc ) - ztrtrd(:,:,jk,jn)
            END DO
            CALL trd_tra( kt, Kmm, Krhs, 'TRC', jn, jptra_zdf, ztrtrd(:,:,:,jn) )
         END DO
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF( sn_cfctl%l_prttrc )   THEN
         WRITE(charout, FMT="('zdf ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc( tab4d=tr(:,:,:,:,Kaa), mask=tmask, clinfo=ctrcnm, clinfo2='trd' )
      END IF
      !
      IF( ln_timing )  CALL timing_stop('trc_zdf')
      !
   END SUBROUTINE trc_zdf
   
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_zdf( kt )
      INTEGER, INTENT(in) :: kt  
      WRITE(*,*) 'trc_zdf: You should not have seen this print! error?', kt
   END SUBROUTINE trc_zdf
#endif
   !!==============================================================================
END MODULE trczdf
