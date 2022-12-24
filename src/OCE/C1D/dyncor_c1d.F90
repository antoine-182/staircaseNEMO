MODULE dyncor_c1d
   !!======================================================================
   !!                     ***  MODULE  dyncor_c1d  ***
   !! Ocean Dynamics :   Coriolis term in 1D configuration
   !!=====================================================================
   !! History :  2.0  !  2004-09  (C. Ethe)  Original code
   !!            3.0  !  2008-04  (G. Madec)  style only
   !!----------------------------------------------------------------------
#if defined key_c1d
   !!----------------------------------------------------------------------
   !!   'key_c1d'                                          1D Configuration
   !!----------------------------------------------------------------------
   !!   cor_c1d       : Coriolis factor at T-point (1D configuration)
   !!   dyn_cor_c1d   : vorticity trend due to Coriolis at T-point
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control

   USE sbcwave        ! Surface Waves (add Stokes-Coriolis force)
   USE sbc_oce , ONLY : ln_stcor    ! use Stoke-Coriolis force
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   cor_c1d      ! called by nemogcm.F90
   PUBLIC   dyn_cor_c1d  ! called by step1d.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dyncor_c1d.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE cor_c1d
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE cor_c1d  ***
      !! 
      !! ** Purpose : set the Coriolis factor at T-point
      !!----------------------------------------------------------------------
      REAL(wp) ::   zphi0, zbeta, zf0   ! local scalars
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cor_c1d : Coriolis factor at T-point'
      IF(lwp) WRITE(numout,*) '~~~~~~~'

      !
   END SUBROUTINE cor_c1d


   SUBROUTINE dyn_cor_c1d( kt, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_cor_c1d  ***
      !! 
      !! ** Purpose :   Compute the now Coriolis trend and add it to 
      !!               the general trend of the momentum equation in 1D case.
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   kt        ! ocean time-step index
      INTEGER                             , INTENT(in   ) ::   Kmm, Krhs ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::   puu, pvv  ! ocean velocities and RHS of momentum equation
      !!
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_cor_c1d : total vorticity trend in 1D'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      ENDIF
      !
      IF( ln_stcor ) THEN
         DO_3D_00_00( 1, jpkm1 )
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + ff_t(ji,jj) * (pvv(ji,jj,jk,Kmm) + vsd(ji,jj,jk))
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - ff_t(ji,jj) * (puu(ji,jj,jk,Kmm) + usd(ji,jj,jk))
         END_3D
      ELSE
         DO_3D_00_00( 1, jpkm1 )
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + ff_t(ji,jj) * pvv(ji,jj,jk,Kmm)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - ff_t(ji,jj) * puu(ji,jj,jk,Kmm)
         END_3D
      END IF
      
      !
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' cor  - Ua: ', mask1=umask,  &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=' Va: '       , mask2=vmask )
      !
   END SUBROUTINE dyn_cor_c1d

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO 1D Configuration
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE cor_c1d              ! Empty routine
       IMPLICIT NONE
   END SUBROUTINE cor_c1d   
   SUBROUTINE dyn_cor_c1d ( kt )      ! Empty routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'dyn_cor_c1d: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_cor_c1d
#endif

   !!=====================================================================
END MODULE dyncor_c1d
