MODULE p5zmort
   !!======================================================================
   !!                         ***  MODULE p5zmort  ***
   !! TOP :   PISCES Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p5z_mort       :   Compute the mortality terms for phytoplankton
   !!   p5z_mort_init  :   Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zlim
   USE p5zlim          !  Phytoplankton limitation terms
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p5z_mort    
   PUBLIC   p5z_mort_init    

   !! * Shared module variables
   REAL(wp), PUBLIC :: wchln    !:
   REAL(wp), PUBLIC :: wchlp   !:
   REAL(wp), PUBLIC :: wchld   !:
   REAL(wp), PUBLIC :: wchldm  !:
   REAL(wp), PUBLIC :: mpratn   !:
   REAL(wp), PUBLIC :: mpratp  !:
   REAL(wp), PUBLIC :: mpratd  !:

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p5zmort.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p5z_mort( kt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      !!---------------------------------------------------------------------

      CALL p5z_nano( Kbb, Krhs )            ! nanophytoplankton
      CALL p5z_pico( Kbb, Krhs )            ! picophytoplankton
      CALL p5z_diat( Kbb, Krhs )            ! diatoms

   END SUBROUTINE p5z_mort


   SUBROUTINE p5z_nano( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph
      REAL(wp) :: zfactfe, zfactch, zfactn, zfactp, zprcaca
      REAL(wp) :: ztortp , zrespp , zmortp
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p5z_nano')
      !
      prodcal(:,:,:) = 0.  !: calcite production variable set to zero
      DO_3D_11_11( 1, jpkm1 )
         zcompaph = MAX( ( tr(ji,jj,jk,jpphy,Kbb) - 1e-9 ), 0.e0 )
         !   Squared mortality of Phyto similar to a sedimentation term during
         !   blooms (Doney et al. 1996)
         !   -----------------------------------------------------------------
         zrespp = wchln * 1.e6 * xstep * xdiss(ji,jj,jk) * zcompaph * tr(ji,jj,jk,jpphy,Kbb)

         !   Phytoplankton linear mortality
         !   ------------------------------
         ztortp = mpratn * xstep  * zcompaph
         zmortp = zrespp + ztortp

         !   Update the arrays TRA which contains the biological sources and sinks

         zfactn  = tr(ji,jj,jk,jpnph,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         zfactp  = tr(ji,jj,jk,jppph,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         zfactfe = tr(ji,jj,jk,jpnfe,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         zfactch = tr(ji,jj,jk,jpnch,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         tr(ji,jj,jk,jpphy,Krhs) = tr(ji,jj,jk,jpphy,Krhs) - zmortp
         tr(ji,jj,jk,jpnph,Krhs) = tr(ji,jj,jk,jpnph,Krhs) - zmortp * zfactn
         tr(ji,jj,jk,jppph,Krhs) = tr(ji,jj,jk,jppph,Krhs) - zmortp * zfactp
         tr(ji,jj,jk,jpnch,Krhs) = tr(ji,jj,jk,jpnch,Krhs) - zmortp * zfactch
         tr(ji,jj,jk,jpnfe,Krhs) = tr(ji,jj,jk,jpnfe,Krhs) - zmortp * zfactfe
         zprcaca = xfracal(ji,jj,jk) * zmortp
         !
         prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         !
         tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) - zprcaca
         tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) - 2. * zprcaca
         tr(ji,jj,jk,jpcal,Krhs) = tr(ji,jj,jk,jpcal,Krhs) + zprcaca
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + zmortp
         tr(ji,jj,jk,jppon,Krhs) = tr(ji,jj,jk,jppon,Krhs) + zmortp * zfactn
         tr(ji,jj,jk,jppop,Krhs) = tr(ji,jj,jk,jppop,Krhs) + zmortp * zfactp
         prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortp
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + zmortp * zfactfe
      END_3D
      !
       IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p5z_nano')
      !
   END SUBROUTINE p5z_nano


   SUBROUTINE p5z_pico( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_pico  ***
      !!
      !! ** Purpose :   Compute the mortality terms for picophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph
      REAL(wp) :: zfactfe, zfactch, zfactn, zfactp
      REAL(wp) :: ztortp , zrespp , zmortp 
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p5z_pico')
      !
      DO_3D_11_11( 1, jpkm1 )
         zcompaph = MAX( ( tr(ji,jj,jk,jppic,Kbb) - 1e-9 ), 0.e0 )
         !  Squared mortality of Phyto similar to a sedimentation term during
         !  blooms (Doney et al. 1996)
         !  -----------------------------------------------------------------
         zrespp = wchlp * 1.e6 * xstep * xdiss(ji,jj,jk) * zcompaph * tr(ji,jj,jk,jppic,Kbb)

         !     Phytoplankton mortality 
         ztortp = mpratp * xstep  * zcompaph
         zmortp = zrespp + ztortp

         !   Update the arrays TRA which contains the biological sources and sinks

         zfactn = tr(ji,jj,jk,jpnpi,Kbb)/(tr(ji,jj,jk,jppic,Kbb)+rtrn)
         zfactp = tr(ji,jj,jk,jpppi,Kbb)/(tr(ji,jj,jk,jppic,Kbb)+rtrn)
         zfactfe = tr(ji,jj,jk,jppfe,Kbb)/(tr(ji,jj,jk,jppic,Kbb)+rtrn)
         zfactch = tr(ji,jj,jk,jppch,Kbb)/(tr(ji,jj,jk,jppic,Kbb)+rtrn)
         tr(ji,jj,jk,jppic,Krhs) = tr(ji,jj,jk,jppic,Krhs) - zmortp
         tr(ji,jj,jk,jpnpi,Krhs) = tr(ji,jj,jk,jpnpi,Krhs) - zmortp * zfactn
         tr(ji,jj,jk,jpppi,Krhs) = tr(ji,jj,jk,jpppi,Krhs) - zmortp * zfactp
         tr(ji,jj,jk,jppch,Krhs) = tr(ji,jj,jk,jppch,Krhs) - zmortp * zfactch
         tr(ji,jj,jk,jppfe,Krhs) = tr(ji,jj,jk,jppfe,Krhs) - zmortp * zfactfe
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + zmortp
         tr(ji,jj,jk,jppon,Krhs) = tr(ji,jj,jk,jppon,Krhs) + zmortp * zfactn
         tr(ji,jj,jk,jppop,Krhs) = tr(ji,jj,jk,jppop,Krhs) + zmortp * zfactp
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + zmortp * zfactfe
         prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortp
      END_3D
      !
       IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('pico')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p5z_pico')
      !
   END SUBROUTINE p5z_pico


   SUBROUTINE p5z_diat( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_diat  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diatoms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zfactfe,zfactsi,zfactch, zfactn, zfactp, zcompadi
      REAL(wp) ::  zrespp2, ztortp2, zmortp2
      REAL(wp) ::  zlim2, zlim1
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p5z_diat')
      !

      DO_3D_11_11( 1, jpkm1 )

         zcompadi = MAX( ( tr(ji,jj,jk,jpdia,Kbb) - 1E-9), 0. )

         !   Aggregation term for diatoms is increased in case of nutrient
         !   stress as observed in reality. The stressed cells become more
         !   sticky and coagulate to sink quickly out of the euphotic zone
         !   -------------------------------------------------------------
         !  Phytoplankton squared mortality
         !  -------------------------------
         zlim2   = xlimdia(ji,jj,jk) * xlimdia(ji,jj,jk)
         zlim1   = 0.25 * ( 1. - zlim2 ) / ( 0.25 + zlim2 ) 
         zrespp2 = 1.e6 * xstep * (  wchld + wchldm * zlim1 ) * xdiss(ji,jj,jk) * zcompadi * tr(ji,jj,jk,jpdia,Kbb)

         !  Phytoplankton linear mortality 
         !  ------------------------------
         ztortp2 = mpratd * xstep  * zcompadi
         zmortp2 = zrespp2 + ztortp2

         !   Update the arrays tr(:,:,:,:,Krhs) which contains the biological sources and sinks
         !   ---------------------------------------------------------------------
         zfactn  = tr(ji,jj,jk,jpndi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         zfactp  = tr(ji,jj,jk,jppdi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         zfactch = tr(ji,jj,jk,jpdch,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         zfactfe = tr(ji,jj,jk,jpdfe,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         zfactsi = tr(ji,jj,jk,jpdsi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         tr(ji,jj,jk,jpdia,Krhs) = tr(ji,jj,jk,jpdia,Krhs) - zmortp2 
         tr(ji,jj,jk,jpndi,Krhs) = tr(ji,jj,jk,jpndi,Krhs) - zmortp2 * zfactn
         tr(ji,jj,jk,jppdi,Krhs) = tr(ji,jj,jk,jppdi,Krhs) - zmortp2 * zfactp
         tr(ji,jj,jk,jpdch,Krhs) = tr(ji,jj,jk,jpdch,Krhs) - zmortp2 * zfactch
         tr(ji,jj,jk,jpdfe,Krhs) = tr(ji,jj,jk,jpdfe,Krhs) - zmortp2 * zfactfe
         tr(ji,jj,jk,jpdsi,Krhs) = tr(ji,jj,jk,jpdsi,Krhs) - zmortp2 * zfactsi
         tr(ji,jj,jk,jpgsi,Krhs) = tr(ji,jj,jk,jpgsi,Krhs) + zmortp2 * zfactsi
         tr(ji,jj,jk,jpgoc,Krhs) = tr(ji,jj,jk,jpgoc,Krhs) + zrespp2 
         tr(ji,jj,jk,jpgon,Krhs) = tr(ji,jj,jk,jpgon,Krhs) + zrespp2 * zfactn
         tr(ji,jj,jk,jpgop,Krhs) = tr(ji,jj,jk,jpgop,Krhs) + zrespp2 * zfactp
         tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zrespp2 * zfactfe
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + ztortp2
         tr(ji,jj,jk,jppon,Krhs) = tr(ji,jj,jk,jppon,Krhs) + ztortp2 * zfactn
         tr(ji,jj,jk,jppop,Krhs) = tr(ji,jj,jk,jppop,Krhs) + ztortp2 * zfactp
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + ztortp2 * zfactfe
         prodpoc(ji,jj,jk)   = prodpoc(ji,jj,jk) + ztortp2
         prodgoc(ji,jj,jk)   = prodgoc(ji,jj,jk) + zrespp2
      END_3D
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p5z_diat')
      !
   END SUBROUTINE p5z_diat


   SUBROUTINE p5z_mort_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p5z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the nampismort namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampismort
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namp5zmort/ wchln, wchlp, wchld, wchldm, mpratn, mpratp, mpratd
      !!----------------------------------------------------------------------

      READ  ( numnatp_ref, namp5zmort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namp5zmort in reference namelist' )

      READ  ( numnatp_cfg, namp5zmort, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namp5zmort in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp5zmort )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mortality, namp5zmort'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    quadratic mortality of phytoplankton      wchln     =', wchln
         WRITE(numout,*) '    quadratic mortality of picophyto.         wchlp     =', wchlp
         WRITE(numout,*) '    quadratic mortality of diatoms            wchld     =', wchld
         WRITE(numout,*) '    Additional quadratic mortality of diatoms wchldm    =', wchldm
         WRITE(numout,*) '    nanophyto. mortality rate                 mpratn    =', mpratn
         WRITE(numout,*) '    picophyto. mortality rate                 mpratp    =', mpratp
         WRITE(numout,*) '    Diatoms mortality rate                    mpratd    =', mpratd
      ENDIF

   END SUBROUTINE p5z_mort_init

   !!======================================================================
END MODULE p5zmort
