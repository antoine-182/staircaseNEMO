MODULE p4zmort
   !!======================================================================
   !!                         ***  MODULE p4zmort  ***
   !! TOP :   PISCES Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
   !!   p4z_mort       : Compute the mortality terms for phytoplankton
   !!   p4z_mort_init  : Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zprod         ! Primary productivity 
   USE p4zlim          ! Phytoplankton limitation terms
   USE prtctl_trc      ! print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_mort    
   PUBLIC   p4z_mort_init    

   REAL(wp), PUBLIC ::   wchl     !:
   REAL(wp), PUBLIC ::   wchld    !:
   REAL(wp), PUBLIC ::   wchldm   !:
   REAL(wp), PUBLIC ::   mprat    !:
   REAL(wp), PUBLIC ::   mprat2   !:

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmort.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_mort( kt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to initialize and compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      !!---------------------------------------------------------------------
      !
      CALL p4z_nano( Kbb, Krhs )            ! nanophytoplankton
      !
      CALL p4z_diat( Kbb, Krhs )            ! diatoms
      !
   END SUBROUTINE p4z_mort


   SUBROUTINE p4z_nano( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zsizerat, zcompaph
      REAL(wp) ::   zfactfe, zfactch, zprcaca, zfracal
      REAL(wp) ::   ztortp , zrespp , zmortp 
      CHARACTER (len=25) ::   charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_nano')
      !
      prodcal(:,:,:) = 0._wp   ! calcite production variable set to zero
      DO_3D_11_11( 1, jpkm1 )
         zcompaph = MAX( ( tr(ji,jj,jk,jpphy,Kbb) - 1e-8 ), 0.e0 )
         !     When highly limited by macronutrients, very small cells 
         !     dominate the community. As a consequence, aggregation
         !     due to turbulence is negligible. Mortality is also set
         !     to 0
         zsizerat = MIN(1., MAX( 0., (quotan(ji,jj,jk) - 0.2) / 0.3) ) * tr(ji,jj,jk,jpphy,Kbb)
         !     Squared mortality of Phyto similar to a sedimentation term during
         !     blooms (Doney et al. 1996)
         zrespp = wchl * 1.e6 * xstep * xdiss(ji,jj,jk) * zcompaph * zsizerat 

         !     Phytoplankton mortality. This mortality loss is slightly
         !     increased when nutrients are limiting phytoplankton growth
         !     as observed for instance in case of iron limitation.
         ztortp = mprat * xstep * zcompaph / ( xkmort + tr(ji,jj,jk,jpphy,Kbb) ) * zsizerat

         zmortp = zrespp + ztortp

         !   Update the arrays TRA which contains the biological sources and sinks

         zfactfe = tr(ji,jj,jk,jpnfe,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         zfactch = tr(ji,jj,jk,jpnch,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         tr(ji,jj,jk,jpphy,Krhs) = tr(ji,jj,jk,jpphy,Krhs) - zmortp
         tr(ji,jj,jk,jpnch,Krhs) = tr(ji,jj,jk,jpnch,Krhs) - zmortp * zfactch
         tr(ji,jj,jk,jpnfe,Krhs) = tr(ji,jj,jk,jpnfe,Krhs) - zmortp * zfactfe
         zprcaca = xfracal(ji,jj,jk) * zmortp
         !
         prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         !
         zfracal = 0.5 * xfracal(ji,jj,jk)
         tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) - zprcaca
         tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) - 2. * zprcaca
         tr(ji,jj,jk,jpcal,Krhs) = tr(ji,jj,jk,jpcal,Krhs) + zprcaca
         tr(ji,jj,jk,jpgoc,Krhs) = tr(ji,jj,jk,jpgoc,Krhs) + zfracal * zmortp
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + ( 1. - zfracal ) * zmortp
         prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + ( 1. - zfracal ) * zmortp
         prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zfracal * zmortp
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + ( 1. - zfracal ) * zmortp * zfactfe
         tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zfracal * zmortp * zfactfe
      END_3D
      !
       IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_nano')
      !
   END SUBROUTINE p4z_nano


   SUBROUTINE p4z_diat( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_diat  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diatoms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zfactfe,zfactsi,zfactch, zcompadi
      REAL(wp) ::   zrespp2, ztortp2, zmortp2
      REAL(wp) ::   zlim2, zlim1
      CHARACTER (len=25) ::   charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_diat')
      !
      !    Aggregation term for diatoms is increased in case of nutrient
      !    stress as observed in reality. The stressed cells become more
      !    sticky and coagulate to sink quickly out of the euphotic zone
      !     ------------------------------------------------------------

      DO_3D_11_11( 1, jpkm1 )

         zcompadi = MAX( ( tr(ji,jj,jk,jpdia,Kbb) - 1e-9), 0. )

         !    Aggregation term for diatoms is increased in case of nutrient
         !    stress as observed in reality. The stressed cells become more
         !    sticky and coagulate to sink quickly out of the euphotic zone
         !     ------------------------------------------------------------
         !  Phytoplankton respiration 
         !     ------------------------
         zlim2   = xlimdia(ji,jj,jk) * xlimdia(ji,jj,jk)
         zlim1   = 0.25 * ( 1. - zlim2 ) / ( 0.25 + zlim2 ) 
         zrespp2 = 1.e6 * xstep * (  wchld + wchldm * zlim1 ) * xdiss(ji,jj,jk) * zcompadi * tr(ji,jj,jk,jpdia,Kbb)

         !     Phytoplankton mortality. 
         !     ------------------------
         ztortp2 = mprat2 * xstep * tr(ji,jj,jk,jpdia,Kbb)  / ( xkmort + tr(ji,jj,jk,jpdia,Kbb) ) * zcompadi 

         zmortp2 = zrespp2 + ztortp2

         !   Update the arrays tr(:,:,:,:,Krhs) which contains the biological sources and sinks
         !   ---------------------------------------------------------------------
         zfactch = tr(ji,jj,jk,jpdch,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         zfactfe = tr(ji,jj,jk,jpdfe,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         zfactsi = tr(ji,jj,jk,jpdsi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         tr(ji,jj,jk,jpdia,Krhs) = tr(ji,jj,jk,jpdia,Krhs) - zmortp2 
         tr(ji,jj,jk,jpdch,Krhs) = tr(ji,jj,jk,jpdch,Krhs) - zmortp2 * zfactch
         tr(ji,jj,jk,jpdfe,Krhs) = tr(ji,jj,jk,jpdfe,Krhs) - zmortp2 * zfactfe
         tr(ji,jj,jk,jpdsi,Krhs) = tr(ji,jj,jk,jpdsi,Krhs) - zmortp2 * zfactsi
         tr(ji,jj,jk,jpgsi,Krhs) = tr(ji,jj,jk,jpgsi,Krhs) + zmortp2 * zfactsi
         tr(ji,jj,jk,jpgoc,Krhs) = tr(ji,jj,jk,jpgoc,Krhs) + zrespp2 + 0.5 * ztortp2
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + 0.5 * ztortp2
         prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + 0.5 * ztortp2
         prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zrespp2 + 0.5 * ztortp2
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + 0.5 * ztortp2 * zfactfe
         tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + ( zrespp2 + 0.5 * ztortp2 ) * zfactfe
      END_3D
      !
      IF(sn_cfctl%l_prttrc) THEN      ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_diat')
      !
   END SUBROUTINE p4z_diat


   SUBROUTINE p4z_mort_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the nampismort namelist and check the parameters
      !!              called at the first timestep
      !!
      !! ** input   :   Namelist nampismort
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zmort/ wchl, wchld, wchldm, mprat, mprat2
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) 
         WRITE(numout,*) 'p4z_mort_init : Initialization of phytoplankton mortality parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      READ  ( numnatp_ref, namp4zmort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zmort in reference namelist' )
      READ  ( numnatp_cfg, namp4zmort, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zmort in configuration namelist' )
      IF(lwm) WRITE( numonp, namp4zmort )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zmort'
         WRITE(numout,*) '      quadratic mortality of phytoplankton        wchl   =', wchl
         WRITE(numout,*) '      maximum quadratic mortality of diatoms      wchld  =', wchld
         WRITE(numout,*) '      maximum quadratic mortality of diatoms      wchldm =', wchldm
         WRITE(numout,*) '      phytoplankton mortality rate                mprat  =', mprat
         WRITE(numout,*) '      Diatoms mortality rate                      mprat2 =', mprat2
      ENDIF
      !
   END SUBROUTINE p4z_mort_init

   !!======================================================================
END MODULE p4zmort
