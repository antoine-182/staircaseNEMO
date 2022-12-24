MODULE p4zmeso
   !!======================================================================
   !!                         ***  MODULE p4zmeso  ***
   !! TOP :   PISCES Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
   !!   p4z_meso       : Compute the sources/sinks for mesozooplankton
   !!   p4z_meso_init  : Initialization of the parameters for mesozooplankton
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zprod         ! production
   USE prtctl_trc      ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_meso              ! called in p4zbio.F90
   PUBLIC   p4z_meso_init         ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC ::  part2        !: part of calcite not dissolved in mesozoo guts
   REAL(wp), PUBLIC ::  xpref2d      !: mesozoo preference for diatoms
   REAL(wp), PUBLIC ::  xpref2n      !: mesozoo preference for nanophyto
   REAL(wp), PUBLIC ::  xpref2z      !: mesozoo preference for microzooplankton
   REAL(wp), PUBLIC ::  xpref2c      !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xthresh2zoo  !: zoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2dia  !: diatoms feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2phy  !: nanophyto feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2poc  !: poc feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2     !: feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  resrat2      !: exsudation rate of mesozooplankton
   REAL(wp), PUBLIC ::  mzrat2       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat2     !: maximal mesozoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz2      !: non assimilated fraction of P by mesozoo 
   REAL(wp), PUBLIC ::  unass2       !: Efficicency of mesozoo growth 
   REAL(wp), PUBLIC ::  sigma2       !: Fraction of mesozoo excretion as DOM 
   REAL(wp), PUBLIC ::  epsher2      !: growth efficiency
   REAL(wp), PUBLIC ::  epsher2min   !: minimum growth efficiency at high food for grazing 2
   REAL(wp), PUBLIC ::  grazflux     !: mesozoo flux feeding rate

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmeso.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_meso( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step and ???
      INTEGER, INTENT(in)  ::  Kbb, Krhs ! time level indices
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaph, zcompapoc, zcompaz, zcompam
      REAL(wp) :: zgraze2 , zdenom, zdenom2
      REAL(wp) :: zfact   , zfood, zfoodlim, zproport, zbeta
      REAL(wp) :: zmortzgoc, zfrac, zfracfe, zratio, zratio2, zfracal, zgrazcal
      REAL(wp) :: zepsherf, zepshert, zepsherv, zgrarsig, zgraztotc, zgraztotn, zgraztotf
      REAL(wp) :: zgrarem2, zgrafer2, zgrapoc2, zprcaca, zmortz, zgrasrat, zgrasratn
      REAL(wp) :: zrespz, ztortz, zgrazd, zgrazz, zgrazpof
      REAL(wp) :: zgrazn, zgrazpoc, zgraznf, zgrazf
      REAL(wp) :: zgrazfffp, zgrazfffg, zgrazffep, zgrazffeg
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zgrazing2, zfezoo2, zz2ligprod
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_meso')
      !
      DO_3D_11_11( 1, jpkm1 )
         zcompam   = MAX( ( tr(ji,jj,jk,jpmes,Kbb) - 1.e-9 ), 0.e0 )
         zfact     = xstep * tgfunc2(ji,jj,jk) * zcompam

         !  Respiration rates of both zooplankton
         !  -------------------------------------
         zrespz    = resrat2 * zfact * ( tr(ji,jj,jk,jpmes,Kbb) / ( xkmort + tr(ji,jj,jk,jpmes,Kbb) )  &
         &           + 3. * nitrfac(ji,jj,jk) )

         !  Zooplankton mortality. A square function has been selected with
         !  no real reason except that it seems to be more stable and may mimic predation
         !  ---------------------------------------------------------------
         ztortz    = mzrat2 * 1.e6 * zfact * tr(ji,jj,jk,jpmes,Kbb)  * (1. - nitrfac(ji,jj,jk) )
         !
         zcompadi  = MAX( ( tr(ji,jj,jk,jpdia,Kbb) - xthresh2dia ), 0.e0 )
         zcompaz   = MAX( ( tr(ji,jj,jk,jpzoo,Kbb) - xthresh2zoo ), 0.e0 )
         zcompapoc = MAX( ( tr(ji,jj,jk,jppoc,Kbb) - xthresh2poc ), 0.e0 )
         ! Size effect of nanophytoplankton on grazing : the smaller it is, the less prone
         ! it is to predation by mesozooplankton
         ! -------------------------------------------------------------------------------
         zcompaph  = MAX( ( tr(ji,jj,jk,jpphy,Kbb) - xthresh2phy ), 0.e0 ) &
            &      * MIN(1., MAX( 0., ( quotan(ji,jj,jk) - 0.2) / 0.3 ) )

         !   Mesozooplankton grazing
         !   ------------------------
         zfood     = xpref2d * zcompadi + xpref2z * zcompaz + xpref2n * zcompaph + xpref2c * zcompapoc 
         zfoodlim  = MAX( 0., zfood - MIN( 0.5 * zfood, xthresh2 ) )
         zdenom    = zfoodlim / ( xkgraz2 + zfoodlim )
         zdenom2   = zdenom / ( zfood + rtrn )
         zgraze2   = grazrat2 * xstep * tgfunc2(ji,jj,jk) * tr(ji,jj,jk,jpmes,Kbb) * (1. - nitrfac(ji,jj,jk)) 

         zgrazd    = zgraze2  * xpref2d  * zcompadi  * zdenom2 
         zgrazz    = zgraze2  * xpref2z  * zcompaz   * zdenom2 
         zgrazn    = zgraze2  * xpref2n  * zcompaph  * zdenom2 
         zgrazpoc  = zgraze2  * xpref2c  * zcompapoc * zdenom2 

         zgraznf   = zgrazn   * tr(ji,jj,jk,jpnfe,Kbb) / ( tr(ji,jj,jk,jpphy,Kbb) + rtrn)
         zgrazf    = zgrazd   * tr(ji,jj,jk,jpdfe,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn)
         zgrazpof  = zgrazpoc * tr(ji,jj,jk,jpsfe,Kbb) / ( tr(ji,jj,jk,jppoc,Kbb) + rtrn)

         !  Mesozooplankton flux feeding on GOC
         !  ----------------------------------
         zgrazffeg = grazflux  * xstep * wsbio4(ji,jj,jk)      &
         &           * tgfunc2(ji,jj,jk) * tr(ji,jj,jk,jpgoc,Kbb) * tr(ji,jj,jk,jpmes,Kbb) &
         &           * (1. - nitrfac(ji,jj,jk))
         zgrazfffg = zgrazffeg * tr(ji,jj,jk,jpbfe,Kbb) / (tr(ji,jj,jk,jpgoc,Kbb) + rtrn)
         zgrazffep = grazflux  * xstep *  wsbio3(ji,jj,jk)     &
         &           * tgfunc2(ji,jj,jk) * tr(ji,jj,jk,jppoc,Kbb) * tr(ji,jj,jk,jpmes,Kbb) &
         &           * (1. - nitrfac(ji,jj,jk))
         zgrazfffp = zgrazffep * tr(ji,jj,jk,jpsfe,Kbb) / (tr(ji,jj,jk,jppoc,Kbb) + rtrn)
         !
         zgraztotc = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffep + zgrazffeg
         ! Compute the proportion of filter feeders
         zproport  = (zgrazffep + zgrazffeg)/(rtrn + zgraztotc)
         ! Compute fractionation of aggregates. It is assumed that 
         ! diatoms based aggregates are more prone to fractionation
         ! since they are more porous (marine snow instead of fecal pellets)
         zratio    = tr(ji,jj,jk,jpgsi,Kbb) / ( tr(ji,jj,jk,jpgoc,Kbb) + rtrn )
         zratio2   = zratio * zratio
         zfrac     = zproport * grazflux  * xstep * wsbio4(ji,jj,jk)      &
         &          * tr(ji,jj,jk,jpgoc,Kbb) * tr(ji,jj,jk,jpmes,Kbb)          &
         &          * ( 0.2 + 3.8 * zratio2 / ( 1.**2 + zratio2 ) )
         zfracfe   = zfrac * tr(ji,jj,jk,jpbfe,Kbb) / (tr(ji,jj,jk,jpgoc,Kbb) + rtrn)

         zgrazffep = zproport * zgrazffep
         zgrazffeg = zproport * zgrazffeg
         zgrazfffp = zproport * zgrazfffp
         zgrazfffg = zproport * zgrazfffg
         zgraztotc = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffep + zgrazffeg
         zgraztotn = zgrazd * quotad(ji,jj,jk) + zgrazz + zgrazn * quotan(ji,jj,jk)   &
         &   + zgrazpoc + zgrazffep + zgrazffeg
         zgraztotf = zgrazf + zgraznf + zgrazz * ferat3 + zgrazpof + zgrazfffp + zgrazfffg

         ! Total grazing ( grazing by microzoo is already computed in p4zmicro )
         zgrazing2(ji,jj,jk) = zgraztotc

         !    Mesozooplankton efficiency
         !    --------------------------
         zgrasrat  =  ( zgraztotf + rtrn )/ ( zgraztotc + rtrn )
         zgrasratn =  ( zgraztotn + rtrn )/ ( zgraztotc + rtrn )
         zepshert  = MIN( 1., zgrasratn, zgrasrat / ferat3)
         zbeta     = MAX(0., (epsher2 - epsher2min) )
         zepsherf  = epsher2min + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta ) 
         zepsherv  = zepsherf * zepshert 

         zgrarem2  = zgraztotc * ( 1. - zepsherv - unass2 ) &
         &         + ( 1. - epsher2 - unass2 ) / ( 1. - epsher2 ) * ztortz
         zgrafer2  = zgraztotc * MAX( 0. , ( 1. - unass2 ) * zgrasrat - ferat3 * zepsherv )    &
         &         + ferat3 * ( ( 1. - epsher2 - unass2 ) /( 1. - epsher2 ) * ztortz )
         zgrapoc2  = zgraztotc * unass2

         !   Update the arrays TRA which contain the biological sources and sinks
         zgrarsig  = zgrarem2 * sigma2
         tr(ji,jj,jk,jppo4,Krhs) = tr(ji,jj,jk,jppo4,Krhs) + zgrarsig
         tr(ji,jj,jk,jpnh4,Krhs) = tr(ji,jj,jk,jpnh4,Krhs) + zgrarsig
         tr(ji,jj,jk,jpdoc,Krhs) = tr(ji,jj,jk,jpdoc,Krhs) + zgrarem2 - zgrarsig
         !
         IF( ln_ligand ) THEN 
            tr(ji,jj,jk,jplgw,Krhs) = tr(ji,jj,jk,jplgw,Krhs) + (zgrarem2 - zgrarsig) * ldocz
            zz2ligprod(ji,jj,jk) = (zgrarem2 - zgrarsig) * ldocz
         ENDIF
         !
         tr(ji,jj,jk,jpoxy,Krhs) = tr(ji,jj,jk,jpoxy,Krhs) - o2ut * zgrarsig
         tr(ji,jj,jk,jpfer,Krhs) = tr(ji,jj,jk,jpfer,Krhs) + zgrafer2
         zfezoo2(ji,jj,jk)   = zgrafer2
         tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) + zgrarsig
         tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) + rno3 * zgrarsig              

         zmortz = ztortz + zrespz
         zmortzgoc = unass2 / ( 1. - epsher2 ) * ztortz + zrespz
         tr(ji,jj,jk,jpmes,Krhs) = tr(ji,jj,jk,jpmes,Krhs) - zmortz + zepsherv * zgraztotc 
         tr(ji,jj,jk,jpdia,Krhs) = tr(ji,jj,jk,jpdia,Krhs) - zgrazd
         tr(ji,jj,jk,jpzoo,Krhs) = tr(ji,jj,jk,jpzoo,Krhs) - zgrazz
         tr(ji,jj,jk,jpphy,Krhs) = tr(ji,jj,jk,jpphy,Krhs) - zgrazn
         tr(ji,jj,jk,jpnch,Krhs) = tr(ji,jj,jk,jpnch,Krhs) - zgrazn * tr(ji,jj,jk,jpnch,Kbb) / ( tr(ji,jj,jk,jpphy,Kbb) + rtrn )
         tr(ji,jj,jk,jpdch,Krhs) = tr(ji,jj,jk,jpdch,Krhs) - zgrazd * tr(ji,jj,jk,jpdch,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         tr(ji,jj,jk,jpdsi,Krhs) = tr(ji,jj,jk,jpdsi,Krhs) - zgrazd * tr(ji,jj,jk,jpdsi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         tr(ji,jj,jk,jpgsi,Krhs) = tr(ji,jj,jk,jpgsi,Krhs) + zgrazd * tr(ji,jj,jk,jpdsi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         tr(ji,jj,jk,jpnfe,Krhs) = tr(ji,jj,jk,jpnfe,Krhs) - zgraznf
         tr(ji,jj,jk,jpdfe,Krhs) = tr(ji,jj,jk,jpdfe,Krhs) - zgrazf

         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) - zgrazpoc - zgrazffep + zfrac
         prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zfrac
         conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazpoc - zgrazffep
         tr(ji,jj,jk,jpgoc,Krhs) = tr(ji,jj,jk,jpgoc,Krhs) + zmortzgoc - zgrazffeg + zgrapoc2 - zfrac
         prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zmortzgoc + zgrapoc2
         consgoc(ji,jj,jk) = consgoc(ji,jj,jk) - zgrazffeg - zfrac
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) - zgrazpof - zgrazfffp + zfracfe
         tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + ferat3 * zmortzgoc - zgrazfffg     &
           &                + zgraztotf * unass2 - zfracfe
         zfracal = tr(ji,jj,jk,jpcal,Kbb) / (tr(ji,jj,jk,jppoc,Kbb) + tr(ji,jj,jk,jpgoc,Kbb) + rtrn )
         zgrazcal = (zgrazffeg + zgrazpoc) * (1. - part2) * zfracal
         ! calcite production
         zprcaca = xfracal(ji,jj,jk) * zgrazn
         prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         !
         zprcaca = part2 * zprcaca
         tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) + zgrazcal - zprcaca
         tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) - 2. * ( zgrazcal + zprcaca )
         tr(ji,jj,jk,jpcal,Krhs) = tr(ji,jj,jk,jpcal,Krhs) - zgrazcal + zprcaca
      END_3D
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
        CALL iom_put( "PCAL"  , prodcal(:,:,:) * 1.e+3  * rfact2r * tmask(:,:,:) )  !  Calcite production 
        IF( iom_use("GRAZ2") ) THEN  !   Total grazing of phyto by zooplankton
           zgrazing2(:,:,jpk) = 0._wp ;  CALL iom_put( "GRAZ2" , zgrazing2(:,:,:) * 1.e+3  * rfact2r * tmask(:,:,:) ) 
         ENDIF
         IF( iom_use("FEZOO2") ) THEN  
           zfezoo2 (:,:,jpk) = 0._wp ; CALL iom_put( "FEZOO2", zfezoo2(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:) )
         ENDIF
         IF( ln_ligand ) THEN
            zz2ligprod(:,:,jpk) = 0._wp ; CALL iom_put( "LPRODZ2", zz2ligprod(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)  )
         ENDIF
      ENDIF
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_meso')
      !
   END SUBROUTINE p4z_meso


   SUBROUTINE p4z_meso_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zmes/ part2, grazrat2, resrat2, mzrat2, xpref2n, xpref2d, xpref2z,   &
         &                xpref2c, xthresh2dia, xthresh2phy, xthresh2zoo, xthresh2poc, &
         &                xthresh2, xkgraz2, epsher2, epsher2min, sigma2, unass2, grazflux
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) 
         WRITE(numout,*) 'p4z_meso_init : Initialization of mesozooplankton parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      READ  ( numnatp_ref, namp4zmes, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zmes in reference namelist' )
      READ  ( numnatp_cfg, namp4zmes, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zmes in configuration namelist' )
      IF(lwm) WRITE( numonp, namp4zmes )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zmes'
         WRITE(numout,*) '      part of calcite not dissolved in mesozoo guts  part2        =', part2
         WRITE(numout,*) '      mesozoo preference for phyto                   xpref2n      =', xpref2n
         WRITE(numout,*) '      mesozoo preference for diatoms                 xpref2d      =', xpref2d
         WRITE(numout,*) '      mesozoo preference for zoo                     xpref2z      =', xpref2z
         WRITE(numout,*) '      mesozoo preference for poc                     xpref2c      =', xpref2c
         WRITE(numout,*) '      microzoo feeding threshold  for mesozoo        xthresh2zoo  =', xthresh2zoo
         WRITE(numout,*) '      diatoms feeding threshold  for mesozoo         xthresh2dia  =', xthresh2dia
         WRITE(numout,*) '      nanophyto feeding threshold for mesozoo        xthresh2phy  =', xthresh2phy
         WRITE(numout,*) '      poc feeding threshold for mesozoo              xthresh2poc  =', xthresh2poc
         WRITE(numout,*) '      feeding threshold for mesozooplankton          xthresh2     =', xthresh2
         WRITE(numout,*) '      exsudation rate of mesozooplankton             resrat2      =', resrat2
         WRITE(numout,*) '      mesozooplankton mortality rate                 mzrat2       =', mzrat2
         WRITE(numout,*) '      maximal mesozoo grazing rate                   grazrat2     =', grazrat2
         WRITE(numout,*) '      mesozoo flux feeding rate                      grazflux     =', grazflux
         WRITE(numout,*) '      non assimilated fraction of P by mesozoo       unass2       =', unass2
         WRITE(numout,*) '      Efficiency of Mesozoo growth                   epsher2      =', epsher2
         WRITE(numout,*) '      Minimum Efficiency of Mesozoo growth           epsher2min  =', epsher2min
         WRITE(numout,*) '      Fraction of mesozoo excretion as DOM           sigma2       =', sigma2
         WRITE(numout,*) '      half sturation constant for grazing 2          xkgraz2      =', xkgraz2
      ENDIF
      !
   END SUBROUTINE p4z_meso_init

   !!======================================================================
END MODULE p4zmeso
