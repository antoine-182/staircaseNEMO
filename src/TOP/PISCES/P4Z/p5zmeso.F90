MODULE p5zmeso
   !!======================================================================
   !!                         ***  MODULE p5zmeso  ***
   !! TOP :   PISCES Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p5z_meso       :   Compute the sources/sinks for mesozooplankton
   !!   p5z_meso_init  :   Initialization of the parameters for mesozooplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p5z_meso              ! called in p5zbio.F90
   PUBLIC   p5z_meso_init         ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part2        !: part of calcite not dissolved in mesozoo guts
   REAL(wp), PUBLIC ::  xpref2c      !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xpref2n      !: mesozoo preference for nanophyto
   REAL(wp), PUBLIC ::  xpref2z      !: mesozoo preference for zooplankton
   REAL(wp), PUBLIC ::  xpref2d      !: mesozoo preference for Diatoms 
   REAL(wp), PUBLIC ::  xpref2m      !: mesozoo preference for mesozoo
   REAL(wp), PUBLIC ::  xthresh2zoo  !: zoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2dia  !: diatoms feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2phy  !: nanophyto feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2poc  !: poc feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2mes  !: mesozoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2     !: feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  resrat2      !: exsudation rate of mesozooplankton
   REAL(wp), PUBLIC ::  mzrat2       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat2     !: maximal mesozoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz2      !: Half-saturation constant of assimilation
   REAL(wp), PUBLIC ::  unass2c      !: Non-assimilated fraction of food
   REAL(wp), PUBLIC ::  unass2n      !: Non-assimilated fraction of food
   REAL(wp), PUBLIC ::  unass2p      !: Non-assimilated fraction of food
   REAL(wp), PUBLIC ::  epsher2      !: Growth efficiency of mesozoo
   REAL(wp), PUBLIC ::  epsher2min   !: Minimum growth efficiency of mesozoo
   REAL(wp), PUBLIC ::  ssigma2      !: Fraction excreted as semi-labile DOM
   REAL(wp), PUBLIC ::  srespir2     !: Active respiration
   REAL(wp), PUBLIC ::  grazflux     !: mesozoo flux feeding rate
   LOGICAL,  PUBLIC ::  bmetexc2     !: Use of excess carbon for respiration

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p5zmeso.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p5z_meso( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt    ! ocean time step
      INTEGER, INTENT(in)  ::  Kbb, Krhs  ! time level indices
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaph, zcompapoc, zcompaz, zcompam, zcompames
      REAL(wp) :: zgraze2, zdenom, zfact, zfood, zfoodlim, zproport
      REAL(wp) :: zmortzgoc, zfracc, zfracn, zfracp, zfracfe, zratio, zratio2
      REAL(wp) :: zepsherf, zepshert, zepsherv, zrespirc, zrespirn, zrespirp, zbasresb, zbasresi
      REAL(wp) :: zgraztotc, zgraztotn, zgraztotp, zgraztotf, zbasresn, zbasresp, zbasresf
      REAL(wp) :: zgradoc, zgradon, zgradop, zgratmp, zgradoct, zgradont, zgrareft, zgradopt
      REAL(wp) :: zgrapoc, zgrapon, zgrapop, zgrapof, zprcaca, zmortz
      REAL(wp) :: zexcess, zgrarem, zgraren, zgrarep, zgraref
      REAL(wp) :: zbeta, zrespz, ztortz, zgrasratp, zgrasratn, zgrasratf
      REAL(wp) :: ztmp1, ztmp2, ztmp3, ztmp4, ztmp5, ztmptot
      REAL(wp) :: zgrazdc, zgrazz, zgrazm, zgrazpof, zgrazcal, zfracal
      REAL(wp) :: zgraznc, zgrazpoc, zgrazpon, zgrazpop, zgraznf, zgrazdf
      REAL(wp) :: zgraznp, zgraznn, zgrazdn, zgrazdp
      REAL(wp) :: zgrazfffp, zgrazfffg, zgrazffep, zgrazffeg
      REAL(wp) :: zgrazffnp, zgrazffng, zgrazffpp, zgrazffpg
      CHARACTER (len=25) :: charout
      REAL(wp) :: zrfact2, zmetexcess
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zgrazing, zfezoo2,  zz2ligprod

      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p5z_meso')
      !
      zmetexcess = 0.0
      IF ( bmetexc2 ) zmetexcess = 1.0

      DO_3D_11_11( 1, jpkm1 )
         zcompam   = MAX( ( tr(ji,jj,jk,jpmes,Kbb) - 1.e-9 ), 0.e0 )
         zfact     = xstep * tgfunc2(ji,jj,jk) * zcompam

         !   Michaelis-Menten mortality rates of mesozooplankton
         !   ---------------------------------------------------
         zrespz   = resrat2 * zfact * ( tr(ji,jj,jk,jpmes,Kbb) / ( xkmort + tr(ji,jj,jk,jpmes,Kbb) )  &
         &          + 3. * nitrfac(ji,jj,jk) )

         !   Zooplankton mortality. A square function has been selected with
         !   no real reason except that it seems to be more stable and may mimic predation
         !   ---------------------------------------------------------------
         ztortz   = mzrat2 * 1.e6 * zfact * tr(ji,jj,jk,jpmes,Kbb) * (1. - nitrfac(ji,jj,jk))

         !   Computation of the abundance of the preys
         !   A threshold can be specified in the namelist
         !   --------------------------------------------
         zcompadi  = MAX( ( tr(ji,jj,jk,jpdia,Kbb) - xthresh2dia ), 0.e0 )
         zcompaz   = MAX( ( tr(ji,jj,jk,jpzoo,Kbb) - xthresh2zoo ), 0.e0 )
         zcompaph  = MAX( ( tr(ji,jj,jk,jpphy,Kbb) - xthresh2phy ), 0.e0 )
         zcompapoc = MAX( ( tr(ji,jj,jk,jppoc,Kbb) - xthresh2poc ), 0.e0 )
         zcompames = MAX( ( tr(ji,jj,jk,jpmes,Kbb) - xthresh2mes ), 0.e0 )

         !   Mesozooplankton grazing
         !   ------------------------
         zfood     = xpref2d * zcompadi + xpref2z * zcompaz + xpref2n * zcompaph + xpref2c * zcompapoc   &
         &           + xpref2m * zcompames 
         zfoodlim  = MAX( 0., zfood - MIN( 0.5 * zfood, xthresh2 ) )
         zdenom    = zfoodlim / ( xkgraz2 + zfoodlim )
         zgraze2   = grazrat2 * xstep * tgfunc2(ji,jj,jk) * tr(ji,jj,jk,jpmes,Kbb) * (1. - nitrfac(ji,jj,jk)) 

         !   An active switching parameterization is used here.
         !   We don't use the KTW parameterization proposed by 
         !   Vallina et al. because it tends to produce to steady biomass
         !   composition and the variance of Chl is too low as it grazes
         !   too strongly on winning organisms. Thus, instead of a square
         !   a 1.5 power value is used which decreases the pressure on the
         !   most abundant species
         !   ------------------------------------------------------------  
         ztmp1 = xpref2n * zcompaph**1.5
         ztmp2 = xpref2m * zcompames**1.5
         ztmp3 = xpref2c * zcompapoc**1.5
         ztmp4 = xpref2d * zcompadi**1.5
         ztmp5 = xpref2z * zcompaz**1.5
         ztmptot = ztmp1 + ztmp2 + ztmp3 + ztmp4 + ztmp5 + rtrn
         ztmp1 = ztmp1 / ztmptot
         ztmp2 = ztmp2 / ztmptot
         ztmp3 = ztmp3 / ztmptot
         ztmp4 = ztmp4 / ztmptot
         ztmp5 = ztmp5 / ztmptot

         !   Mesozooplankton regular grazing on the different preys
         !   ------------------------------------------------------
         zgrazdc   = zgraze2 * ztmp4 * zdenom
         zgrazdn   = zgrazdc * tr(ji,jj,jk,jpndi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn)
         zgrazdp   = zgrazdc * tr(ji,jj,jk,jppdi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn)
         zgrazdf   = zgrazdc * tr(ji,jj,jk,jpdfe,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn)
         zgrazz    = zgraze2 * ztmp5 * zdenom
         zgrazm    = zgraze2 * ztmp2 * zdenom
         zgraznc   = zgraze2 * ztmp1 * zdenom
         zgraznn   = zgraznc * tr(ji,jj,jk,jpnph,Kbb) / ( tr(ji,jj,jk,jpphy,Kbb) + rtrn)
         zgraznp   = zgraznc * tr(ji,jj,jk,jppph,Kbb) / ( tr(ji,jj,jk,jpphy,Kbb) + rtrn)
         zgraznf   = zgraznc * tr(ji,jj,jk,jpnfe,Kbb) / ( tr(ji,jj,jk,jpphy,Kbb) + rtrn)
         zgrazpoc  = zgraze2 * ztmp3 * zdenom
         zgrazpon  = zgrazpoc * tr(ji,jj,jk,jppon,Kbb) / ( tr(ji,jj,jk,jppoc,Kbb) + rtrn)
         zgrazpop  = zgrazpoc * tr(ji,jj,jk,jppop,Kbb) / ( tr(ji,jj,jk,jppoc,Kbb) + rtrn)
         zgrazpof  = zgrazpoc * tr(ji,jj,jk,jpsfe,Kbb) / ( tr(ji,jj,jk,jppoc,Kbb) + rtrn)

         !   Mesozooplankton flux feeding on GOC
         !   ----------------------------------
         zgrazffeg = grazflux  * xstep * wsbio4(ji,jj,jk)      &
         &           * tgfunc2(ji,jj,jk) * tr(ji,jj,jk,jpgoc,Kbb) * tr(ji,jj,jk,jpmes,Kbb)  &
         &           * (1. - nitrfac(ji,jj,jk))
         zgrazfffg = zgrazffeg * tr(ji,jj,jk,jpbfe,Kbb) / (tr(ji,jj,jk,jpgoc,Kbb) + rtrn)
         zgrazffng = zgrazffeg * tr(ji,jj,jk,jpgon,Kbb) / (tr(ji,jj,jk,jpgoc,Kbb) + rtrn)
         zgrazffpg = zgrazffeg * tr(ji,jj,jk,jpgop,Kbb) / (tr(ji,jj,jk,jpgoc,Kbb) + rtrn)
         zgrazffep = grazflux  * xstep *  wsbio3(ji,jj,jk)     &
         &           * tgfunc2(ji,jj,jk) * tr(ji,jj,jk,jppoc,Kbb) * tr(ji,jj,jk,jpmes,Kbb)   &
         &           * (1. - nitrfac(ji,jj,jk))
         zgrazfffp = zgrazffep * tr(ji,jj,jk,jpsfe,Kbb) / (tr(ji,jj,jk,jppoc,Kbb) + rtrn)
         zgrazffnp = zgrazffep * tr(ji,jj,jk,jppon,Kbb) / (tr(ji,jj,jk,jppoc,Kbb) + rtrn)
         zgrazffpp = zgrazffep * tr(ji,jj,jk,jppop,Kbb) / (tr(ji,jj,jk,jppoc,Kbb) + rtrn)
         !
         zgraztotc  = zgrazdc + zgrazz + zgraznc + zgrazm + zgrazpoc + zgrazffep + zgrazffeg

         !   Compute the proportion of filter feeders
         !   ----------------------------------------  
         zproport  = (zgrazffep + zgrazffeg)/(rtrn + zgraztotc)

         !   Compute fractionation of aggregates. It is assumed that 
         !   diatoms based aggregates are more prone to fractionation
         !   since they are more porous (marine snow instead of fecal pellets)
         !   ----------------------------------------------------------------
         zratio    = tr(ji,jj,jk,jpgsi,Kbb) / ( tr(ji,jj,jk,jpgoc,Kbb) + rtrn )
         zratio2   = zratio * zratio
         zfracc    = zproport * grazflux  * xstep * wsbio4(ji,jj,jk)      &
         &          * tr(ji,jj,jk,jpgoc,Kbb) * tr(ji,jj,jk,jpmes,Kbb)          &
         &          * ( 0.2 + 3.8 * zratio2 / ( 1.**2 + zratio2 ) )
         zfracfe   = zfracc * tr(ji,jj,jk,jpbfe,Kbb) / (tr(ji,jj,jk,jpgoc,Kbb) + rtrn)
         zfracn    = zfracc * tr(ji,jj,jk,jpgon,Kbb) / (tr(ji,jj,jk,jpgoc,Kbb) + rtrn)
         zfracp    = zfracc * tr(ji,jj,jk,jpgop,Kbb) / (tr(ji,jj,jk,jpgoc,Kbb) + rtrn)

         zgrazffep = zproport * zgrazffep   ;   zgrazffeg = zproport * zgrazffeg
         zgrazfffp = zproport * zgrazfffp   ;   zgrazfffg = zproport * zgrazfffg
         zgrazffnp = zproport * zgrazffnp   ;   zgrazffng = zproport * zgrazffng
         zgrazffpp = zproport * zgrazffpp   ;   zgrazffpg = zproport * zgrazffpg

         zgraztotc  = zgrazdc + zgrazz + zgraznc + zgrazm + zgrazpoc + zgrazffep + zgrazffeg
         zgraztotf  = zgrazdf + zgraznf + ( zgrazz + zgrazm ) * ferat3 + zgrazpof &
         &            + zgrazfffp + zgrazfffg
         zgraztotn  = zgrazdn + (zgrazm + zgrazz) * no3rat3 + zgraznn + zgrazpon  &
         &            + zgrazffnp + zgrazffng
         zgraztotp  = zgrazdp + (zgrazz + zgrazm) * po4rat3 + zgraznp + zgrazpop  &
         &            + zgrazffpp + zgrazffpg


         ! Total grazing ( grazing by microzoo is already computed in p5zmicro )
         zgrazing(ji,jj,jk) = zgraztotc

         !   Stoichiometruc ratios of the food ingested by zooplanton 
         !   --------------------------------------------------------
         zgrasratf  =  (zgraztotf + rtrn) / ( zgraztotc + rtrn )
         zgrasratn  =  (zgraztotn + rtrn) / ( zgraztotc + rtrn )
         zgrasratp  =  (zgraztotp + rtrn) / ( zgraztotc + rtrn )

         !   Growth efficiency is made a function of the quality 
         !   and the quantity of the preys
         !   ---------------------------------------------------
         zepshert  = MIN( 1., zgrasratn/ no3rat3, zgrasratp/ po4rat3, zgrasratf / ferat3)
         zbeta     = MAX(0., (epsher2 - epsher2min) )
         zepsherf  = epsher2min + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
         zepsherv  = zepsherf * zepshert

         !   Respiration of mesozooplankton
         !   Excess carbon in the food is used preferentially
         !   ----------------  ------------------------------
         zexcess  = zgraztotc * zepsherf * (1.0 - zepshert) * zmetexcess 
         zbasresb = MAX(0., zrespz - zexcess)
         zbasresi = zexcess + MIN(0., zrespz - zexcess)
         zrespirc = srespir2 * zepsherv * zgraztotc + zbasresb

         !   When excess carbon is used, the other elements in excess
         !   are also used proportionally to their abundance
         !   --------------------------------------------------------
         zexcess  = ( zgrasratn/ no3rat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
         zbasresn = zbasresi * zexcess * zgrasratn
         zexcess  = ( zgrasratp/ po4rat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
         zbasresp = zbasresi * zexcess * zgrasratp
         zexcess  = ( zgrasratf/ ferat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
         zbasresf = zbasresi * zexcess * zgrasratf

         !   Voiding of the excessive elements as organic matter
         !   --------------------------------------------------------
         zgradoct = (1. - unass2c - zepsherv) * zgraztotc - zbasresi
         zgradont = (1. - unass2n) * zgraztotn - zepsherv * no3rat3 * zgraztotc - zbasresn
         zgradopt = (1. - unass2p) * zgraztotp - zepsherv * po4rat3 * zgraztotc - zbasresp
         zgrareft = (1. - unass2c) * zgraztotf - zepsherv * ferat3 * zgraztotc - zbasresf
         ztmp1   = ( 1. - epsher2 - unass2c ) /( 1. - 0.8 * epsher2 ) * ztortz
         zgradoc = (zgradoct + ztmp1) * ssigma2
         zgradon = (zgradont + no3rat3 * ztmp1) * ssigma2
         zgradop = (zgradopt + po4rat3 * ztmp1) * ssigma2
         zgratmp = 0.2 * epsher2 /( 1. - 0.8 * epsher2 ) * ztortz

         !  Since only semilabile DOM is represented in PISCES
         !  part of DOM is in fact labile and is then released
         !  as dissolved inorganic compounds (ssigma2)
         !  --------------------------------------------------
         zgrarem = zgratmp + ( zgradoct + ztmp1 ) * (1.0 - ssigma2)
         zgraren = no3rat3 * zgratmp + ( zgradont + no3rat3 * ztmp1 ) * (1.0 - ssigma2)
         zgrarep = po4rat3 * zgratmp + ( zgradopt + po4rat3 * ztmp1 ) * (1.0 - ssigma2)
         zgraref = zgrareft + ferat3 * ( ztmp1 + zgratmp )

         !   Defecation as a result of non assimilated products
         !   --------------------------------------------------
         zgrapoc  = zgraztotc * unass2c + unass2c / ( 1. - 0.8 * epsher2 ) * ztortz
         zgrapon  = zgraztotn * unass2n + no3rat3 * unass2n / ( 1. - 0.8 * epsher2 ) * ztortz
         zgrapop  = zgraztotp * unass2p + po4rat3 * unass2p / ( 1. - 0.8 * epsher2 ) * ztortz
         zgrapof  = zgraztotf * unass2c + ferat3  * unass2c / ( 1. - 0.8 * epsher2 ) * ztortz

         !  Addition of respiration to the release of inorganic nutrients
         !  -------------------------------------------------------------
         zgrarem = zgrarem + zbasresi + zrespirc
         zgraren = zgraren + zbasresn + zrespirc * no3rat3
         zgrarep = zgrarep + zbasresp + zrespirc * po4rat3
         zgraref = zgraref + zbasresf + zrespirc * ferat3

         !   Update the arrays TRA which contain the biological sources and
         !   sinks
         !   --------------------------------------------------------------
         tr(ji,jj,jk,jppo4,Krhs) = tr(ji,jj,jk,jppo4,Krhs) + zgrarep 
         tr(ji,jj,jk,jpnh4,Krhs) = tr(ji,jj,jk,jpnh4,Krhs) + zgraren
         tr(ji,jj,jk,jpdoc,Krhs) = tr(ji,jj,jk,jpdoc,Krhs) + zgradoc
         !
         IF( ln_ligand ) THEN
            tr(ji,jj,jk,jplgw,Krhs)  = tr(ji,jj,jk,jplgw,Krhs) + zgradoc * ldocz
            zz2ligprod(ji,jj,jk) = zgradoc * ldocz
         ENDIF
         !
         tr(ji,jj,jk,jpdon,Krhs) = tr(ji,jj,jk,jpdon,Krhs) + zgradon
         tr(ji,jj,jk,jpdop,Krhs) = tr(ji,jj,jk,jpdop,Krhs) + zgradop
         tr(ji,jj,jk,jpoxy,Krhs) = tr(ji,jj,jk,jpoxy,Krhs) - o2ut * zgrarem
         tr(ji,jj,jk,jpfer,Krhs) = tr(ji,jj,jk,jpfer,Krhs) + zgraref
         zfezoo2(ji,jj,jk)   = zgraref
         tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) + zgrarem
         tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) + rno3 * zgraren
         tr(ji,jj,jk,jpmes,Krhs) = tr(ji,jj,jk,jpmes,Krhs) + zepsherv * zgraztotc - zrespirc   &
         &                     - ztortz - zgrazm
         tr(ji,jj,jk,jpdia,Krhs) = tr(ji,jj,jk,jpdia,Krhs) - zgrazdc
         tr(ji,jj,jk,jpndi,Krhs) = tr(ji,jj,jk,jpndi,Krhs) - zgrazdn
         tr(ji,jj,jk,jppdi,Krhs) = tr(ji,jj,jk,jppdi,Krhs) - zgrazdp
         tr(ji,jj,jk,jpdfe,Krhs) = tr(ji,jj,jk,jpdfe,Krhs) - zgrazdf
         tr(ji,jj,jk,jpzoo,Krhs) = tr(ji,jj,jk,jpzoo,Krhs) - zgrazz
         tr(ji,jj,jk,jpphy,Krhs) = tr(ji,jj,jk,jpphy,Krhs) - zgraznc
         tr(ji,jj,jk,jpnph,Krhs) = tr(ji,jj,jk,jpnph,Krhs) - zgraznn
         tr(ji,jj,jk,jppph,Krhs) = tr(ji,jj,jk,jppph,Krhs) - zgraznp
         tr(ji,jj,jk,jpnfe,Krhs) = tr(ji,jj,jk,jpnfe,Krhs) - zgraznf
         tr(ji,jj,jk,jpnch,Krhs) = tr(ji,jj,jk,jpnch,Krhs) - zgraznc * tr(ji,jj,jk,jpnch,Kbb) / ( tr(ji,jj,jk,jpphy,Kbb) + rtrn )
         tr(ji,jj,jk,jpdch,Krhs) = tr(ji,jj,jk,jpdch,Krhs) - zgrazdc * tr(ji,jj,jk,jpdch,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         tr(ji,jj,jk,jpdsi,Krhs) = tr(ji,jj,jk,jpdsi,Krhs) - zgrazdc * tr(ji,jj,jk,jpdsi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )
         tr(ji,jj,jk,jpgsi,Krhs) = tr(ji,jj,jk,jpgsi,Krhs) + zgrazdc * tr(ji,jj,jk,jpdsi,Kbb) / ( tr(ji,jj,jk,jpdia,Kbb) + rtrn )

         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) - zgrazpoc - zgrazffep + zfracc
         prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zfracc
         conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazpoc - zgrazffep
         tr(ji,jj,jk,jppon,Krhs) = tr(ji,jj,jk,jppon,Krhs) - zgrazpon - zgrazffnp + zfracn
         tr(ji,jj,jk,jppop,Krhs) = tr(ji,jj,jk,jppop,Krhs) - zgrazpop - zgrazffpp + zfracp
         tr(ji,jj,jk,jpgoc,Krhs) = tr(ji,jj,jk,jpgoc,Krhs) - zgrazffeg + zgrapoc - zfracc
         prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zgrapoc
         consgoc(ji,jj,jk) = consgoc(ji,jj,jk) - zgrazffeg - zfracc
         tr(ji,jj,jk,jpgon,Krhs) = tr(ji,jj,jk,jpgon,Krhs) - zgrazffng + zgrapon - zfracn
         tr(ji,jj,jk,jpgop,Krhs) = tr(ji,jj,jk,jpgop,Krhs) - zgrazffpg + zgrapop - zfracp
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) - zgrazpof - zgrazfffp + zfracfe
         tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) - zgrazfffg + zgrapof - zfracfe
         zfracal = tr(ji,jj,jk,jpcal,Kbb) / ( tr(ji,jj,jk,jpgoc,Kbb) + rtrn )
         zgrazcal = zgrazffeg * (1. - part2) * zfracal

         !  calcite production
         !  ------------------
         zprcaca = xfracal(ji,jj,jk) * zgraznc
         prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         zprcaca = part2 * zprcaca
         tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) + zgrazcal - zprcaca
         tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) + 2. * ( zgrazcal - zprcaca )
         tr(ji,jj,jk,jpcal,Krhs) = tr(ji,jj,jk,jpcal,Krhs) - zgrazcal + zprcaca
      END_3D
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
        CALL iom_put( "PCAL"  , prodcal(:,:,:) * 1.e+3  * rfact2r * tmask(:,:,:) )  !  Calcite production 
        IF( iom_use("GRAZ2") ) THEN  !   Total grazing of phyto by zooplankton
           zgrazing(:,:,jpk) = 0._wp ;  CALL iom_put( "GRAZ2" , zgrazing(:,:,:) * 1.e+3  * rfact2r * tmask(:,:,:) ) 
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
      IF( ln_timing )   CALL timing_stop('p5z_meso')
      !
   END SUBROUTINE p5z_meso


   SUBROUTINE p5z_meso_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p5z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namp5zmes/part2, bmetexc2, grazrat2, resrat2, mzrat2, xpref2c, xpref2n, xpref2z, &
         &                xpref2m, xpref2d, xthresh2dia, xthresh2phy, xthresh2zoo, xthresh2poc, &
         &                xthresh2mes, xthresh2, xkgraz2, epsher2, epsher2min, ssigma2, unass2c, &
         &                unass2n, unass2p, srespir2, grazflux
      !!----------------------------------------------------------------------
      !
      READ  ( numnatp_ref, namp5zmes, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampismes in reference namelist' )
      !
      READ  ( numnatp_cfg, namp5zmes, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'nampismes in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp5zmes )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' ' 
         WRITE(numout,*) ' Namelist parameters for mesozooplankton, namp5zmes'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in mesozoo guts  part2       = ', part2
         WRITE(numout,*) '    mesozoo preference for nano.                   xpref2n     = ', xpref2n
         WRITE(numout,*) '    mesozoo preference for diatoms                 xpref2d     = ', xpref2d
         WRITE(numout,*) '    mesozoo preference for zoo                     xpref2z     = ', xpref2z
         WRITE(numout,*) '    mesozoo preference for mesozoo                 xpref2m     = ', xpref2m
         WRITE(numout,*) '    mesozoo preference for poc                     xpref2c     = ', xpref2c
         WRITE(numout,*) '    microzoo feeding threshold  for mesozoo        xthresh2zoo = ', xthresh2zoo
         WRITE(numout,*) '    diatoms feeding threshold  for mesozoo         xthresh2dia = ', xthresh2dia
         WRITE(numout,*) '    nanophyto feeding threshold for mesozoo        xthresh2phy = ', xthresh2phy
         WRITE(numout,*) '    poc feeding threshold for mesozoo              xthresh2poc = ', xthresh2poc
         WRITE(numout,*) '    mesozoo feeding threshold for mesozoo          xthresh2mes = ', xthresh2mes
         WRITE(numout,*) '    feeding threshold for mesozooplankton          xthresh2    = ', xthresh2
         WRITE(numout,*) '    exsudation rate of mesozooplankton             resrat2     = ', resrat2
         WRITE(numout,*) '    mesozooplankton mortality rate                 mzrat2      = ', mzrat2
         WRITE(numout,*) '    maximal mesozoo grazing rate                   grazrat2    = ', grazrat2
         WRITE(numout,*) '    mesozoo flux feeding rate                      grazflux    = ', grazflux
         WRITE(numout,*) '    C egested fraction of food by mesozoo          unass2c     = ', unass2c
         WRITE(numout,*) '    N egested fraction of food by mesozoo          unass2n     = ', unass2n
         WRITE(numout,*) '    P egested fraction of food by mesozoo          unass2p     = ', unass2p
         WRITE(numout,*) '    Efficicency of Mesozoo growth                  epsher2     = ', epsher2
         WRITE(numout,*) '    Minimum Efficiency of Mesozoo growth           epsher2min  =', epsher2min
         WRITE(numout,*) '    Fraction excreted as semi-labile DOM           ssigma2     = ', ssigma2
         WRITE(numout,*) '    Active respiration                             srespir2    = ', srespir2
         WRITE(numout,*) '    half sturation constant for grazing 2          xkgraz2     = ', xkgraz2
         WRITE(numout,*) '    Use excess carbon for respiration              bmetexc2    = ', bmetexc2
      ENDIF
      !
   END SUBROUTINE p5z_meso_init

   !!======================================================================
END MODULE p5zmeso
