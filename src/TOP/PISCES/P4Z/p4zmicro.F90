MODULE p4zmicro
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :   PISCES Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
   !!   p4z_micro      : Compute the sources/sinks for microzooplankton
   !!   p4z_micro_init : Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zlim          ! Co-limitations
   USE p4zprod         ! production
   USE iom             ! I/O manager
   USE prtctl_trc      ! print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC ::   part        !: part of calcite not dissolved in microzoo guts
   REAL(wp), PUBLIC ::   xprefc      !: microzoo preference for POC 
   REAL(wp), PUBLIC ::   xprefn      !: microzoo preference for nanophyto
   REAL(wp), PUBLIC ::   xprefd      !: microzoo preference for diatoms
   REAL(wp), PUBLIC ::   xthreshdia  !: diatoms feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthreshphy  !: nanophyto threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthreshpoc  !: poc threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthresh     !: feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::   resrat      !: exsudation rate of microzooplankton
   REAL(wp), PUBLIC ::   mzrat       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::   grazrat     !: maximal microzoo grazing rate
   REAL(wp), PUBLIC ::   xkgraz      !: Half-saturation constant of assimilation
   REAL(wp), PUBLIC ::   unass       !: Non-assimilated part of food
   REAL(wp), PUBLIC ::   sigma1      !: Fraction of microzoo excretion as DOM 
   REAL(wp), PUBLIC ::   epsher      !: growth efficiency for grazing 1 
   REAL(wp), PUBLIC ::   epshermin   !: minimum growth efficiency for grazing 1

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmicro.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_micro( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time step
      INTEGER, INTENT(in) ::   knt   ! ??? 
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaz , zcompaph, zcompapoc
      REAL(wp) :: zgraze  , zdenom, zdenom2
      REAL(wp) :: zfact   , zfood, zfoodlim, zbeta
      REAL(wp) :: zepsherf, zepshert, zepsherv, zgrarsig, zgraztotc, zgraztotn, zgraztotf
      REAL(wp) :: zgrarem, zgrafer, zgrapoc, zprcaca, zmortz
      REAL(wp) :: zrespz, ztortz, zgrasrat, zgrasratn
      REAL(wp) :: zgrazp, zgrazm, zgrazsd
      REAL(wp) :: zgrazmf, zgrazsf, zgrazpf
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zgrazing, zfezoo, zzligprod
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_micro')
      !
      DO_3D_11_11( 1, jpkm1 )
         zcompaz = MAX( ( tr(ji,jj,jk,jpzoo,Kbb) - 1.e-9 ), 0.e0 )
         zfact   = xstep * tgfunc2(ji,jj,jk) * zcompaz

         !  Respiration rates of both zooplankton
         !  -------------------------------------
         zrespz = resrat * zfact * tr(ji,jj,jk,jpzoo,Kbb) / ( xkmort + tr(ji,jj,jk,jpzoo,Kbb) )  &
            &   + resrat * zfact * 3. * nitrfac(ji,jj,jk)

         !  Zooplankton mortality. A square function has been selected with
         !  no real reason except that it seems to be more stable and may mimic predation.
         !  ---------------------------------------------------------------
         ztortz = mzrat * 1.e6 * zfact * tr(ji,jj,jk,jpzoo,Kbb) * (1. - nitrfac(ji,jj,jk))

         zcompadi  = MIN( MAX( ( tr(ji,jj,jk,jpdia,Kbb) - xthreshdia ), 0.e0 ), xsizedia )
         zcompaph  = MAX( ( tr(ji,jj,jk,jpphy,Kbb) - xthreshphy ), 0.e0 )
         zcompapoc = MAX( ( tr(ji,jj,jk,jppoc,Kbb) - xthreshpoc ), 0.e0 )
         
         !     Microzooplankton grazing
         !     ------------------------
         zfood     = xprefn * zcompaph + xprefc * zcompapoc + xprefd * zcompadi
         zfoodlim  = MAX( 0. , zfood - min(xthresh,0.5*zfood) )
         zdenom    = zfoodlim / ( xkgraz + zfoodlim )
         zdenom2   = zdenom / ( zfood + rtrn )
         zgraze    = grazrat * xstep * tgfunc2(ji,jj,jk) * tr(ji,jj,jk,jpzoo,Kbb) * (1. - nitrfac(ji,jj,jk))

         zgrazp    = zgraze  * xprefn * zcompaph  * zdenom2 
         zgrazm    = zgraze  * xprefc * zcompapoc * zdenom2 
         zgrazsd   = zgraze  * xprefd * zcompadi  * zdenom2 

         zgrazpf   = zgrazp  * tr(ji,jj,jk,jpnfe,Kbb) / (tr(ji,jj,jk,jpphy,Kbb) + rtrn)
         zgrazmf   = zgrazm  * tr(ji,jj,jk,jpsfe,Kbb) / (tr(ji,jj,jk,jppoc,Kbb) + rtrn)
         zgrazsf   = zgrazsd * tr(ji,jj,jk,jpdfe,Kbb) / (tr(ji,jj,jk,jpdia,Kbb) + rtrn)
         !
         zgraztotc = zgrazp  + zgrazm  + zgrazsd 
         zgraztotf = zgrazpf + zgrazsf + zgrazmf 
         zgraztotn = zgrazp * quotan(ji,jj,jk) + zgrazm + zgrazsd * quotad(ji,jj,jk)

         ! Grazing by microzooplankton
         zgrazing(ji,jj,jk) = zgraztotc

         !    Various remineralization and excretion terms
         !    --------------------------------------------
         zgrasrat  = ( zgraztotf + rtrn ) / ( zgraztotc + rtrn )
         zgrasratn = ( zgraztotn + rtrn ) / ( zgraztotc + rtrn )
         zepshert  =  MIN( 1., zgrasratn, zgrasrat / ferat3)
         zbeta     = MAX(0., (epsher - epshermin) )
         zepsherf  = epshermin + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
         zepsherv  = zepsherf * zepshert 

         zgrafer   = zgraztotc * MAX( 0. , ( 1. - unass ) * zgrasrat - ferat3 * zepsherv ) 
         zgrarem   = zgraztotc * ( 1. - zepsherv - unass )
         zgrapoc   = zgraztotc * unass

         !  Update of the TRA arrays
         !  ------------------------
         zgrarsig  = zgrarem * sigma1
         tr(ji,jj,jk,jppo4,Krhs) = tr(ji,jj,jk,jppo4,Krhs) + zgrarsig
         tr(ji,jj,jk,jpnh4,Krhs) = tr(ji,jj,jk,jpnh4,Krhs) + zgrarsig
         tr(ji,jj,jk,jpdoc,Krhs) = tr(ji,jj,jk,jpdoc,Krhs) + zgrarem - zgrarsig
         !
         IF( ln_ligand ) THEN
            tr(ji,jj,jk,jplgw,Krhs) = tr(ji,jj,jk,jplgw,Krhs) + (zgrarem - zgrarsig) * ldocz
            zzligprod(ji,jj,jk) = (zgrarem - zgrarsig) * ldocz
         ENDIF
         !
         tr(ji,jj,jk,jpoxy,Krhs) = tr(ji,jj,jk,jpoxy,Krhs) - o2ut * zgrarsig
         tr(ji,jj,jk,jpfer,Krhs) = tr(ji,jj,jk,jpfer,Krhs) + zgrafer
         zfezoo(ji,jj,jk)    = zgrafer
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + zgrapoc
         prodpoc(ji,jj,jk)   = prodpoc(ji,jj,jk) + zgrapoc
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + zgraztotf * unass
         tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) + zgrarsig
         tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) + rno3 * zgrarsig
         !   Update the arrays TRA which contain the biological sources and sinks
         !   --------------------------------------------------------------------
         zmortz = ztortz + zrespz
         tr(ji,jj,jk,jpzoo,Krhs) = tr(ji,jj,jk,jpzoo,Krhs) - zmortz + zepsherv * zgraztotc 
         tr(ji,jj,jk,jpphy,Krhs) = tr(ji,jj,jk,jpphy,Krhs) - zgrazp
         tr(ji,jj,jk,jpdia,Krhs) = tr(ji,jj,jk,jpdia,Krhs) - zgrazsd
         tr(ji,jj,jk,jpnch,Krhs) = tr(ji,jj,jk,jpnch,Krhs) - zgrazp  * tr(ji,jj,jk,jpnch,Kbb)/(tr(ji,jj,jk,jpphy,Kbb)+rtrn)
         tr(ji,jj,jk,jpdch,Krhs) = tr(ji,jj,jk,jpdch,Krhs) - zgrazsd * tr(ji,jj,jk,jpdch,Kbb)/(tr(ji,jj,jk,jpdia,Kbb)+rtrn)
         tr(ji,jj,jk,jpdsi,Krhs) = tr(ji,jj,jk,jpdsi,Krhs) - zgrazsd * tr(ji,jj,jk,jpdsi,Kbb)/(tr(ji,jj,jk,jpdia,Kbb)+rtrn)
         tr(ji,jj,jk,jpgsi,Krhs) = tr(ji,jj,jk,jpgsi,Krhs) + zgrazsd * tr(ji,jj,jk,jpdsi,Kbb)/(tr(ji,jj,jk,jpdia,Kbb)+rtrn)
         tr(ji,jj,jk,jpnfe,Krhs) = tr(ji,jj,jk,jpnfe,Krhs) - zgrazpf
         tr(ji,jj,jk,jpdfe,Krhs) = tr(ji,jj,jk,jpdfe,Krhs) - zgrazsf
         tr(ji,jj,jk,jppoc,Krhs) = tr(ji,jj,jk,jppoc,Krhs) + zmortz - zgrazm
         prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortz
         conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazm
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + ferat3 * zmortz - zgrazmf
         !
         ! calcite production
         zprcaca = xfracal(ji,jj,jk) * zgrazp
         prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         !
         zprcaca = part * zprcaca
         tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) - zprcaca
         tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) - 2. * zprcaca
         tr(ji,jj,jk,jpcal,Krhs) = tr(ji,jj,jk,jpcal,Krhs) + zprcaca
      END_3D
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
        IF( iom_use("GRAZ1") ) THEN  !   Total grazing of phyto by zooplankton
           zgrazing(:,:,jpk) = 0._wp   ; CALL iom_put( "GRAZ1" , zgrazing(:,:,:) * 1.e+3  * rfact2r * tmask(:,:,:) ) 
         ENDIF
         IF( iom_use("FEZOO") ) THEN  
           zfezoo (:,:,jpk) = 0._wp    ; CALL iom_put( "FEZOO", zfezoo(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:) )
         ENDIF
         IF( ln_ligand ) THEN
            zzligprod(:,:,jpk) = 0._wp ; CALL iom_put( "LPRODZ", zzligprod(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:))
         ENDIF
      ENDIF
      !
      IF(sn_cfctl%l_prttrc) THEN      ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_micro')
      !
   END SUBROUTINE p4z_micro


   SUBROUTINE p4z_micro_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiszoo
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zzoo/ part, grazrat, resrat, mzrat, xprefn, xprefc, &
         &                xprefd,  xthreshdia,  xthreshphy,  xthreshpoc, &
         &                xthresh, xkgraz, epsher, epshermin, sigma1, unass
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) 
         WRITE(numout,*) 'p4z_micro_init : Initialization of microzooplankton parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF
      !
      READ  ( numnatp_ref, namp4zzoo, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zzoo in reference namelist' )
      READ  ( numnatp_cfg, namp4zzoo, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zzoo in configuration namelist' )
      IF(lwm) WRITE( numonp, namp4zzoo )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zzoo'
         WRITE(numout,*) '      part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(numout,*) '      microzoo preference for POC                     xprefc      =', xprefc
         WRITE(numout,*) '      microzoo preference for nano                    xprefn      =', xprefn
         WRITE(numout,*) '      microzoo preference for diatoms                 xprefd      =', xprefd
         WRITE(numout,*) '      diatoms feeding threshold  for microzoo         xthreshdia  =', xthreshdia
         WRITE(numout,*) '      nanophyto feeding threshold for microzoo        xthreshphy  =', xthreshphy
         WRITE(numout,*) '      poc feeding threshold for microzoo              xthreshpoc  =', xthreshpoc
         WRITE(numout,*) '      feeding threshold for microzooplankton          xthresh     =', xthresh
         WRITE(numout,*) '      exsudation rate of microzooplankton             resrat      =', resrat
         WRITE(numout,*) '      microzooplankton mortality rate                 mzrat       =', mzrat
         WRITE(numout,*) '      maximal microzoo grazing rate                   grazrat     =', grazrat
         WRITE(numout,*) '      non assimilated fraction of P by microzoo       unass       =', unass
         WRITE(numout,*) '      Efficicency of microzoo growth                  epsher      =', epsher
         WRITE(numout,*) '      Minimum efficicency of microzoo growth          epshermin   =', epshermin
         WRITE(numout,*) '      Fraction of microzoo excretion as DOM           sigma1      =', sigma1
         WRITE(numout,*) '      half sturation constant for grazing 1           xkgraz      =', xkgraz
      ENDIF
      !
   END SUBROUTINE p4z_micro_init

   !!======================================================================
END MODULE p4zmicro
