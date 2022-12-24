MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/dissolution of organic compounds
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/dissolution of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zche          !  chemical model
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zlim
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_rem_alloc

   REAL(wp), PUBLIC ::   xremikc    !: remineralisation rate of DOC 
   REAL(wp), PUBLIC ::   xremikn    !: remineralisation rate of DON 
   REAL(wp), PUBLIC ::   xremikp    !: remineralisation rate of DOP 
   REAL(wp), PUBLIC ::   xremik     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::   nitrif     !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::   xsirem     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::   xsiremlab  !: fast remineralisation rate of POC 
   REAL(wp), PUBLIC ::   xsilab     !: fraction of labile biogenic silica 
   REAL(wp), PUBLIC ::   feratb     !: Fe/C quota in bacteria
   REAL(wp), PUBLIC ::   xkferb     !: Half-saturation constant for bacteria Fe/C

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr   !: denitrification array

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zrem.F90 12377 2020-02-12 14:39:06Z acc $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt         ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremik, zremikc, zremikn, zremikp, zsiremin, zfact 
      REAL(wp) ::   zsatur, zsatur2, znusil, znusil2, zdep, zdepmin, zfactdep
      REAL(wp) ::   zbactfer, zolimit, zonitr, zrfact2
      REAL(wp) ::   zammonic, zoxyremc, zoxyremn, zoxyremp
      REAL(wp) ::   zosil, ztem, zdenitnh4, zolimic, zolimin, zolimip, zdenitrn, zdenitrp
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj    ) :: ztempbac
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zdepbac, zolimi, zdepprod, zfacsi, zfacsib, zdepeff, zfebact
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_rem')
      !
      ! Initialisation of arrys
      zdepprod(:,:,:) = 1._wp
      zdepeff (:,:,:) = 0.3_wp
      ztempbac(:,:)   = 0._wp
      zfacsib(:,:,:)  = xsilab / ( 1.0 - xsilab )
      zfebact(:,:,:)  = 0._wp
      zfacsi(:,:,:)   = xsilab

      ! Computation of the mean phytoplankton concentration as
      ! a crude estimate of the bacterial biomass
      ! this parameterization has been deduced from a model version
      ! that was modeling explicitely bacteria
      ! -------------------------------------------------------
      DO_3D_11_11( 1, jpkm1 )
         zdep = MAX( hmld(ji,jj), heup(ji,jj) )
         IF( gdept(ji,jj,jk,Kmm) < zdep ) THEN
            zdepbac(ji,jj,jk) = MIN( 0.7 * ( tr(ji,jj,jk,jpzoo,Kbb) + 2.* tr(ji,jj,jk,jpmes,Kbb) ), 4.e-6 )
            ztempbac(ji,jj)   = zdepbac(ji,jj,jk)
         ELSE
            zdepmin = MIN( 1., zdep / gdept(ji,jj,jk,Kmm) )
            zdepbac (ji,jj,jk) = zdepmin**0.683 * ztempbac(ji,jj)
            zdepprod(ji,jj,jk) = zdepmin**0.273
            zdepeff (ji,jj,jk) = zdepeff(ji,jj,jk) * zdepmin**0.3
         ENDIF
      END_3D

      IF( ln_p4z ) THEN
         DO_3D_11_11( 1, jpkm1 )
            ! DOC ammonification. Depends on depth, phytoplankton biomass
            ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
            zremik = xremik * xstep / 1.e-6 * xlimbac(ji,jj,jk) * zdepbac(ji,jj,jk) 
            zremik = MAX( zremik, 2.74e-4 * xstep )
            ! Ammonification in oxic waters with oxygen consumption
            ! -----------------------------------------------------
            zolimit = zremik * ( 1.- nitrfac(ji,jj,jk) ) * tr(ji,jj,jk,jpdoc,Kbb) 
            zolimi(ji,jj,jk) = MIN( ( tr(ji,jj,jk,jpoxy,Kbb) - rtrn ) / o2ut, zolimit ) 
            ! Ammonification in suboxic waters with denitrification
            ! -------------------------------------------------------
            zammonic = zremik * nitrfac(ji,jj,jk) * tr(ji,jj,jk,jpdoc,Kbb)
            denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) )
            denitr(ji,jj,jk)  = MIN( ( tr(ji,jj,jk,jpno3,Kbb) - rtrn ) / rdenit, denitr(ji,jj,jk) )
            zoxyremc          = zammonic - denitr(ji,jj,jk)
            !
            zolimi (ji,jj,jk) = MAX( 0.e0, zolimi (ji,jj,jk) )
            denitr (ji,jj,jk) = MAX( 0.e0, denitr (ji,jj,jk) )
            zoxyremc          = MAX( 0.e0, zoxyremc )

            !
            tr(ji,jj,jk,jppo4,Krhs) = tr(ji,jj,jk,jppo4,Krhs) + zolimi (ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc
            tr(ji,jj,jk,jpnh4,Krhs) = tr(ji,jj,jk,jpnh4,Krhs) + zolimi (ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc
            tr(ji,jj,jk,jpno3,Krhs) = tr(ji,jj,jk,jpno3,Krhs) - denitr (ji,jj,jk) * rdenit
            tr(ji,jj,jk,jpdoc,Krhs) = tr(ji,jj,jk,jpdoc,Krhs) - zolimi (ji,jj,jk) - denitr(ji,jj,jk) - zoxyremc
            tr(ji,jj,jk,jpoxy,Krhs) = tr(ji,jj,jk,jpoxy,Krhs) - zolimi (ji,jj,jk) * o2ut
            tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) + zolimi (ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc
            tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) + rno3 * ( zolimi(ji,jj,jk) + zoxyremc    &
            &                     + ( rdenit + 1.) * denitr(ji,jj,jk) )
         END_3D
      ELSE
         DO_3D_11_11( 1, jpkm1 )
            ! DOC ammonification. Depends on depth, phytoplankton biomass
            ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
            ! -----------------------------------------------------------------
            zremik = xstep / 1.e-6 * MAX(0.01, xlimbac(ji,jj,jk)) * zdepbac(ji,jj,jk) 
            zremik = MAX( zremik, 2.74e-4 * xstep / xremikc )

            zremikc = xremikc * zremik
            zremikn = xremikn / xremikc
            zremikp = xremikp / xremikc

            ! Ammonification in oxic waters with oxygen consumption
            ! -----------------------------------------------------
            zolimit = zremikc * ( 1.- nitrfac(ji,jj,jk) ) * tr(ji,jj,jk,jpdoc,Kbb) 
            zolimic = MAX( 0.e0, MIN( ( tr(ji,jj,jk,jpoxy,Kbb) - rtrn ) / o2ut, zolimit ) ) 
            zolimi(ji,jj,jk) = zolimic
            zolimin = zremikn * zolimic * tr(ji,jj,jk,jpdon,Kbb) / ( tr(ji,jj,jk,jpdoc,Kbb) + rtrn )
            zolimip = zremikp * zolimic * tr(ji,jj,jk,jpdop,Kbb) / ( tr(ji,jj,jk,jpdoc,Kbb) + rtrn ) 

            ! Ammonification in suboxic waters with denitrification
            ! -------------------------------------------------------
            zammonic = zremikc * nitrfac(ji,jj,jk) * tr(ji,jj,jk,jpdoc,Kbb)
            denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) )
            denitr(ji,jj,jk)  = MAX(0., MIN(  ( tr(ji,jj,jk,jpno3,Kbb) - rtrn ) / rdenit, denitr(ji,jj,jk) ) )
            zoxyremc          = MAX(0., zammonic - denitr(ji,jj,jk))
            zdenitrn  = zremikn * denitr(ji,jj,jk) * tr(ji,jj,jk,jpdon,Kbb) / ( tr(ji,jj,jk,jpdoc,Kbb) + rtrn )
            zdenitrp  = zremikp * denitr(ji,jj,jk) * tr(ji,jj,jk,jpdop,Kbb) / ( tr(ji,jj,jk,jpdoc,Kbb) + rtrn )
            zoxyremn  = zremikn * zoxyremc * tr(ji,jj,jk,jpdon,Kbb) / ( tr(ji,jj,jk,jpdoc,Kbb) + rtrn )
            zoxyremp  = zremikp * zoxyremc * tr(ji,jj,jk,jpdop,Kbb) / ( tr(ji,jj,jk,jpdoc,Kbb) + rtrn )

            tr(ji,jj,jk,jppo4,Krhs) = tr(ji,jj,jk,jppo4,Krhs) + zolimip + zdenitrp + zoxyremp
            tr(ji,jj,jk,jpnh4,Krhs) = tr(ji,jj,jk,jpnh4,Krhs) + zolimin + zdenitrn + zoxyremn
            tr(ji,jj,jk,jpno3,Krhs) = tr(ji,jj,jk,jpno3,Krhs) - denitr(ji,jj,jk) * rdenit
            tr(ji,jj,jk,jpdoc,Krhs) = tr(ji,jj,jk,jpdoc,Krhs) - zolimic - denitr(ji,jj,jk) - zoxyremc
            tr(ji,jj,jk,jpdon,Krhs) = tr(ji,jj,jk,jpdon,Krhs) - zolimin - zdenitrn - zoxyremn
            tr(ji,jj,jk,jpdop,Krhs) = tr(ji,jj,jk,jpdop,Krhs) - zolimip - zdenitrp - zoxyremp
            tr(ji,jj,jk,jpoxy,Krhs) = tr(ji,jj,jk,jpoxy,Krhs) - zolimic * o2ut
            tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) + zolimic + denitr(ji,jj,jk) + zoxyremc
            tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) + rno3 * ( zolimin + zoxyremn + ( rdenit + 1.) * zdenitrn )
         END_3D
         !
      ENDIF


      DO_3D_11_11( 1, jpkm1 )
         ! NH4 nitrification to NO3. Ceased for oxygen concentrations
         ! below 2 umol/L. Inhibited at strong light 
         ! ----------------------------------------------------------
         zonitr  = nitrif * xstep * tr(ji,jj,jk,jpnh4,Kbb) * ( 1.- nitrfac(ji,jj,jk) )  &
         &         / ( 1.+ emoy(ji,jj,jk) ) * ( 1. + fr_i(ji,jj) * emoy(ji,jj,jk) ) 
         zdenitnh4 = nitrif * xstep * tr(ji,jj,jk,jpnh4,Kbb) * nitrfac(ji,jj,jk)
         zdenitnh4 = MIN(  ( tr(ji,jj,jk,jpno3,Kbb) - rtrn ) / rdenita, zdenitnh4 ) 
         ! Update of the tracers trends
         ! ----------------------------
         tr(ji,jj,jk,jpnh4,Krhs) = tr(ji,jj,jk,jpnh4,Krhs) - zonitr - zdenitnh4
         tr(ji,jj,jk,jpno3,Krhs) = tr(ji,jj,jk,jpno3,Krhs) + zonitr - rdenita * zdenitnh4
         tr(ji,jj,jk,jpoxy,Krhs) = tr(ji,jj,jk,jpoxy,Krhs) - o2nit * zonitr
         tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) - 2 * rno3 * zonitr + rno3 * ( rdenita - 1. ) * zdenitnh4
      END_3D

       IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO_3D_11_11( 1, jpkm1 )

         ! Bacterial uptake of iron. No iron is available in DOC. So
         ! Bacteries are obliged to take up iron from the water. Some
         ! studies (especially at Papa) have shown this uptake to be significant
         ! ----------------------------------------------------------
         zbactfer = feratb *  rfact2 * 0.6_wp / rday * tgfunc(ji,jj,jk) * xlimbacl(ji,jj,jk)     &
            &              * tr(ji,jj,jk,jpfer,Kbb) / ( xkferb + tr(ji,jj,jk,jpfer,Kbb) )    &
            &              * zdepprod(ji,jj,jk) * zdepeff(ji,jj,jk) * zdepbac(ji,jj,jk)
         tr(ji,jj,jk,jpfer,Krhs) = tr(ji,jj,jk,jpfer,Krhs) - zbactfer*0.33
         tr(ji,jj,jk,jpsfe,Krhs) = tr(ji,jj,jk,jpsfe,Krhs) + zbactfer*0.25
         tr(ji,jj,jk,jpbfe,Krhs) = tr(ji,jj,jk,jpbfe,Krhs) + zbactfer*0.08
         zfebact(ji,jj,jk)   = zbactfer * 0.33
         blim(ji,jj,jk)      = xlimbacl(ji,jj,jk)  * zdepbac(ji,jj,jk) / 1.e-6 * zdepprod(ji,jj,jk)
      END_3D

       IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
       ENDIF

      ! Initialization of the array which contains the labile fraction
      ! of bSi. Set to a constant in the upper ocean
      ! ---------------------------------------------------------------

      DO_3D_11_11( 1, jpkm1 )
         zdep     = MAX( hmld(ji,jj), heup_01(ji,jj) )
         zsatur   = MAX( rtrn, ( sio3eq(ji,jj,jk) - tr(ji,jj,jk,jpsil,Kbb) ) / ( sio3eq(ji,jj,jk) + rtrn ) )
         zsatur2  = ( 1. + ts(ji,jj,jk,jp_tem,Kmm) / 400.)**37
         znusil   = 0.225  * ( 1. + ts(ji,jj,jk,jp_tem,Kmm) / 15.) * zsatur + 0.775 * zsatur2 * zsatur**9.25
         ! Remineralization rate of BSi depedant on T and saturation
         ! ---------------------------------------------------------
         IF ( gdept(ji,jj,jk,Kmm) > zdep ) THEN
            zfacsib(ji,jj,jk) = zfacsib(ji,jj,jk-1) * EXP( -0.5 * ( xsiremlab - xsirem )  &
            &                   * znusil * e3t(ji,jj,jk,Kmm) / wsbio4(ji,jj,jk) )
            zfacsi(ji,jj,jk)  = zfacsib(ji,jj,jk) / ( 1.0 + zfacsib(ji,jj,jk) )
            zfacsib(ji,jj,jk) = zfacsib(ji,jj,jk) * EXP( -0.5 * ( xsiremlab - xsirem )    &
            &                   * znusil * e3t(ji,jj,jk,Kmm) / wsbio4(ji,jj,jk) )
         ENDIF
         zsiremin = ( xsiremlab * zfacsi(ji,jj,jk) + xsirem * ( 1. - zfacsi(ji,jj,jk) ) ) * xstep * znusil
         zosil    = zsiremin * tr(ji,jj,jk,jpgsi,Kbb)
         !
         tr(ji,jj,jk,jpgsi,Krhs) = tr(ji,jj,jk,jpgsi,Krhs) - zosil
         tr(ji,jj,jk,jpsil,Krhs) = tr(ji,jj,jk,jpsil,Krhs) + zosil
      END_3D

      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tr(:,:,:,:,Krhs), mask=tmask, clinfo=ctrcnm)
       ENDIF

      IF( knt == nrdttrc ) THEN
          zrfact2 = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "REMIN" ) )  THEN !  Remineralisation rate
             zolimi(:,:,jpk) = 0. ; CALL iom_put( "REMIN"  , zolimi(:,:,:) * tmask(:,:,:) * zrfact2  )
          ENDIF
          CALL iom_put( "DENIT"  , denitr(:,:,:) * rdenit * rno3 * tmask(:,:,:) * zrfact2 ) ! Denitrification 
          IF( iom_use( "BACT" ) )  THEN ! Bacterial biomass
             zdepbac(:,:,jpk) = 0.  ;   CALL iom_put( "BACT", zdepbac(:,:,:) * 1.E6 * tmask(:,:,:) )
          ENDIF
          CALL iom_put( "FEBACT" , zfebact(:,:,:) * 1E9 * tmask(:,:,:) * zrfact2  )
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_rem')
      !
   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisrem/ xremik, nitrif, xsirem, xsiremlab, xsilab, feratb, xkferb, & 
         &                xremikc, xremikn, xremikp
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_rem_init : Initialization of remineralization parameters'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      READ  ( numnatp_ref, nampisrem, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisrem in reference namelist' )
      READ  ( numnatp_cfg, nampisrem, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisrem in configuration namelist' )
      IF(lwm) WRITE( numonp, nampisrem )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist parameters for remineralization, nampisrem'
         IF( ln_p4z ) THEN
            WRITE(numout,*) '      remineralization rate of DOC              xremik    =', xremik
         ELSE
            WRITE(numout,*) '      remineralization rate of DOC              xremikc   =', xremikc
            WRITE(numout,*) '      remineralization rate of DON              xremikn   =', xremikn
            WRITE(numout,*) '      remineralization rate of DOP              xremikp   =', xremikp
         ENDIF
         WRITE(numout,*) '      remineralization rate of Si               xsirem    =', xsirem
         WRITE(numout,*) '      fast remineralization rate of Si          xsiremlab =', xsiremlab
         WRITE(numout,*) '      fraction of labile biogenic silica        xsilab    =', xsilab
         WRITE(numout,*) '      NH4 nitrification rate                    nitrif    =', nitrif
         WRITE(numout,*) '      Bacterial Fe/C ratio                      feratb    =', feratb
         WRITE(numout,*) '      Half-saturation constant for bact. Fe/C   xkferb    =', xkferb
      ENDIF
      !
      denitr(:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_rem_init


   INTEGER FUNCTION p4z_rem_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( denitr(jpi,jpj,jpk), STAT=p4z_rem_alloc )
      !
      IF( p4z_rem_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_rem_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_rem_alloc

   !!======================================================================
END MODULE p4zrem
