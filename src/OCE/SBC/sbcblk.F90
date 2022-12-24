MODULE sbcblk
   !!======================================================================
   !!                       ***  MODULE  sbcblk  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!                         Aerodynamic Bulk Formulas
   !!                        SUCCESSOR OF "sbcblk_core"
   !!=====================================================================
   !! History :  1.0  !  2004-08  (U. Schweckendiek)  Original CORE code
   !!            2.0  !  2005-04  (L. Brodeau, A.M. Treguier)  improved CORE bulk and its user interface
   !!            3.0  !  2006-06  (G. Madec)  sbc rewritting
   !!             -   !  2006-12  (L. Brodeau)  Original code for turb_core
   !!            3.2  !  2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.3  !  2010-10  (S. Masson)  add diurnal cycle
   !!            3.4  !  2011-11  (C. Harris)  Fill arrays required by CICE
   !!            3.7  !  2014-06  (L. Brodeau)  simplification and optimization of CORE bulk
   !!            4.0  !  2016-06  (L. Brodeau)  sbcblk_core becomes sbcblk and is not restricted to the CORE algorithm anymore
   !!                 !                        ==> based on AeroBulk (https://github.com/brodeau/aerobulk/)
   !!            4.0  !  2016-10  (G. Madec)  introduce a sbc_blk_init routine
   !!            4.0  !  2016-10  (M. Vancoppenolle)  Introduce conduction flux emulator (M. Vancoppenolle)
   !!            4.0  !  2019-03  (F. Lemarié & G. Samson)  add ABL compatibility (ln_abl=TRUE)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_blk_init  : initialisation of the chosen bulk formulation as ocean surface boundary condition
   !!   sbc_blk       : bulk formulation as ocean surface boundary condition
   !!   blk_oce_1     : computes pieces of momentum, heat and freshwater fluxes over ocean for ABL model  (ln_abl=TRUE)
   !!   blk_oce_2     : finalizes momentum, heat and freshwater fluxes computation over ocean after the ABL step  (ln_abl=TRUE)
   !!             sea-ice case only :
   !!   blk_ice_1   : provide the air-ice stress
   !!   blk_ice_2   : provide the heat and mass fluxes at air-ice interface
   !!   blk_ice_qcn   : provide ice surface temperature and snow/ice conduction flux (emulating conduction flux)
   !!   Cdn10_Lupkes2012 : Lupkes et al. (2012) air-ice drag
   !!   Cdn10_Lupkes2015 : Lupkes et al. (2015) air-ice drag
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE fldread        ! read input fields
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE cyclone        ! Cyclone 10m wind form trac of cyclone centres
   USE sbcdcy         ! surface boundary condition: diurnal cycle
   USE sbcwave , ONLY :   cdn_wave ! wave module
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE lib_fortran    ! to use key_nosignedzero
#if defined key_si3
   USE ice     , ONLY :   jpl, a_i_b, at_i_b, rn_cnd_s, hfx_err_dif
   USE icethd_dh      ! for CALL ice_thd_snwblow
#endif
   USE sbcblk_algo_ncar     ! => turb_ncar     : NCAR - CORE (Large & Yeager, 2009)
   USE sbcblk_algo_coare3p0 ! => turb_coare3p0 : COAREv3.0 (Fairall et al. 2003)
   USE sbcblk_algo_coare3p6 ! => turb_coare3p6 : COAREv3.6 (Fairall et al. 2018 + Edson et al. 2013)
   USE sbcblk_algo_ecmwf    ! => turb_ecmwf    : ECMWF (IFS cycle 45r1)
   !
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control

   USE sbcblk_phy     ! a catalog of functions for physical/meteorological parameters in the marine boundary layer, rho_air, q_sat, etc...


   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_blk_init  ! called in sbcmod
   PUBLIC   sbc_blk       ! called in sbcmod
   PUBLIC   blk_oce_1     ! called in sbcabl
   PUBLIC   blk_oce_2     ! called in sbcabl
#if defined key_si3
   PUBLIC   blk_ice_1     ! routine called in icesbc
   PUBLIC   blk_ice_2     ! routine called in icesbc
   PUBLIC   blk_ice_qcn   ! routine called in icesbc
#endif

   INTEGER , PUBLIC            ::   jpfld         ! maximum number of files to read
   INTEGER , PUBLIC, PARAMETER ::   jp_wndi = 1   ! index of 10m wind velocity (i-component) (m/s)    at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_wndj = 2   ! index of 10m wind velocity (j-component) (m/s)    at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_tair = 3   ! index of 10m air temperature             (Kelvin)
   INTEGER , PUBLIC, PARAMETER ::   jp_humi = 4   ! index of specific humidity               ( % )
   INTEGER , PUBLIC, PARAMETER ::   jp_qsr  = 5   ! index of solar heat                      (W/m2)
   INTEGER , PUBLIC, PARAMETER ::   jp_qlw  = 6   ! index of Long wave                       (W/m2)
   INTEGER , PUBLIC, PARAMETER ::   jp_prec = 7   ! index of total precipitation (rain+snow) (Kg/m2/s)
   INTEGER , PUBLIC, PARAMETER ::   jp_snow = 8   ! index of snow (solid prcipitation)       (kg/m2/s)
   INTEGER , PUBLIC, PARAMETER ::   jp_slp  = 9   ! index of sea level pressure              (Pa)
   INTEGER , PUBLIC, PARAMETER ::   jp_hpgi =10   ! index of ABL geostrophic wind or hpg (i-component) (m/s) at T-point
   INTEGER , PUBLIC, PARAMETER ::   jp_hpgj =11   ! index of ABL geostrophic wind or hpg (j-component) (m/s) at T-point

   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:) ::   sf   ! structure of input atmospheric fields (file informations, fields read)

   !                           !!* Namelist namsbc_blk : bulk parameters
   LOGICAL  ::   ln_NCAR        ! "NCAR"      algorithm   (Large and Yeager 2008)
   LOGICAL  ::   ln_COARE_3p0   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
   LOGICAL  ::   ln_COARE_3p6   ! "COARE 3.6" algorithm   (Edson et al. 2013)
   LOGICAL  ::   ln_ECMWF       ! "ECMWF"     algorithm   (IFS cycle 45r1)
   !
   LOGICAL  ::   ln_Cd_L12      ! ice-atm drag = F( ice concentration )                        (Lupkes et al. JGR2012)
   LOGICAL  ::   ln_Cd_L15      ! ice-atm drag = F( ice concentration, atmospheric stability ) (Lupkes et al. JGR2015)
   !
   REAL(wp)         ::   rn_pfac   ! multiplication factor for precipitation
   REAL(wp), PUBLIC ::   rn_efac   ! multiplication factor for evaporation
   REAL(wp), PUBLIC ::   rn_vfac   ! multiplication factor for ice/ocean velocity in the calculation of wind stress
   REAL(wp)         ::   rn_zqt    ! z(q,t) : height of humidity and temperature measurements
   REAL(wp)         ::   rn_zu     ! z(u)   : height of wind measurements
   !
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   Cd_ice , Ch_ice , Ce_ice   ! transfert coefficients over ice
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   Cdn_oce, Chn_oce, Cen_oce  ! neutral coeffs over ocean (L15 bulk scheme)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   t_zu, q_zu                 ! air temp. and spec. hum. at wind speed height (L15 bulk scheme)

   LOGICAL  ::   ln_skin_cs     ! use the cool-skin (only available in ECMWF and COARE algorithms) !LB
   LOGICAL  ::   ln_skin_wl     ! use the warm-layer parameterization (only available in ECMWF and COARE algorithms) !LB
   LOGICAL  ::   ln_humi_sph    ! humidity read in files ("sn_humi") is specific humidity [kg/kg] if .true. !LB
   LOGICAL  ::   ln_humi_dpt    ! humidity read in files ("sn_humi") is dew-point temperature [K] if .true. !LB
   LOGICAL  ::   ln_humi_rlh    ! humidity read in files ("sn_humi") is relative humidity     [%] if .true. !LB
   !
   INTEGER  ::   nhumi          ! choice of the bulk algorithm
   !                            ! associated indices:
   INTEGER, PARAMETER :: np_humi_sph = 1
   INTEGER, PARAMETER :: np_humi_dpt = 2
   INTEGER, PARAMETER :: np_humi_rlh = 3

   INTEGER  ::   nblk           ! choice of the bulk algorithm
   !                            ! associated indices:
   INTEGER, PARAMETER ::   np_NCAR      = 1   ! "NCAR" algorithm        (Large and Yeager 2008)
   INTEGER, PARAMETER ::   np_COARE_3p0 = 2   ! "COARE 3.0" algorithm   (Fairall et al. 2003)
   INTEGER, PARAMETER ::   np_COARE_3p6 = 3   ! "COARE 3.6" algorithm   (Edson et al. 2013)
   INTEGER, PARAMETER ::   np_ECMWF     = 4   ! "ECMWF" algorithm       (IFS cycle 45r1)

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcblk.F90 12489 2020-02-28 15:55:11Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_blk_alloc()
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE sbc_blk_alloc ***
      !!-------------------------------------------------------------------
      ALLOCATE( t_zu(jpi,jpj)   , q_zu(jpi,jpj)   ,                                      &
         &      Cdn_oce(jpi,jpj), Chn_oce(jpi,jpj), Cen_oce(jpi,jpj),                    &
         &      Cd_ice (jpi,jpj), Ch_ice (jpi,jpj), Ce_ice (jpi,jpj), STAT=sbc_blk_alloc )
      !
      CALL mpp_sum ( 'sbcblk', sbc_blk_alloc )
      IF( sbc_blk_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_alloc: failed to allocate arrays' )
   END FUNCTION sbc_blk_alloc


   SUBROUTINE sbc_blk_init
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk_init  ***
      !!
      !! ** Purpose :   choose and initialize a bulk formulae formulation
      !!
      !! ** Method  :
      !!
      !!----------------------------------------------------------------------
      INTEGER  ::   jfpr                  ! dummy loop indice and argument
      INTEGER  ::   ios, ierror, ioptio   ! Local integer
      !!
      CHARACTER(len=100)            ::   cn_dir                ! Root directory for location of atmospheric forcing files
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::   slf_i        ! array of namelist informations on the fields to read
      TYPE(FLD_N) ::   sn_wndi, sn_wndj, sn_humi, sn_qsr       ! informations about the fields to be read
      TYPE(FLD_N) ::   sn_qlw , sn_tair, sn_prec, sn_snow      !       "                        "
      TYPE(FLD_N) ::   sn_slp , sn_hpgi, sn_hpgj               !       "                        "
      NAMELIST/namsbc_blk/ sn_wndi, sn_wndj, sn_humi, sn_qsr, sn_qlw ,                &   ! input fields
         &                 sn_tair, sn_prec, sn_snow, sn_slp, sn_hpgi, sn_hpgj,       &
         &                 ln_NCAR, ln_COARE_3p0, ln_COARE_3p6, ln_ECMWF,             &   ! bulk algorithm
         &                 cn_dir , rn_zqt, rn_zu,                                    &
         &                 rn_pfac, rn_efac, rn_vfac, ln_Cd_L12, ln_Cd_L15,           &
         &                 ln_skin_cs, ln_skin_wl, ln_humi_sph, ln_humi_dpt, ln_humi_rlh  ! cool-skin / warm-layer !LB
      !!---------------------------------------------------------------------
      !
      !                                      ! allocate sbc_blk_core array
      IF( sbc_blk_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_blk : unable to allocate standard arrays' )
      !
      !                             !** read bulk namelist
      READ  ( numnam_ref, namsbc_blk, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_blk in reference namelist' )
      !
      READ  ( numnam_cfg, namsbc_blk, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_blk in configuration namelist' )
      !
      IF(lwm) WRITE( numond, namsbc_blk )
      !
      !                             !** initialization of the chosen bulk formulae (+ check)
      !                                   !* select the bulk chosen in the namelist and check the choice
      ioptio = 0
      IF( ln_NCAR      ) THEN
         nblk =  np_NCAR        ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_COARE_3p0 ) THEN
         nblk =  np_COARE_3p0   ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_COARE_3p6 ) THEN
         nblk =  np_COARE_3p6   ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_ECMWF     ) THEN
         nblk =  np_ECMWF       ;   ioptio = ioptio + 1
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'sbc_blk_init: Choose one and only one bulk algorithm' )

      !                             !** initialization of the cool-skin / warm-layer parametrization
      IF( ln_skin_cs .OR. ln_skin_wl ) THEN
         !! Some namelist sanity tests:
         IF( ln_NCAR )      &
            & CALL ctl_stop( 'sbc_blk_init: Cool-skin/warm-layer param. not compatible with NCAR algorithm' )
         IF( nn_fsbc /= 1 ) &
            & CALL ctl_stop( 'sbc_blk_init: Please set "nn_fsbc" to 1 when using cool-skin/warm-layer param.')
      END IF

      IF( ln_skin_wl ) THEN
         !! Check if the frequency of downwelling solar flux input makes sense and if ln_dm2dc=T if it is daily!
         IF( (sn_qsr%freqh  < 0.).OR.(sn_qsr%freqh  > 24.) ) &
            & CALL ctl_stop( 'sbc_blk_init: Warm-layer param. (ln_skin_wl) not compatible with freq. of solar flux > daily' )
         IF( (sn_qsr%freqh == 24.).AND.(.NOT. ln_dm2dc) ) &
            & CALL ctl_stop( 'sbc_blk_init: Please set ln_dm2dc=T for warm-layer param. (ln_skin_wl) to work properly' )
      END IF

      ioptio = 0
      IF( ln_humi_sph ) THEN
         nhumi =  np_humi_sph    ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_humi_dpt ) THEN
         nhumi =  np_humi_dpt    ;   ioptio = ioptio + 1
      ENDIF
      IF( ln_humi_rlh ) THEN
         nhumi =  np_humi_rlh    ;   ioptio = ioptio + 1
      ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'sbc_blk_init: Choose one and only one type of air humidity' )
      !
      IF( ln_dm2dc ) THEN                 !* check: diurnal cycle on Qsr
         IF( sn_qsr%freqh /= 24. )   CALL ctl_stop( 'sbc_blk_init: ln_dm2dc=T only with daily short-wave input' )
         IF( sn_qsr%ln_tint ) THEN
            CALL ctl_warn( 'sbc_blk_init: ln_dm2dc=T daily qsr time interpolation done by sbcdcy module',   &
               &           '              ==> We force time interpolation = .false. for qsr' )
            sn_qsr%ln_tint = .false.
         ENDIF
      ENDIF
      !                                   !* set the bulk structure
      !                                      !- store namelist information in an array
      IF( ln_blk ) jpfld = 9
      IF( ln_abl ) jpfld = 11
      ALLOCATE( slf_i(jpfld) )
      !
      slf_i(jp_wndi) = sn_wndi   ;   slf_i(jp_wndj) = sn_wndj
      slf_i(jp_qsr ) = sn_qsr    ;   slf_i(jp_qlw ) = sn_qlw
      slf_i(jp_tair) = sn_tair   ;   slf_i(jp_humi) = sn_humi
      slf_i(jp_prec) = sn_prec   ;   slf_i(jp_snow) = sn_snow
      slf_i(jp_slp ) = sn_slp
      IF( ln_abl ) THEN
         slf_i(jp_hpgi) = sn_hpgi   ;   slf_i(jp_hpgj) = sn_hpgj
      END IF
      !
      !                                      !- allocate the bulk structure
      ALLOCATE( sf(jpfld), STAT=ierror )
      IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_blk_init: unable to allocate sf structure' )
      !
      !                                      !- fill the bulk structure with namelist informations
      CALL fld_fill( sf, slf_i, cn_dir, 'sbc_blk_init', 'surface boundary condition -- bulk formulae', 'namsbc_blk' )
      !
      DO jfpr= 1, jpfld
         !
         IF( TRIM(sf(jfpr)%clrootname) == 'NOT USED' ) THEN    !--  not used field  --!   (only now allocated and set to zero)
            ALLOCATE( sf(jfpr)%fnow(jpi,jpj,1) )
            sf(jfpr)%fnow(:,:,1) = 0._wp
         ELSE                                                  !-- used field  --!
            IF(   ln_abl    .AND.                                                      &
               &    ( jfpr == jp_wndi .OR. jfpr == jp_wndj .OR. jfpr == jp_humi .OR.   &
               &      jfpr == jp_hpgi .OR. jfpr == jp_hpgj .OR. jfpr == jp_tair     )  ) THEN   ! ABL: some fields are 3D input
               ALLOCATE( sf(jfpr)%fnow(jpi,jpj,jpka) )
               IF( sf(jfpr)%ln_tint )   ALLOCATE( sf(jfpr)%fdta(jpi,jpj,jpka,2) )
            ELSE                                                                                ! others or Bulk fields are 2D fiels
               ALLOCATE( sf(jfpr)%fnow(jpi,jpj,1) )
               IF( sf(jfpr)%ln_tint )   ALLOCATE( sf(jfpr)%fdta(jpi,jpj,1,2) )
            ENDIF
            !
            IF( sf(jfpr)%freqh > 0. .AND. MOD( NINT(3600. * sf(jfpr)%freqh), nn_fsbc * NINT(rn_Dt) ) /= 0 )   &
               &  CALL ctl_warn( 'sbc_blk_init: sbcmod timestep rn_Dt*nn_fsbc is NOT a submultiple of atmospheric forcing frequency.',   &
               &                 '               This is not ideal. You should consider changing either rn_Dt or nn_fsbc value...' )
         ENDIF
      END DO
      !
      IF( ln_wave ) THEN
         !Activated wave module but neither drag nor stokes drift activated
         IF( .NOT.(ln_cdgw .OR. ln_sdw .OR. ln_tauwoc .OR. ln_stcor ) )   THEN
            CALL ctl_stop( 'STOP',  'Ask for wave coupling but ln_cdgw=F, ln_sdw=F, ln_tauwoc=F, ln_stcor=F' )
            !drag coefficient read from wave model definable only with mfs bulk formulae and core
         ELSEIF(ln_cdgw .AND. .NOT. ln_NCAR )       THEN
            CALL ctl_stop( 'drag coefficient read from wave model definable only with NCAR and CORE bulk formulae')
         ELSEIF(ln_stcor .AND. .NOT. ln_sdw)                             THEN
            CALL ctl_stop( 'Stokes-Coriolis term calculated only if activated Stokes Drift ln_sdw=T')
         ENDIF
      ELSE
         IF( ln_cdgw .OR. ln_sdw .OR. ln_tauwoc .OR. ln_stcor )                &
            &   CALL ctl_stop( 'Not Activated Wave Module (ln_wave=F) but asked coupling ',    &
            &                  'with drag coefficient (ln_cdgw =T) '  ,                        &
            &                  'or Stokes Drift (ln_sdw=T) ' ,                                 &
            &                  'or ocean stress modification due to waves (ln_tauwoc=T) ',      &
            &                  'or Stokes-Coriolis term (ln_stcori=T)'  )
      ENDIF
      !
      IF( ln_abl ) THEN       ! ABL: read 3D fields for wind, temperature, humidity and pressure gradient
         rn_zqt = ght_abl(2)          ! set the bulk altitude to ABL first level
         rn_zu  = ght_abl(2)
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   ABL formulation: overwrite rn_zqt & rn_zu with ABL first level altitude'
      ENDIF
      !
      ! set transfer coefficients to default sea-ice values
      Cd_ice(:,:) = rCd_ice
      Ch_ice(:,:) = rCd_ice
      Ce_ice(:,:) = rCd_ice
      !
      IF(lwp) THEN                     !** Control print
         !
         WRITE(numout,*)                  !* namelist
         WRITE(numout,*) '   Namelist namsbc_blk (other than data information):'
         WRITE(numout,*) '      "NCAR"      algorithm   (Large and Yeager 2008)     ln_NCAR      = ', ln_NCAR
         WRITE(numout,*) '      "COARE 3.0" algorithm   (Fairall et al. 2003)       ln_COARE_3p0 = ', ln_COARE_3p0
         WRITE(numout,*) '      "COARE 3.6" algorithm (Fairall 2018 + Edson al 2013)ln_COARE_3p6 = ', ln_COARE_3p6
         WRITE(numout,*) '      "ECMWF"     algorithm   (IFS cycle 45r1)            ln_ECMWF     = ', ln_ECMWF
         WRITE(numout,*) '      Air temperature and humidity reference height (m)   rn_zqt       = ', rn_zqt
         WRITE(numout,*) '      Wind vector reference height (m)                    rn_zu        = ', rn_zu
         WRITE(numout,*) '      factor applied on precipitation (total & snow)      rn_pfac      = ', rn_pfac
         WRITE(numout,*) '      factor applied on evaporation                       rn_efac      = ', rn_efac
         WRITE(numout,*) '      factor applied on ocean/ice velocity                rn_vfac      = ', rn_vfac
         WRITE(numout,*) '         (form absolute (=0) to relative winds(=1))'
         WRITE(numout,*) '      use ice-atm drag from Lupkes2012                    ln_Cd_L12    = ', ln_Cd_L12
         WRITE(numout,*) '      use ice-atm drag from Lupkes2015                    ln_Cd_L15    = ', ln_Cd_L15
         !
         WRITE(numout,*)
         SELECT CASE( nblk )              !* Print the choice of bulk algorithm
         CASE( np_NCAR      )   ;   WRITE(numout,*) '   ==>>>   "NCAR" algorithm        (Large and Yeager 2008)'
         CASE( np_COARE_3p0 )   ;   WRITE(numout,*) '   ==>>>   "COARE 3.0" algorithm   (Fairall et al. 2003)'
         CASE( np_COARE_3p6 )   ;   WRITE(numout,*) '   ==>>>   "COARE 3.6" algorithm (Fairall 2018+Edson et al. 2013)'
         CASE( np_ECMWF     )   ;   WRITE(numout,*) '   ==>>>   "ECMWF" algorithm       (IFS cycle 45r1)'
         END SELECT
         !
         WRITE(numout,*)
         WRITE(numout,*) '      use cool-skin  parameterization (SSST)  ln_skin_cs  = ', ln_skin_cs
         WRITE(numout,*) '      use warm-layer parameterization (SSST)  ln_skin_wl  = ', ln_skin_wl
         !
         WRITE(numout,*)
         SELECT CASE( nhumi )              !* Print the choice of air humidity
         CASE( np_humi_sph )   ;   WRITE(numout,*) '   ==>>>   air humidity is SPECIFIC HUMIDITY     [kg/kg]'
         CASE( np_humi_dpt )   ;   WRITE(numout,*) '   ==>>>   air humidity is DEW-POINT TEMPERATURE [K]'
         CASE( np_humi_rlh )   ;   WRITE(numout,*) '   ==>>>   air humidity is RELATIVE HUMIDITY     [%]'
         END SELECT
         !
      ENDIF
      !
   END SUBROUTINE sbc_blk_init


   SUBROUTINE sbc_blk( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_blk  ***
      !!
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!              (momentum, heat, freshwater and runoff)
      !!
      !! ** Method  :
      !!              (1) READ each fluxes in NetCDF files:
      !!      the wind velocity (i-component) at z=rn_zu  (m/s) at T-point
      !!      the wind velocity (j-component) at z=rn_zu  (m/s) at T-point
      !!      the specific humidity           at z=rn_zqt (kg/kg)
      !!      the air temperature             at z=rn_zqt (Kelvin)
      !!      the solar heat                              (W/m2)
      !!      the Long wave                               (W/m2)
      !!      the total precipitation (rain+snow)         (Kg/m2/s)
      !!      the snow (solid precipitation)              (kg/m2/s)
      !!      ABL dynamical forcing (i/j-components of either hpg or geostrophic winds)
      !!              (2) CALL blk_oce_1 and blk_oce_2
      !!
      !!      C A U T I O N : never mask the surface stress fields
      !!                      the stress is assumed to be in the (i,j) mesh referential
      !!
      !! ** Action  :   defined at each time-step at the air-sea interface
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        wind speed  module at T-point over free ocean or leads in presence of sea-ice
      !!              - qns, qsr    non-solar and solar heat fluxes
      !!              - emp         upward mass flux (evapo. - precip.)
      !!              - sfx         salt flux due to freezing/melting (non-zero only if ice is present)
      !!
      !! ** References :   Large & Yeager, 2004 / Large & Yeager, 2008
      !!                   Brodeau et al. Ocean Modelling 2010
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) ::   zssq, zcd_du, zsen, zevp
      REAL(wp) :: ztmp
      !!----------------------------------------------------------------------
      !
      CALL fld_read( kt, nn_fsbc, sf )             ! input fields provided at the current time-step

      ! Sanity/consistence test on humidity at first time step to detect potential screw-up:
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*) ''
#if defined key_agrif
         IF(lwp) WRITE(numout,*) ' === AGRIF => Sanity/consistence test on air humidity SKIPPED! :( ==='
#else
         ztmp = SUM(tmask(:,:,1)) ! number of ocean points on local proc domain
         IF( ztmp > 8._wp ) THEN ! test only on proc domains with at least 8 ocean points!
            ztmp = SUM(sf(jp_humi)%fnow(:,:,1)*tmask(:,:,1))/ztmp ! mean humidity over ocean on proc
            SELECT CASE( nhumi )
            CASE( np_humi_sph ) ! specific humidity => expect: 0. <= something < 0.065 [kg/kg] (0.061 is saturation at 45degC !!!)
               IF(  (ztmp < 0._wp) .OR. (ztmp > 0.065)  ) ztmp = -1._wp
            CASE( np_humi_dpt ) ! dew-point temperature => expect: 110. <= something < 320. [K]
               IF( (ztmp < 110._wp).OR.(ztmp > 320._wp) ) ztmp = -1._wp
            CASE( np_humi_rlh ) ! relative humidity => expect: 0. <= something < 100. [%]
               IF(  (ztmp < 0._wp) .OR.(ztmp > 100._wp) ) ztmp = -1._wp
            END SELECT
            IF(ztmp < 0._wp) THEN
               IF (lwp) WRITE(numout,'("   Mean humidity value found on proc #",i6.6," is: ",f10.5)') narea, ztmp
               CALL ctl_stop( 'STOP', 'Something is wrong with air humidity!!!', &
                  &   ' ==> check the unit in your input files'       , &
                  &   ' ==> check consistence of namelist choice: specific? relative? dew-point?', &
                  &   ' ==> ln_humi_sph -> [kg/kg] | ln_humi_rlh -> [%] | ln_humi_dpt -> [K] !!!' )
            END IF
         END IF
         IF(lwp) WRITE(numout,*) ' === Sanity/consistence test on air humidity sucessfuly passed! ==='
#endif
         IF(lwp) WRITE(numout,*) ''
      END IF !IF( kt == nit000 )
      !                                            ! compute the surface ocean fluxes using bulk formulea
      IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN
         CALL blk_oce_1( kt, sf(jp_wndi)%fnow(:,:,1), sf(jp_wndj)%fnow(:,:,1),   &   !   <<= in
            &                sf(jp_tair)%fnow(:,:,1), sf(jp_humi)%fnow(:,:,1),   &   !   <<= in
            &                sf(jp_slp )%fnow(:,:,1), sst_m, ssu_m, ssv_m,       &   !   <<= in
            &                sf(jp_qsr )%fnow(:,:,1), sf(jp_qlw )%fnow(:,:,1),   &   !   <<= in (wl/cs)
            &                tsk_m, zssq, zcd_du, zsen, zevp )                       !   =>> out

         CALL blk_oce_2(     sf(jp_tair)%fnow(:,:,1), sf(jp_qsr )%fnow(:,:,1),   &   !   <<= in
            &                sf(jp_qlw )%fnow(:,:,1), sf(jp_prec)%fnow(:,:,1),   &   !   <<= in
            &                sf(jp_snow)%fnow(:,:,1), tsk_m,                     &   !   <<= in
            &                zsen, zevp )                                            !   <=> in out
      ENDIF
      !
#if defined key_cice
      IF( MOD( kt - 1, nn_fsbc ) == 0 )   THEN
         qlw_ice(:,:,1)   = sf(jp_qlw )%fnow(:,:,1)
         IF( ln_dm2dc ) THEN
            qsr_ice(:,:,1) = sbc_dcy( sf(jp_qsr)%fnow(:,:,1) )
         ELSE
            qsr_ice(:,:,1) =          sf(jp_qsr)%fnow(:,:,1)
         ENDIF
         tatm_ice(:,:)    = sf(jp_tair)%fnow(:,:,1)

         SELECT CASE( nhumi )
         CASE( np_humi_sph )
            qatm_ice(:,:) =           sf(jp_humi)%fnow(:,:,1)
         CASE( np_humi_dpt )
            qatm_ice(:,:) = q_sat(    sf(jp_humi)%fnow(:,:,1), sf(jp_slp)%fnow(:,:,1) )
         CASE( np_humi_rlh )
            qatm_ice(:,:) = q_air_rh( 0.01_wp*sf(jp_humi)%fnow(:,:,1), sf(jp_tair)%fnow(:,:,1), sf(jp_slp)%fnow(:,:,1)) !LB: 0.01 => RH is % percent in file
         END SELECT

         tprecip(:,:)     = sf(jp_prec)%fnow(:,:,1) * rn_pfac
         sprecip(:,:)     = sf(jp_snow)%fnow(:,:,1) * rn_pfac
         wndi_ice(:,:)    = sf(jp_wndi)%fnow(:,:,1)
         wndj_ice(:,:)    = sf(jp_wndj)%fnow(:,:,1)
      ENDIF
#endif
      !
   END SUBROUTINE sbc_blk


   SUBROUTINE blk_oce_1( kt, pwndi, pwndj , ptair, phumi, &  ! inp
      &                  pslp , pst   , pu   , pv,        &  ! inp
      &                  pqsr , pqlw  ,                   &  ! inp
      &                  ptsk, pssq , pcd_du, psen , pevp   )  ! out
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_oce_1  ***
      !!
      !! ** Purpose :   if ln_blk=T, computes surface momentum, heat and freshwater fluxes
      !!                if ln_abl=T, computes Cd x |U|, Ch x |U|, Ce x |U| for ABL integration
      !!
      !! ** Method  :   bulk formulae using atmospheric fields from :
      !!                if ln_blk=T, atmospheric fields read in sbc_read
      !!                if ln_abl=T, the ABL model at previous time-step
      !!
      !! ** Outputs : - pssq    : surface humidity used to compute latent heat flux (kg/kg)
      !!              - pcd_du  : Cd x |dU| at T-points  (m/s)
      !!              - psen    : Ch x |dU| at T-points  (m/s)
      !!              - pevp    : Ce x |dU| at T-points  (m/s)
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in   )                 ::   kt     ! time step index
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pwndi  ! atmospheric wind at U-point              [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pwndj  ! atmospheric wind at V-point              [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   phumi  ! specific humidity at T-points            [kg/kg]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   ptair  ! potential temperature at T-points        [Kelvin]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pslp   ! sea-level pressure                       [Pa]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pst    ! surface temperature                      [Celsius]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pu     ! surface current at U-point (i-component) [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pv     ! surface current at V-point (j-component) [m/s]
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pqsr   !
      REAL(wp), INTENT(in   ), DIMENSION(:,:) ::   pqlw   !
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   ptsk   ! skin temp. (or SST if CS & WL not used)  [Celsius]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pssq   ! specific humidity at pst                 [kg/kg]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pcd_du ! Cd x |dU| at T-points                    [m/s]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   psen   ! Ch x |dU| at T-points                    [m/s]
      REAL(wp), INTENT(  out), DIMENSION(:,:) ::   pevp   ! Ce x |dU| at T-points                    [m/s]
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zztmp                ! local variable
      REAL(wp), DIMENSION(jpi,jpj) ::   zwnd_i, zwnd_j    ! wind speed components at T-point
      REAL(wp), DIMENSION(jpi,jpj) ::   zU_zu             ! bulk wind speed at height zu  [m/s]
      REAL(wp), DIMENSION(jpi,jpj) ::   ztpot             ! potential temperature of air at z=rn_zqt [K]
      REAL(wp), DIMENSION(jpi,jpj) ::   zqair             ! specific humidity     of air at z=rn_zqt [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj) ::   zcd_oce           ! momentum transfert coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj) ::   zch_oce           ! sensible heat transfert coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj) ::   zce_oce           ! latent   heat transfert coefficient over ocean
      REAL(wp), DIMENSION(jpi,jpj) ::   zqla              ! latent heat flux
      REAL(wp), DIMENSION(jpi,jpj) ::   zztmp1, zztmp2
      !!---------------------------------------------------------------------
      !
      ! local scalars ( place there for vector optimisation purposes)
      !                           ! Temporary conversion from Celcius to Kelvin (and set minimum value far above 0 K)
      ptsk(:,:) = pst(:,:) + rt0  ! by default: skin temperature = "bulk SST" (will remain this way if NCAR algorithm used!)

      ! ----------------------------------------------------------------------------- !
      !      0   Wind components and module at T-point relative to the moving ocean   !
      ! ----------------------------------------------------------------------------- !

      ! ... components ( U10m - U_oce ) at T-point (unmasked)
#if defined key_cyclone
      zwnd_i(:,:) = 0._wp
      zwnd_j(:,:) = 0._wp
      CALL wnd_cyc( kt, zwnd_i, zwnd_j )    ! add analytical tropical cyclone (Vincent et al. JGR 2012)
      DO_2D_00_00
         pwndi(ji,jj) = pwndi(ji,jj) + zwnd_i(ji,jj)
         pwndj(ji,jj) = pwndj(ji,jj) + zwnd_j(ji,jj)
      END_2D
#endif
      DO_2D_00_00
         zwnd_i(ji,jj) = (  pwndi(ji,jj) - rn_vfac * 0.5 * ( pu(ji-1,jj  ) + pu(ji,jj) )  )
         zwnd_j(ji,jj) = (  pwndj(ji,jj) - rn_vfac * 0.5 * ( pv(ji  ,jj-1) + pv(ji,jj) )  )
      END_2D
      CALL lbc_lnk_multi( 'sbcblk', zwnd_i, 'T', -1., zwnd_j, 'T', -1. )
      ! ... scalar wind ( = | U10m - U_oce | ) at T-point (masked)
      wndm(:,:) = SQRT(  zwnd_i(:,:) * zwnd_i(:,:)   &
         &             + zwnd_j(:,:) * zwnd_j(:,:)  ) * tmask(:,:,1)

      ! ----------------------------------------------------------------------------- !
      !      I   Solar FLUX                                                           !
      ! ----------------------------------------------------------------------------- !

      ! ocean albedo assumed to be constant + modify now Qsr to include the diurnal cycle                    ! Short Wave
      zztmp = 1. - albo
      IF( ln_dm2dc ) THEN
         qsr(:,:) = zztmp * sbc_dcy( sf(jp_qsr)%fnow(:,:,1) ) * tmask(:,:,1)
      ELSE
         qsr(:,:) = zztmp *          sf(jp_qsr)%fnow(:,:,1)   * tmask(:,:,1)
      ENDIF


      ! ----------------------------------------------------------------------------- !
      !     II   Turbulent FLUXES                                                     !
      ! ----------------------------------------------------------------------------- !

      ! specific humidity at SST
      pssq(:,:) = rdct_qsat_salt * q_sat( ptsk(:,:), pslp(:,:) )

      IF( ln_skin_cs .OR. ln_skin_wl ) THEN
         !! Backup "bulk SST" and associated spec. hum.
         zztmp1(:,:) = ptsk(:,:)
         zztmp2(:,:) = pssq(:,:)
      ENDIF

      ! specific humidity of air at "rn_zqt" m above the sea
      SELECT CASE( nhumi )
      CASE( np_humi_sph )
         zqair(:,:) = phumi(:,:)      ! what we read in file is already a spec. humidity!
      CASE( np_humi_dpt )
         !IF(lwp) WRITE(numout,*) ' *** blk_oce => computing q_air out of d_air and slp !' !LBrm
         zqair(:,:) = q_sat( phumi(:,:), pslp(:,:) )
      CASE( np_humi_rlh )
         !IF(lwp) WRITE(numout,*) ' *** blk_oce => computing q_air out of RH, t_air and slp !' !LBrm
         zqair(:,:) = q_air_rh( 0.01_wp*phumi(:,:), ptair(:,:), pslp(:,:) ) !LB: 0.01 => RH is % percent in file
      END SELECT
      !
      ! potential temperature of air at "rn_zqt" m above the sea
      IF( ln_abl ) THEN
         ztpot = ptair(:,:)
      ELSE
         ! Estimate of potential temperature at z=rn_zqt, based on adiabatic lapse-rate
         !    (see Josey, Gulev & Yu, 2013) / doi=10.1016/B978-0-12-391851-2.00005-2
         !    (since reanalysis products provide T at z, not theta !)
         !#LB: because AGRIF hates functions that return something else than a scalar, need to
         !     use scalar version of gamma_moist() ...
         DO_2D_11_11
            ztpot(ji,jj) = ptair(ji,jj) + gamma_moist( ptair(ji,jj), zqair(ji,jj) ) * rn_zqt
         END_2D
      ENDIF



      !! Time to call the user-selected bulk parameterization for
      !!  ==  transfer coefficients  ==!   Cd, Ch, Ce at T-point, and more...
      SELECT CASE( nblk )

      CASE( np_NCAR      )
         CALL turb_ncar    ( rn_zqt, rn_zu, ptsk, ztpot, pssq, zqair, wndm,                              &
            &                zcd_oce, zch_oce, zce_oce, t_zu, q_zu, zU_zu, cdn_oce, chn_oce, cen_oce )

      CASE( np_COARE_3p0 )
         CALL turb_coare3p0 ( kt, rn_zqt, rn_zu, ptsk, ztpot, pssq, zqair, wndm, ln_skin_cs, ln_skin_wl, &
            &                zcd_oce, zch_oce, zce_oce, t_zu, q_zu, zU_zu, cdn_oce, chn_oce, cen_oce,   &
            &                Qsw=qsr(:,:), rad_lw=pqlw(:,:), slp=pslp(:,:) )

      CASE( np_COARE_3p6 )
         CALL turb_coare3p6 ( kt, rn_zqt, rn_zu, ptsk, ztpot, pssq, zqair, wndm, ln_skin_cs, ln_skin_wl, &
            &                zcd_oce, zch_oce, zce_oce, t_zu, q_zu, zU_zu, cdn_oce, chn_oce, cen_oce,   &
            &                Qsw=qsr(:,:), rad_lw=pqlw(:,:), slp=pslp(:,:) )

      CASE( np_ECMWF     )
         CALL turb_ecmwf   ( kt, rn_zqt, rn_zu, ptsk, ztpot, pssq, zqair, wndm, ln_skin_cs, ln_skin_wl,  &
            &                zcd_oce, zch_oce, zce_oce, t_zu, q_zu, zU_zu, cdn_oce, chn_oce, cen_oce,   &
            &                Qsw=qsr(:,:), rad_lw=pqlw(:,:), slp=pslp(:,:) )

      CASE DEFAULT
         CALL ctl_stop( 'STOP', 'sbc_oce: non-existing bulk formula selected' )

      END SELECT

      IF( ln_skin_cs .OR. ln_skin_wl ) THEN
         !! ptsk and pssq have been updated!!!
         !!
         !! In the presence of sea-ice we forget about the cool-skin/warm-layer update of ptsk and pssq:
         WHERE ( fr_i(:,:) > 0.001_wp )
            ! sea-ice present, we forget about the update, using what we backed up before call to turb_*()
            ptsk(:,:) = zztmp1(:,:)
            pssq(:,:) = zztmp2(:,:)
         END WHERE
      END IF

      !  Turbulent fluxes over ocean  => BULK_FORMULA @ sbcblk_phy.F90
      ! -------------------------------------------------------------

      IF( ln_abl ) THEN         !==  ABL formulation  ==!   multiplication by rho_air and turbulent fluxes computation done in ablstp
         !! FL do we need this multiplication by tmask ... ???
         DO_2D_11_11
            zztmp = zU_zu(ji,jj) !* tmask(ji,jj,1)
            wndm(ji,jj)   = zztmp                   ! Store zU_zu in wndm to compute ustar2 in ablmod
            pcd_du(ji,jj) = zztmp * zcd_oce(ji,jj)
            psen(ji,jj)   = zztmp * zch_oce(ji,jj)
            pevp(ji,jj)   = zztmp * zce_oce(ji,jj)
         END_2D
      ELSE                      !==  BLK formulation  ==!   turbulent fluxes computation
         CALL BULK_FORMULA( rn_zu, ptsk(:,:), pssq(:,:), t_zu(:,:), q_zu(:,:), &
            &               zcd_oce(:,:), zch_oce(:,:), zce_oce(:,:),          &
            &               wndm(:,:), zU_zu(:,:), pslp(:,:),                  &
            &               taum(:,:), psen(:,:), zqla(:,:),                   &
            &               pEvap=pevp(:,:), prhoa=rhoa(:,:), pfact_evap=rn_efac )

         zqla(:,:) = zqla(:,:) * tmask(:,:,1)
         psen(:,:) = psen(:,:) * tmask(:,:,1)
         taum(:,:) = taum(:,:) * tmask(:,:,1)
         pevp(:,:) = pevp(:,:) * tmask(:,:,1)

         ! Tau i and j component on T-grid points, using array "zcd_oce" as a temporary array...
         zcd_oce = 0._wp
         WHERE ( wndm > 0._wp ) zcd_oce = taum / wndm
         zwnd_i = zcd_oce * zwnd_i
         zwnd_j = zcd_oce * zwnd_j

         CALL iom_put( "taum_oce", taum )   ! output wind stress module

         ! ... utau, vtau at U- and V_points, resp.
         !     Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
         !     Note the use of MAX(tmask(i,j),tmask(i+1,j) is to mask tau over ice shelves
         DO_2D_10_10
            utau(ji,jj) = 0.5 * ( 2. - umask(ji,jj,1) ) * ( zwnd_i(ji,jj) + zwnd_i(ji+1,jj  ) ) &
               &          * MAX(tmask(ji,jj,1),tmask(ji+1,jj,1))
            vtau(ji,jj) = 0.5 * ( 2. - vmask(ji,jj,1) ) * ( zwnd_j(ji,jj) + zwnd_j(ji  ,jj+1) ) &
               &          * MAX(tmask(ji,jj,1),tmask(ji,jj+1,1))
         END_2D
         CALL lbc_lnk_multi( 'sbcblk', utau, 'U', -1., vtau, 'V', -1. )

         IF(sn_cfctl%l_prtctl) THEN
            CALL prt_ctl( tab2d_1=wndm  , clinfo1=' blk_oce_1: wndm   : ')
            CALL prt_ctl( tab2d_1=utau  , clinfo1=' blk_oce_1: utau   : ', mask1=umask,   &
               &          tab2d_2=vtau  , clinfo2='            vtau   : ', mask2=vmask )
         ENDIF
         !
      ENDIF !IF( ln_abl )
      
      ptsk(:,:) = ( ptsk(:,:) - rt0 ) * tmask(:,:,1)  ! Back to Celsius
            
      IF( ln_skin_cs .OR. ln_skin_wl ) THEN
         CALL iom_put( "t_skin" ,  ptsk        )  ! T_skin in Celsius
         CALL iom_put( "dt_skin" , ptsk - pst  )  ! T_skin - SST temperature difference...
      ENDIF

      IF(sn_cfctl%l_prtctl) THEN
         CALL prt_ctl( tab2d_1=pevp  , clinfo1=' blk_oce_1: pevp   : ' )
         CALL prt_ctl( tab2d_1=psen  , clinfo1=' blk_oce_1: psen   : ' )
         CALL prt_ctl( tab2d_1=pssq  , clinfo1=' blk_oce_1: pssq   : ' )
      ENDIF
      !
   END SUBROUTINE blk_oce_1


   SUBROUTINE blk_oce_2( ptair, pqsr, pqlw, pprec,   &   ! <<= in
      &                  psnow, ptsk, psen, pevp     )   ! <<= in
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_oce_2  ***
      !!
      !! ** Purpose :   finalize the momentum, heat and freshwater fluxes computation
      !!                at the ocean surface at each time step knowing Cd, Ch, Ce and
      !!                atmospheric variables (from ABL or external data)
      !!
      !! ** Outputs : - utau    : i-component of the stress at U-point  (N/m2)
      !!              - vtau    : j-component of the stress at V-point  (N/m2)
      !!              - taum    : Wind stress module at T-point         (N/m2)
      !!              - wndm    : Wind speed module at T-point          (m/s)
      !!              - qsr     : Solar heat flux over the ocean        (W/m2)
      !!              - qns     : Non Solar heat flux over the ocean    (W/m2)
      !!              - emp     : evaporation minus precipitation       (kg/m2/s)
      !!---------------------------------------------------------------------
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptair
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   pqsr
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   pqlw
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   pprec
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   psnow
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   ptsk   ! SKIN surface temperature   [Celsius]
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   psen
      REAL(wp), INTENT(in), DIMENSION(:,:) ::   pevp
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) ::   zztmp,zz1,zz2,zz3    ! local variable
      REAL(wp), DIMENSION(jpi,jpj) ::   ztskk             ! skin temp. in Kelvin 
      REAL(wp), DIMENSION(jpi,jpj) ::   zqlw              ! long wave and sensible heat fluxes      
      REAL(wp), DIMENSION(jpi,jpj) ::   zqla              ! latent heat fluxes and evaporation
      !!---------------------------------------------------------------------
      !
      ! local scalars ( place there for vector optimisation purposes)


      ztskk(:,:) = ptsk(:,:) + rt0  ! => ptsk in Kelvin rather than Celsius
      
      ! ----------------------------------------------------------------------------- !
      !     III    Net longwave radiative FLUX                                        !
      ! ----------------------------------------------------------------------------- !

      !! LB: now moved after Turbulent fluxes because must use the skin temperature rather that the SST
      !! (ztskk is skin temperature if ln_skin_cs==.TRUE. .OR. ln_skin_wl==.TRUE.)
      zqlw(:,:) = emiss_w * ( pqlw(:,:) - stefan*ztskk(:,:)*ztskk(:,:)*ztskk(:,:)*ztskk(:,:) ) * tmask(:,:,1)   ! Net radiative longwave flux

      !  Latent flux over ocean
      ! -----------------------

      ! use scalar version of L_vap() for AGRIF compatibility
      DO_2D_11_11
         zqla(ji,jj) = - L_vap( ztskk(ji,jj) ) * pevp(ji,jj)    ! Latent Heat flux !!GS: possibility to add a global qla to avoid recomputation after abl update
      END_2D

      IF(sn_cfctl%l_prtctl) THEN
         CALL prt_ctl( tab2d_1=zqla  , clinfo1=' blk_oce_2: zqla   : ' )
         CALL prt_ctl( tab2d_1=zqlw  , clinfo1=' blk_oce_2: zqlw   : ', tab2d_2=qsr, clinfo2=' qsr : ' )

      ENDIF

      ! ----------------------------------------------------------------------------- !
      !     IV    Total FLUXES                                                       !
      ! ----------------------------------------------------------------------------- !
      !
      emp (:,:) = (  pevp(:,:)                                       &   ! mass flux (evap. - precip.)
         &         - pprec(:,:) * rn_pfac  ) * tmask(:,:,1)
      !
      qns(:,:) = zqlw(:,:) + psen(:,:) + zqla(:,:)                   &   ! Downward Non Solar
         &     - psnow(:,:) * rn_pfac * rLfus                        &   ! remove latent melting heat for solid precip
         &     - pevp(:,:) * ptsk(:,:) * rcp                         &   ! remove evap heat content at SST
         &     + ( pprec(:,:) - psnow(:,:) ) * rn_pfac               &   ! add liquid precip heat content at Tair
         &     * ( ptair(:,:) - rt0 ) * rcp                          &
         &     + psnow(:,:) * rn_pfac                                &   ! add solid  precip heat content at min(Tair,Tsnow)
         &     * ( MIN( ptair(:,:), rt0 ) - rt0 ) * rcpi
      qns(:,:) = qns(:,:) * tmask(:,:,1)
      !
#if defined key_si3
      qns_oce(:,:) = zqlw(:,:) + psen(:,:) + zqla(:,:)                             ! non solar without emp (only needed by SI3)
      qsr_oce(:,:) = qsr(:,:)
#endif
      !
      CALL iom_put( "rho_air"  , rhoa*tmask(:,:,1) )       ! output air density [kg/m^3]
      CALL iom_put( "evap_oce" , pevp )                    ! evaporation
      CALL iom_put( "qlw_oce"  , zqlw )                    ! output downward longwave heat over the ocean
      CALL iom_put( "qsb_oce"  , psen )                    ! output downward sensible heat over the ocean
      CALL iom_put( "qla_oce"  , zqla )                    ! output downward latent   heat over the ocean
      tprecip(:,:) = pprec(:,:) * rn_pfac * tmask(:,:,1)   ! output total precipitation [kg/m2/s]
      sprecip(:,:) = psnow(:,:) * rn_pfac * tmask(:,:,1)   ! output solid precipitation [kg/m2/s]
      CALL iom_put( 'snowpre', sprecip )                   ! Snow
      CALL iom_put( 'precip' , tprecip )                   ! Total precipitation
      !
      IF ( nn_ice == 0 ) THEN
         CALL iom_put( "qemp_oce" , qns-zqlw-psen-zqla )   ! output downward heat content of E-P over the ocean
         CALL iom_put( "qns_oce"  ,   qns  )               ! output downward non solar heat over the ocean
         CALL iom_put( "qsr_oce"  ,   qsr  )               ! output downward solar heat over the ocean
         CALL iom_put( "qt_oce"   ,   qns+qsr )            ! output total downward heat over the ocean
      ENDIF
      !
      IF(sn_cfctl%l_prtctl) THEN
         CALL prt_ctl(tab2d_1=zqlw , clinfo1=' blk_oce_2: zqlw  : ')
         CALL prt_ctl(tab2d_1=zqla , clinfo1=' blk_oce_2: zqla  : ', tab2d_2=qsr  , clinfo2=' qsr   : ')
         CALL prt_ctl(tab2d_1=emp  , clinfo1=' blk_oce_2: emp   : ')
      ENDIF
      !
   END SUBROUTINE blk_oce_2


#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   blk_ice_1   : provide the air-ice stress
   !!   blk_ice_2   : provide the heat and mass fluxes at air-ice interface
   !!   blk_ice_qcn : provide ice surface temperature and snow/ice conduction flux (emulating conduction flux)
   !!   Cdn10_Lupkes2012 : Lupkes et al. (2012) air-ice drag
   !!   Cdn10_Lupkes2015 : Lupkes et al. (2015) air-ice drag
   !!----------------------------------------------------------------------

   SUBROUTINE blk_ice_1( pwndi, pwndj, ptair, phumi, pslp , puice, pvice, ptsui,  &   ! inputs
      &                  putaui, pvtaui, pseni, pevpi, pssqi, pcd_dui             )   ! optional outputs
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_1  ***
      !!
      !! ** Purpose :   provide the surface boundary condition over sea-ice
      !!
      !! ** Method  :   compute momentum using bulk formulation
      !!                formulea, ice variables and read atmospheric fields.
      !!                NB: ice drag coefficient is assumed to be a constant
      !!---------------------------------------------------------------------
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pslp    ! sea-level pressure [Pa]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pwndi   ! atmospheric wind at T-point [m/s]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pwndj   ! atmospheric wind at T-point [m/s]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   ptair   ! atmospheric wind at T-point [m/s]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   phumi   ! atmospheric wind at T-point [m/s]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   puice   ! sea-ice velocity on I or C grid [m/s]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pvice   ! "
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   ptsui   ! sea-ice surface temperature [K]
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   putaui  ! if ln_blk
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   pvtaui  ! if ln_blk
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   pseni   ! if ln_abl
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   pevpi   ! if ln_abl
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   pssqi   ! if ln_abl
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ), OPTIONAL ::   pcd_dui ! if ln_abl
      !
      INTEGER  ::   ji, jj    ! dummy loop indices
      REAL(wp) ::   zwndi_t , zwndj_t             ! relative wind components at T-point
      REAL(wp) ::   zootm_su                      ! sea-ice surface mean temperature
      REAL(wp) ::   zztmp1, zztmp2                ! temporary arrays
      REAL(wp), DIMENSION(jpi,jpj) ::   zcd_dui   ! transfer coefficient for momentum      (tau)
      !!---------------------------------------------------------------------
      !

      ! ------------------------------------------------------------ !
      !    Wind module relative to the moving ice ( U10m - U_ice )   !
      ! ------------------------------------------------------------ !
      ! C-grid ice dynamics :   U & V-points (same as ocean)
      DO_2D_00_00
         zwndi_t = (  pwndi(ji,jj) - rn_vfac * 0.5_wp * ( puice(ji-1,jj  ) + puice(ji,jj) )  )
         zwndj_t = (  pwndj(ji,jj) - rn_vfac * 0.5_wp * ( pvice(ji  ,jj-1) + pvice(ji,jj) )  )
         wndm_ice(ji,jj) = SQRT( zwndi_t * zwndi_t + zwndj_t * zwndj_t ) * tmask(ji,jj,1)
      END_2D
      CALL lbc_lnk( 'sbcblk', wndm_ice, 'T',  1. )
      !
      ! Make ice-atm. drag dependent on ice concentration
      IF    ( ln_Cd_L12 ) THEN   ! calculate new drag from Lupkes(2012) equations
         CALL Cdn10_Lupkes2012( Cd_ice )
         Ch_ice(:,:) = Cd_ice(:,:)       ! momentum and heat transfer coef. are considered identical
         Ce_ice(:,:) = Cd_ice(:,:)
      ELSEIF( ln_Cd_L15 ) THEN   ! calculate new drag from Lupkes(2015) equations
         CALL Cdn10_Lupkes2015( ptsui, pslp, Cd_ice, Ch_ice )
         Ce_ice(:,:) = Ch_ice(:,:)       ! sensible and latent heat transfer coef. are considered identical
      ENDIF

      !! IF ( iom_use("Cd_ice") ) CALL iom_put("Cd_ice", Cd_ice)   ! output value of pure ice-atm. transfer coef.
      !! IF ( iom_use("Ch_ice") ) CALL iom_put("Ch_ice", Ch_ice)   ! output value of pure ice-atm. transfer coef.

      ! local scalars ( place there for vector optimisation purposes)
      !IF (ln_abl) rhoa  (:,:)  = rho_air( ptair(:,:), phumi(:,:), pslp(:,:) ) !!GS: rhoa must be (re)computed here with ABL to avoid division by zero after (TBI)
      zcd_dui(:,:) = wndm_ice(:,:) * Cd_ice(:,:)

      IF( ln_blk ) THEN
         ! ------------------------------------------------------------ !
         !    Wind stress relative to the moving ice ( U10m - U_ice )   !
         ! ------------------------------------------------------------ !
         ! C-grid ice dynamics :   U & V-points (same as ocean)
         DO_2D_00_00
            putaui(ji,jj) = 0.5_wp * (  rhoa(ji+1,jj) * zcd_dui(ji+1,jj)             &
               &                      + rhoa(ji  ,jj) * zcd_dui(ji  ,jj)  )          &
               &         * ( 0.5_wp * ( pwndi(ji+1,jj) + pwndi(ji,jj) ) - rn_vfac * puice(ji,jj) )
            pvtaui(ji,jj) = 0.5_wp * (  rhoa(ji,jj+1) * zcd_dui(ji,jj+1)             &
               &                      + rhoa(ji,jj  ) * zcd_dui(ji,jj  )  )          &
               &         * ( 0.5_wp * ( pwndj(ji,jj+1) + pwndj(ji,jj) ) - rn_vfac * pvice(ji,jj) )
         END_2D
         CALL lbc_lnk_multi( 'sbcblk', putaui, 'U', -1., pvtaui, 'V', -1. )
         !
         IF(sn_cfctl%l_prtctl)  CALL prt_ctl( tab2d_1=putaui  , clinfo1=' blk_ice: putaui : '   &
            &                               , tab2d_2=pvtaui  , clinfo2='          pvtaui : ' )
      ELSE
         zztmp1 = 11637800.0_wp
         zztmp2 =    -5897.8_wp
         DO_2D_11_11
            pcd_dui(ji,jj) = zcd_dui (ji,jj)
            pseni  (ji,jj) = wndm_ice(ji,jj) * Ch_ice(ji,jj)
            pevpi  (ji,jj) = wndm_ice(ji,jj) * Ce_ice(ji,jj)
            zootm_su       = zztmp2 / ptsui(ji,jj)   ! ptsui is in K (it can't be zero ??)
            pssqi  (ji,jj) = zztmp1 * EXP( zootm_su ) / rhoa(ji,jj)
         END_2D
      ENDIF
      !
      IF(sn_cfctl%l_prtctl)  CALL prt_ctl(tab2d_1=wndm_ice  , clinfo1=' blk_ice: wndm_ice : ')
      !
   END SUBROUTINE blk_ice_1


   SUBROUTINE blk_ice_2( ptsu, phs, phi, palb, ptair, phumi, pslp, pqlw, pprec, psnow  )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_2  ***
      !!
      !! ** Purpose :   provide the heat and mass fluxes at air-ice interface
      !!
      !! ** Method  :   compute heat and freshwater exchanged
      !!                between atmosphere and sea-ice using bulk formulation
      !!                formulea, ice variables and read atmmospheric fields.
      !!
      !! caution : the net upward water flux has with mm/day unit
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   ptsu   ! sea ice surface temperature [K]
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   palb   ! ice albedo (all skies)
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   ptair
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   phumi
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   pslp
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   pqlw
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   pprec
      REAL(wp), DIMENSION(:,:  ), INTENT(in)  ::   psnow
      !!
      INTEGER  ::   ji, jj, jl               ! dummy loop indices
      REAL(wp) ::   zst3                     ! local variable
      REAL(wp) ::   zcoef_dqlw, zcoef_dqla   !   -      -
      REAL(wp) ::   zztmp, zztmp2, z1_rLsub  !   -      -
      REAL(wp) ::   zfr1, zfr2               ! local variables
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z1_st         ! inverse of surface temperature
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_qlw         ! long wave heat flux over ice
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_qsb         ! sensible  heat flux over ice
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_dqlw        ! long wave heat sensitivity over ice
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   z_dqsb        ! sensible  heat sensitivity over ice
      REAL(wp), DIMENSION(jpi,jpj)     ::   zevap, zsnw   ! evaporation and snw distribution after wind blowing (SI3)
      REAL(wp), DIMENSION(jpi,jpj)     ::   zqair         ! specific humidity of air at z=rn_zqt [kg/kg] !LB
      REAL(wp), DIMENSION(jpi,jpj)     ::   ztmp, ztmp2
      !!---------------------------------------------------------------------
      !
      zcoef_dqlw = 4._wp * 0.95_wp * stefan             ! local scalars
      zcoef_dqla = -rLsub * 11637800._wp * (-5897.8_wp) !LB: BAD!
      !
      SELECT CASE( nhumi )
      CASE( np_humi_sph )
         zqair(:,:) =  phumi(:,:)      ! what we read in file is already a spec. humidity!
      CASE( np_humi_dpt )
         zqair(:,:) = q_sat( phumi(:,:), pslp )
      CASE( np_humi_rlh )
         zqair(:,:) = q_air_rh( 0.01_wp*phumi(:,:), ptair(:,:), pslp(:,:) ) !LB: 0.01 => RH is % percent in file
      END SELECT
      !
      zztmp = 1. / ( 1. - albo )
      WHERE( ptsu(:,:,:) /= 0._wp )
         z1_st(:,:,:) = 1._wp / ptsu(:,:,:)
      ELSEWHERE
         z1_st(:,:,:) = 0._wp
      END WHERE
      !                                     ! ========================== !
      DO jl = 1, jpl                        !  Loop over ice categories  !
         !                                  ! ========================== !
         DO jj = 1 , jpj
            DO ji = 1, jpi
               ! ----------------------------!
               !      I   Radiative FLUXES   !
               ! ----------------------------!
               zst3 = ptsu(ji,jj,jl) * ptsu(ji,jj,jl) * ptsu(ji,jj,jl)
               ! Short Wave (sw)
               qsr_ice(ji,jj,jl) = zztmp * ( 1. - palb(ji,jj,jl) ) * qsr(ji,jj)
               ! Long  Wave (lw)
               z_qlw(ji,jj,jl) = 0.95 * ( pqlw(ji,jj) - stefan * ptsu(ji,jj,jl) * zst3 ) * tmask(ji,jj,1)
               ! lw sensitivity
               z_dqlw(ji,jj,jl) = zcoef_dqlw * zst3

               ! ----------------------------!
               !     II    Turbulent FLUXES  !
               ! ----------------------------!

               ! ... turbulent heat fluxes with Ch_ice recalculated in blk_ice_1
               ! Sensible Heat
               z_qsb(ji,jj,jl) = rhoa(ji,jj) * rCp_air * Ch_ice(ji,jj) * wndm_ice(ji,jj) * (ptsu(ji,jj,jl) - ptair(ji,jj))
               ! Latent Heat
               zztmp2 = EXP( -5897.8 * z1_st(ji,jj,jl) )
               qla_ice(ji,jj,jl) = rn_efac * MAX( 0.e0, rhoa(ji,jj) * rLsub  * Ce_ice(ji,jj) * wndm_ice(ji,jj) *  &
                  &                ( 11637800. * zztmp2 / rhoa(ji,jj) - zqair(ji,jj) ) )
               ! Latent heat sensitivity for ice (Dqla/Dt)
               IF( qla_ice(ji,jj,jl) > 0._wp ) THEN
                  dqla_ice(ji,jj,jl) = rn_efac * zcoef_dqla * Ce_ice(ji,jj) * wndm_ice(ji,jj) *  &
                     &                 z1_st(ji,jj,jl) * z1_st(ji,jj,jl) * zztmp2
               ELSE
                  dqla_ice(ji,jj,jl) = 0._wp
               ENDIF

               ! Sensible heat sensitivity (Dqsb_ice/Dtn_ice)
               z_dqsb(ji,jj,jl) = rhoa(ji,jj) * rCp_air * Ch_ice(ji,jj) * wndm_ice(ji,jj)

               ! ----------------------------!
               !     III    Total FLUXES     !
               ! ----------------------------!
               ! Downward Non Solar flux
               qns_ice (ji,jj,jl) =     z_qlw (ji,jj,jl) - z_qsb (ji,jj,jl) - qla_ice (ji,jj,jl)
               ! Total non solar heat flux sensitivity for ice
               dqns_ice(ji,jj,jl) = - ( z_dqlw(ji,jj,jl) + z_dqsb(ji,jj,jl) + dqla_ice(ji,jj,jl) )
            END DO
            !
         END DO
         !
      END DO
      !
      tprecip(:,:) = pprec(:,:) * rn_pfac * tmask(:,:,1)  ! total precipitation [kg/m2/s]
      sprecip(:,:) = psnow(:,:) * rn_pfac * tmask(:,:,1)  ! solid precipitation [kg/m2/s]
      CALL iom_put( 'snowpre', sprecip )                  ! Snow precipitation
      CALL iom_put( 'precip' , tprecip )                  ! Total precipitation

      ! --- evaporation --- !
      z1_rLsub = 1._wp / rLsub
      evap_ice (:,:,:) = rn_efac * qla_ice (:,:,:) * z1_rLsub    ! sublimation
      devap_ice(:,:,:) = rn_efac * dqla_ice(:,:,:) * z1_rLsub    ! d(sublimation)/dT
      zevap    (:,:)   = emp(:,:) + tprecip(:,:)   ! evaporation over ocean  !LB: removed rn_efac here, correct???

      ! --- evaporation minus precipitation --- !
      zsnw(:,:) = 0._wp
      CALL ice_thd_snwblow( (1.-at_i_b(:,:)), zsnw )  ! snow distribution over ice after wind blowing
      emp_oce(:,:) = ( 1._wp - at_i_b(:,:) ) * zevap(:,:) - ( tprecip(:,:) - sprecip(:,:) ) - sprecip(:,:) * (1._wp - zsnw )
      emp_ice(:,:) = SUM( a_i_b(:,:,:) * evap_ice(:,:,:), dim=3 ) - sprecip(:,:) * zsnw
      emp_tot(:,:) = emp_oce(:,:) + emp_ice(:,:)

      ! --- heat flux associated with emp --- !
      qemp_oce(:,:) = - ( 1._wp - at_i_b(:,:) ) * zevap(:,:) * sst_m(:,:) * rcp                  & ! evap at sst
         &          + ( tprecip(:,:) - sprecip(:,:) ) * ( ptair(:,:) - rt0 ) * rcp               & ! liquid precip at Tair
         &          +   sprecip(:,:) * ( 1._wp - zsnw ) *                                        & ! solid precip at min(Tair,Tsnow)
         &              ( ( MIN( ptair(:,:), rt0 ) - rt0 ) * rcpi * tmask(:,:,1) - rLfus )
      qemp_ice(:,:) =   sprecip(:,:) * zsnw *                                                    & ! solid precip (only)
         &              ( ( MIN( ptair(:,:), rt0 ) - rt0 ) * rcpi * tmask(:,:,1) - rLfus )

      ! --- total solar and non solar fluxes --- !
      qns_tot(:,:) = ( 1._wp - at_i_b(:,:) ) * qns_oce(:,:) + SUM( a_i_b(:,:,:) * qns_ice(:,:,:), dim=3 )  &
         &           + qemp_ice(:,:) + qemp_oce(:,:)
      qsr_tot(:,:) = ( 1._wp - at_i_b(:,:) ) * qsr_oce(:,:) + SUM( a_i_b(:,:,:) * qsr_ice(:,:,:), dim=3 )

      ! --- heat content of precip over ice in J/m3 (to be used in 1D-thermo) --- !
      qprec_ice(:,:) = rhos * ( ( MIN( ptair(:,:), rt0 ) - rt0 ) * rcpi * tmask(:,:,1) - rLfus )

      ! --- heat content of evap over ice in W/m2 (to be used in 1D-thermo) ---
      DO jl = 1, jpl
         qevap_ice(:,:,jl) = 0._wp ! should be -evap_ice(:,:,jl)*( ( Tice - rt0 ) * rcpi * tmask(:,:,1) )
         !                         ! But we do not have Tice => consider it at 0degC => evap=0
      END DO

      ! --- shortwave radiation transmitted below the surface (W/m2, see Grenfell Maykut 77) --- !
      zfr1 = ( 0.18 * ( 1.0 - cldf_ice ) + 0.35 * cldf_ice )            ! transmission when hi>10cm
      zfr2 = ( 0.82 * ( 1.0 - cldf_ice ) + 0.65 * cldf_ice )            ! zfr2 such that zfr1 + zfr2 to equal 1
      !
      WHERE    ( phs(:,:,:) <= 0._wp .AND. phi(:,:,:) <  0.1_wp )       ! linear decrease from hi=0 to 10cm
         qtr_ice_top(:,:,:) = qsr_ice(:,:,:) * ( zfr1 + zfr2 * ( 1._wp - phi(:,:,:) * 10._wp ) )
      ELSEWHERE( phs(:,:,:) <= 0._wp .AND. phi(:,:,:) >= 0.1_wp )       ! constant (zfr1) when hi>10cm
         qtr_ice_top(:,:,:) = qsr_ice(:,:,:) * zfr1
      ELSEWHERE                                                         ! zero when hs>0
         qtr_ice_top(:,:,:) = 0._wp
      END WHERE
      !

      IF( iom_use('evap_ao_cea') .OR. iom_use('hflx_evap_cea') ) THEN
         ztmp(:,:) = zevap(:,:) * ( 1._wp - at_i_b(:,:) )
         IF( iom_use('evap_ao_cea'  ) )  CALL iom_put( 'evap_ao_cea'  , ztmp(:,:) * tmask(:,:,1) )   ! ice-free oce evap (cell average)
         IF( iom_use('hflx_evap_cea') )  CALL iom_put( 'hflx_evap_cea', ztmp(:,:) * sst_m(:,:) * rcp * tmask(:,:,1) )   ! heat flux from evap (cell average)
      ENDIF
      IF( iom_use('hflx_rain_cea') ) THEN
         ztmp(:,:) = rcp * ( SUM( (ptsu-rt0) * a_i_b, dim=3 ) + sst_m(:,:) * ( 1._wp - at_i_b(:,:) ) )
         IF( iom_use('hflx_rain_cea') )  CALL iom_put( 'hflx_rain_cea', ( tprecip(:,:) - sprecip(:,:) ) * ztmp(:,:) )   ! heat flux from rain (cell average)
      ENDIF
      IF( iom_use('hflx_snow_cea') .OR. iom_use('hflx_snow_ao_cea') .OR. iom_use('hflx_snow_ai_cea')  )  THEN
         WHERE( SUM( a_i_b, dim=3 ) > 1.e-10 )
            ztmp(:,:) = rcpi * SUM( (ptsu-rt0) * a_i_b, dim=3 ) / SUM( a_i_b, dim=3 )
         ELSEWHERE
            ztmp(:,:) = rcp * sst_m(:,:)
         ENDWHERE
         ztmp2(:,:) = sprecip(:,:) * ( ztmp(:,:) - rLfus )
         IF( iom_use('hflx_snow_cea')    ) CALL iom_put('hflx_snow_cea'   , ztmp2(:,:) ) ! heat flux from snow (cell average)
         IF( iom_use('hflx_snow_ao_cea') ) CALL iom_put('hflx_snow_ao_cea', ztmp2(:,:) * ( 1._wp - zsnw(:,:) ) ) ! heat flux from snow (over ocean)
         IF( iom_use('hflx_snow_ai_cea') ) CALL iom_put('hflx_snow_ai_cea', ztmp2(:,:) *           zsnw(:,:)   ) ! heat flux from snow (over ice)
      ENDIF
      !
      IF(sn_cfctl%l_prtctl) THEN
         CALL prt_ctl(tab3d_1=qla_ice , clinfo1=' blk_ice: qla_ice  : ', tab3d_2=z_qsb   , clinfo2=' z_qsb    : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=z_qlw   , clinfo1=' blk_ice: z_qlw    : ', tab3d_2=dqla_ice, clinfo2=' dqla_ice : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=z_dqsb  , clinfo1=' blk_ice: z_dqsb   : ', tab3d_2=z_dqlw  , clinfo2=' z_dqlw   : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=dqns_ice, clinfo1=' blk_ice: dqns_ice : ', tab3d_2=qsr_ice , clinfo2=' qsr_ice  : ', kdim=jpl)
         CALL prt_ctl(tab3d_1=ptsu    , clinfo1=' blk_ice: ptsu     : ', tab3d_2=qns_ice , clinfo2=' qns_ice  : ', kdim=jpl)
         CALL prt_ctl(tab2d_1=tprecip , clinfo1=' blk_ice: tprecip  : ', tab2d_2=sprecip , clinfo2=' sprecip  : ')
      ENDIF
      !
   END SUBROUTINE blk_ice_2


   SUBROUTINE blk_ice_qcn( ld_virtual_itd, ptsu, ptb, phs, phi )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE blk_ice_qcn  ***
      !!
      !! ** Purpose :   Compute surface temperature and snow/ice conduction flux
      !!                to force sea ice / snow thermodynamics
      !!                in the case conduction flux is emulated
      !!
      !! ** Method  :   compute surface energy balance assuming neglecting heat storage
      !!                following the 0-layer Semtner (1976) approach
      !!
      !! ** Outputs : - ptsu    : sea-ice / snow surface temperature (K)
      !!              - qcn_ice : surface inner conduction flux (W/m2)
      !!
      !!---------------------------------------------------------------------
      LOGICAL                   , INTENT(in   ) ::   ld_virtual_itd  ! single-category option
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   ptsu            ! sea ice / snow surface temperature
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   ptb             ! sea ice base temperature
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   phs             ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   phi             ! sea ice thickness
      !
      INTEGER , PARAMETER ::   nit = 10                  ! number of iterations
      REAL(wp), PARAMETER ::   zepsilon = 0.1_wp         ! characteristic thickness for enhanced conduction
      !
      INTEGER  ::   ji, jj, jl           ! dummy loop indices
      INTEGER  ::   iter                 ! local integer
      REAL(wp) ::   zfac, zfac2, zfac3   ! local scalars
      REAL(wp) ::   zkeff_h, ztsu, ztsu0 !
      REAL(wp) ::   zqc, zqnet           !
      REAL(wp) ::   zhe, zqa0            !
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zgfac   ! enhanced conduction factor
      !!---------------------------------------------------------------------

      ! -------------------------------------!
      !      I   Enhanced conduction factor  !
      ! -------------------------------------!
      ! Emulates the enhancement of conduction by unresolved thin ice (ld_virtual_itd = T)
      ! Fichefet and Morales Maqueda, JGR 1997
      !
      zgfac(:,:,:) = 1._wp

      IF( ld_virtual_itd ) THEN
         !
         zfac  = 1._wp /  ( rn_cnd_s + rcnd_i )
         zfac2 = EXP(1._wp) * 0.5_wp * zepsilon
         zfac3 = 2._wp / zepsilon
         !
         DO jl = 1, jpl
            DO_2D_11_11
               zhe = ( rn_cnd_s * phi(ji,jj,jl) + rcnd_i * phs(ji,jj,jl) ) * zfac                            ! Effective thickness
               IF( zhe >=  zfac2 )   zgfac(ji,jj,jl) = MIN( 2._wp, 0.5_wp * ( 1._wp + LOG( zhe * zfac3 ) ) ) ! Enhanced conduction factor
            END_2D
         END DO
         !
      ENDIF

      ! -------------------------------------------------------------!
      !      II   Surface temperature and conduction flux            !
      ! -------------------------------------------------------------!
      !
      zfac = rcnd_i * rn_cnd_s
      !
      DO jl = 1, jpl
         DO_2D_11_11
            !
            zkeff_h = zfac * zgfac(ji,jj,jl) / &                                    ! Effective conductivity of the snow-ice system divided by thickness
               &      ( rcnd_i * phs(ji,jj,jl) + rn_cnd_s * MAX( 0.01, phi(ji,jj,jl) ) )
            ztsu    = ptsu(ji,jj,jl)                                                ! Store current iteration temperature
            ztsu0   = ptsu(ji,jj,jl)                                                ! Store initial surface temperature
            zqa0    = qsr_ice(ji,jj,jl) - qtr_ice_top(ji,jj,jl) + qns_ice(ji,jj,jl) ! Net initial atmospheric heat flux
            !
            DO iter = 1, nit     ! --- Iterative loop
               zqc   = zkeff_h * ( ztsu - ptb(ji,jj) )                              ! Conduction heat flux through snow-ice system (>0 downwards)
               zqnet = zqa0 + dqns_ice(ji,jj,jl) * ( ztsu - ptsu(ji,jj,jl) ) - zqc  ! Surface energy budget
               ztsu  = ztsu - zqnet / ( dqns_ice(ji,jj,jl) - zkeff_h )              ! Temperature update
            END DO
            !
            ptsu   (ji,jj,jl) = MIN( rt0, ztsu )
            qcn_ice(ji,jj,jl) = zkeff_h * ( ptsu(ji,jj,jl) - ptb(ji,jj) )
            qns_ice(ji,jj,jl) = qns_ice(ji,jj,jl) + dqns_ice(ji,jj,jl) * ( ptsu(ji,jj,jl) - ztsu0 )
            qml_ice(ji,jj,jl) = ( qsr_ice(ji,jj,jl) - qtr_ice_top(ji,jj,jl) + qns_ice(ji,jj,jl) - qcn_ice(ji,jj,jl) )  &
               &   * MAX( 0._wp , SIGN( 1._wp, ptsu(ji,jj,jl) - rt0 ) )

            ! --- Diagnose the heat loss due to changing non-solar flux (as in icethd_zdf_bl99) --- !
            hfx_err_dif(ji,jj) = hfx_err_dif(ji,jj) - ( dqns_ice(ji,jj,jl) * ( ptsu(ji,jj,jl) - ztsu0 ) ) * a_i_b(ji,jj,jl)

         END_2D
         !
      END DO
      !
   END SUBROUTINE blk_ice_qcn


   SUBROUTINE Cdn10_Lupkes2012( pcd )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  Cdn10_Lupkes2012  ***
      !!
      !! ** Purpose :    Recompute the neutral air-ice drag referenced at 10m
      !!                 to make it dependent on edges at leads, melt ponds and flows.
      !!                 After some approximations, this can be resumed to a dependency
      !!                 on ice concentration.
      !!
      !! ** Method :     The parameterization is taken from Lupkes et al. (2012) eq.(50)
      !!                 with the highest level of approximation: level4, eq.(59)
      !!                 The generic drag over a cell partly covered by ice can be re-written as follows:
      !!
      !!                 Cd = Cdw * (1-A) + Cdi * A + Ce * (1-A)**(nu+1/(10*beta)) * A**mu
      !!
      !!                    Ce = 2.23e-3       , as suggested by Lupkes (eq. 59)
      !!                    nu = mu = beta = 1 , as suggested by Lupkes (eq. 59)
      !!                    A is the concentration of ice minus melt ponds (if any)
      !!
      !!                 This new drag has a parabolic shape (as a function of A) starting at
      !!                 Cdw(say 1.5e-3) for A=0, reaching 1.97e-3 for A~0.5
      !!                 and going down to Cdi(say 1.4e-3) for A=1
      !!
      !!                 It is theoretically applicable to all ice conditions (not only MIZ)
      !!                 => see Lupkes et al (2013)
      !!
      !! ** References : Lupkes et al. JGR 2012 (theory)
      !!                 Lupkes et al. GRL 2013 (application to GCM)
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pcd
      REAL(wp), PARAMETER ::   zCe   = 2.23e-03_wp
      REAL(wp), PARAMETER ::   znu   = 1._wp
      REAL(wp), PARAMETER ::   zmu   = 1._wp
      REAL(wp), PARAMETER ::   zbeta = 1._wp
      REAL(wp)            ::   zcoef
      !!----------------------------------------------------------------------
      zcoef = znu + 1._wp / ( 10._wp * zbeta )

      ! generic drag over a cell partly covered by ice
      !!Cd(:,:) = Cd_oce(:,:) * ( 1._wp - at_i_b(:,:) ) +  &                        ! pure ocean drag
      !!   &      Cd_ice      *           at_i_b(:,:)   +  &                        ! pure ice drag
      !!   &      zCe         * ( 1._wp - at_i_b(:,:) )**zcoef * at_i_b(:,:)**zmu   ! change due to sea-ice morphology

      ! ice-atm drag
      pcd(:,:) = rCd_ice +  &                                                         ! pure ice drag
         &      zCe     * ( 1._wp - at_i_b(:,:) )**zcoef * at_i_b(:,:)**(zmu-1._wp)  ! change due to sea-ice morphology

   END SUBROUTINE Cdn10_Lupkes2012


   SUBROUTINE Cdn10_Lupkes2015( ptm_su, pslp, pcd, pch )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE  Cdn10_Lupkes2015  ***
      !!
      !! ** pUrpose :    Alternative turbulent transfert coefficients formulation
      !!                 between sea-ice and atmosphere with distinct momentum
      !!                 and heat coefficients depending on sea-ice concentration
      !!                 and atmospheric stability (no meltponds effect for now).
      !!
      !! ** Method :     The parameterization is adapted from Lupkes et al. (2015)
      !!                 and ECHAM6 atmospheric model. Compared to Lupkes2012 scheme,
      !!                 it considers specific skin and form drags (Andreas et al. 2010)
      !!                 to compute neutral transfert coefficients for both heat and
      !!                 momemtum fluxes. Atmospheric stability effect on transfert
      !!                 coefficient is also taken into account following Louis (1979).
      !!
      !! ** References : Lupkes et al. JGR 2015 (theory)
      !!                 Lupkes et al. ECHAM6 documentation 2015 (implementation)
      !!
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   ptm_su ! sea-ice surface temperature [K]
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pslp   ! sea-level pressure [Pa]
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pcd    ! momentum transfert coefficient
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pch    ! heat transfert coefficient
      REAL(wp), DIMENSION(jpi,jpj)            ::   zst, zqo_sat, zqi_sat
      !
      ! ECHAM6 constants
      REAL(wp), PARAMETER ::   z0_skin_ice  = 0.69e-3_wp  ! Eq. 43 [m]
      REAL(wp), PARAMETER ::   z0_form_ice  = 0.57e-3_wp  ! Eq. 42 [m]
      REAL(wp), PARAMETER ::   z0_ice       = 1.00e-3_wp  ! Eq. 15 [m]
      REAL(wp), PARAMETER ::   zce10        = 2.80e-3_wp  ! Eq. 41
      REAL(wp), PARAMETER ::   zbeta        = 1.1_wp      ! Eq. 41
      REAL(wp), PARAMETER ::   zc           = 5._wp       ! Eq. 13
      REAL(wp), PARAMETER ::   zc2          = zc * zc
      REAL(wp), PARAMETER ::   zam          = 2. * zc     ! Eq. 14
      REAL(wp), PARAMETER ::   zah          = 3. * zc     ! Eq. 30
      REAL(wp), PARAMETER ::   z1_alpha     = 1._wp / 0.2_wp  ! Eq. 51
      REAL(wp), PARAMETER ::   z1_alphaf    = z1_alpha    ! Eq. 56
      REAL(wp), PARAMETER ::   zbetah       = 1.e-3_wp    ! Eq. 26
      REAL(wp), PARAMETER ::   zgamma       = 1.25_wp     ! Eq. 26
      REAL(wp), PARAMETER ::   z1_gamma     = 1._wp / zgamma
      REAL(wp), PARAMETER ::   r1_3         = 1._wp / 3._wp
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::   zthetav_os, zthetav_is, zthetav_zu
      REAL(wp) ::   zrib_o, zrib_i
      REAL(wp) ::   zCdn_skin_ice, zCdn_form_ice, zCdn_ice
      REAL(wp) ::   zChn_skin_ice, zChn_form_ice
      REAL(wp) ::   z0w, z0i, zfmi, zfmw, zfhi, zfhw
      REAL(wp) ::   zCdn_form_tmp
      !!----------------------------------------------------------------------

      ! Momentum Neutral Transfert Coefficients (should be a constant)
      zCdn_form_tmp = zce10 * ( LOG( 10._wp / z0_form_ice + 1._wp ) / LOG( rn_zu / z0_form_ice + 1._wp ) )**2   ! Eq. 40
      zCdn_skin_ice = ( vkarmn                                      / LOG( rn_zu / z0_skin_ice + 1._wp ) )**2   ! Eq. 7
      zCdn_ice      = zCdn_skin_ice   ! Eq. 7
      !zCdn_ice     = 1.89e-3         ! old ECHAM5 value (cf Eq. 32)

      ! Heat Neutral Transfert Coefficients
      zChn_skin_ice = vkarmn**2 / ( LOG( rn_zu / z0_ice + 1._wp ) * LOG( rn_zu * z1_alpha / z0_skin_ice + 1._wp ) )   ! Eq. 50 + Eq. 52

      ! Atmospheric and Surface Variables
      zst(:,:)     = sst_m(:,:) + rt0                                        ! convert SST from Celcius to Kelvin
      zqo_sat(:,:) = rdct_qsat_salt * q_sat( zst(:,:)   , pslp(:,:) )   ! saturation humidity over ocean [kg/kg]
      zqi_sat(:,:) =                  q_sat( ptm_su(:,:), pslp(:,:) )   ! saturation humidity over ice   [kg/kg]
      !
      DO_2D_00_00
         ! Virtual potential temperature [K]
         zthetav_os = zst(ji,jj)    * ( 1._wp + rctv0 * zqo_sat(ji,jj) )   ! over ocean
         zthetav_is = ptm_su(ji,jj) * ( 1._wp + rctv0 * zqi_sat(ji,jj) )   ! ocean ice
         zthetav_zu = t_zu (ji,jj)  * ( 1._wp + rctv0 * q_zu(ji,jj)    )   ! at zu

         ! Bulk Richardson Number (could use Ri_bulk function from aerobulk instead)
         zrib_o = grav / zthetav_os * ( zthetav_zu - zthetav_os ) * rn_zu / MAX( 0.5, wndm(ji,jj)     )**2   ! over ocean
         zrib_i = grav / zthetav_is * ( zthetav_zu - zthetav_is ) * rn_zu / MAX( 0.5, wndm_ice(ji,jj) )**2   ! over ice

         ! Momentum and Heat Neutral Transfert Coefficients
         zCdn_form_ice = zCdn_form_tmp * at_i_b(ji,jj) * ( 1._wp - at_i_b(ji,jj) )**zbeta  ! Eq. 40
         zChn_form_ice = zCdn_form_ice / ( 1._wp + ( LOG( z1_alphaf ) / vkarmn ) * SQRT( zCdn_form_ice ) )               ! Eq. 53

         ! Momentum and Heat Stability functions (possibility to use psi_m_ecmwf instead ?)
         z0w = rn_zu * EXP( -1._wp * vkarmn / SQRT( Cdn_oce(ji,jj) ) ) ! over water
         z0i = z0_skin_ice                                             ! over ice
         IF( zrib_o <= 0._wp ) THEN
            zfmw = 1._wp - zam * zrib_o / ( 1._wp + 3._wp * zc2 * Cdn_oce(ji,jj) * SQRT( -zrib_o * ( rn_zu / z0w + 1._wp ) ) )  ! Eq. 10
            zfhw = ( 1._wp + ( zbetah * ( zthetav_os - zthetav_zu )**r1_3 / ( Chn_oce(ji,jj) * MAX(0.01, wndm(ji,jj)) )   &     ! Eq. 26
               &             )**zgamma )**z1_gamma
         ELSE
            zfmw = 1._wp / ( 1._wp + zam * zrib_o / SQRT( 1._wp + zrib_o ) )   ! Eq. 12
            zfhw = 1._wp / ( 1._wp + zah * zrib_o / SQRT( 1._wp + zrib_o ) )   ! Eq. 28
         ENDIF

         IF( zrib_i <= 0._wp ) THEN
            zfmi = 1._wp - zam * zrib_i / (1._wp + 3._wp * zc2 * zCdn_ice * SQRT( -zrib_i * ( rn_zu / z0i + 1._wp)))   ! Eq.  9
            zfhi = 1._wp - zah * zrib_i / (1._wp + 3._wp * zc2 * zCdn_ice * SQRT( -zrib_i * ( rn_zu / z0i + 1._wp)))   ! Eq. 25
         ELSE
            zfmi = 1._wp / ( 1._wp + zam * zrib_i / SQRT( 1._wp + zrib_i ) )   ! Eq. 11
            zfhi = 1._wp / ( 1._wp + zah * zrib_i / SQRT( 1._wp + zrib_i ) )   ! Eq. 27
         ENDIF

         ! Momentum Transfert Coefficients (Eq. 38)
         pcd(ji,jj) = zCdn_skin_ice *   zfmi +  &
            &        zCdn_form_ice * ( zfmi * at_i_b(ji,jj) + zfmw * ( 1._wp - at_i_b(ji,jj) ) ) / MAX( 1.e-06, at_i_b(ji,jj) )

         ! Heat Transfert Coefficients (Eq. 49)
         pch(ji,jj) = zChn_skin_ice *   zfhi +  &
            &        zChn_form_ice * ( zfhi * at_i_b(ji,jj) + zfhw * ( 1._wp - at_i_b(ji,jj) ) ) / MAX( 1.e-06, at_i_b(ji,jj) )
         !
      END_2D
      CALL lbc_lnk_multi( 'sbcblk', pcd, 'T',  1., pch, 'T', 1. )
      !
   END SUBROUTINE Cdn10_Lupkes2015

#endif

   !!======================================================================
END MODULE sbcblk
