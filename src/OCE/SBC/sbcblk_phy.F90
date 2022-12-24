MODULE sbcblk_phy
   !!======================================================================
   !!                       ***  MODULE  sbcblk_phy  ***
   !! A set of functions to compute air themodynamics parameters
   !!                     needed by Aerodynamic Bulk Formulas
   !!=====================================================================
   !!            4.0  !  2019 L. Brodeau from AeroBulk package (https://github.com/brodeau/aerobulk/)
   !!----------------------------------------------------------------------

   !!   virt_temp     : virtual (aka sensible) temperature (potential or absolute)
   !!   rho_air       : density of (moist) air (depends on T_air, q_air and SLP
   !!   visc_air      : kinematic viscosity (aka Nu_air) of air from temperature
   !!   L_vap         : latent heat of vaporization of water as a function of temperature
   !!   cp_air        : specific heat of (moist) air (depends spec. hum. q_air)
   !!   gamma_moist   : adiabatic lapse-rate of moist air
   !!   One_on_L      : 1. / ( Monin-Obukhov length )
   !!   Ri_bulk       : bulk Richardson number aka BRN
   !!   q_sat         : saturation humidity as a function of SLP and temperature
   !!   q_air_rh      : specific humidity as a function of RH (fraction, not %), t_air and SLP

   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants

   IMPLICIT NONE
   PRIVATE

   !!  (mainly removed from sbcblk.F90)
   REAL(wp), PARAMETER, PUBLIC :: rCp_dry = 1005.0_wp   !: Specic heat of dry air, constant pressure      [J/K/kg]
   REAL(wp), PARAMETER, PUBLIC :: rCp_vap = 1860.0_wp   !: Specic heat of water vapor, constant pressure  [J/K/kg]
   REAL(wp), PARAMETER, PUBLIC :: R_dry   = 287.05_wp   !: Specific gas constant for dry air              [J/K/kg]
   REAL(wp), PARAMETER, PUBLIC :: R_vap   = 461.495_wp  !: Specific gas constant for water vapor          [J/K/kg]
   REAL(wp), PARAMETER, PUBLIC :: reps0   = R_dry/R_vap !: ratio of gas constant for dry air and water vapor => ~ 0.622
   REAL(wp), PARAMETER, PUBLIC :: rctv0   = R_vap/R_dry !: for virtual temperature (== (1-eps)/eps) => ~ 0.608
   REAL(wp), PARAMETER, PUBLIC :: rCp_air = 1000.5_wp   !: specific heat of air (only used for ice fluxes now...)
   REAL(wp), PARAMETER, PUBLIC :: rCd_ice = 1.4e-3_wp   !: transfer coefficient over ice
   REAL(wp), PARAMETER, PUBLIC :: albo    = 0.066_wp    !: ocean albedo assumed to be constant
   !
   REAL(wp), PARAMETER, PUBLIC :: rho0_a  = 1.2_wp      !: Approx. of density of air                       [kg/m^3]
   REAL(wp), PARAMETER, PUBLIC :: rho0_w  = 1025._wp    !: Density of sea-water  (ECMWF->1025)             [kg/m^3]
   REAL(wp), PARAMETER, PUBLIC :: radrw = rho0_a/rho0_w !: Density ratio
   REAL(wp), PARAMETER, PUBLIC :: sq_radrw = SQRT(rho0_a/rho0_w)
   REAL(wp), PARAMETER, PUBLIC :: rCp0_w  = 4190._wp    !: Specific heat capacity of seawater (ECMWF 4190) [J/K/kg]
   REAL(wp), PARAMETER, PUBLIC :: rnu0_w  = 1.e-6_wp    !: kinetic viscosity of water                      [m^2/s]
   REAL(wp), PARAMETER, PUBLIC :: rk0_w   = 0.6_wp      !: thermal conductivity of water (at 20C)          [W/m/K]
   !
   REAL(wp), PARAMETER, PUBLIC :: emiss_w = 1._wp       !: Surface emissivity (black-body long-wave radiation) of sea-water []
   !                                                    !: Theoretically close to 0.97! Yet, taken equal as 1 to account for
   !                                                    !: the small fraction of downwelling longwave reflected at the
   !                                                    !: surface (Lind & Katsaros, 1986)
   REAL(wp), PARAMETER, PUBLIC :: rdct_qsat_salt = 0.98_wp  !: reduction factor on specific humidity at saturation (q_sat(T_s)) due to salt
   REAL(wp), PARAMETER, PUBLIC :: rtt0 = 273.16_wp        !: triple point of temperature    [K]
   !
   REAL(wp), PARAMETER, PUBLIC :: rcst_cs = -16._wp*9.80665_wp*rho0_w*rCp0_w*rnu0_w*rnu0_w*rnu0_w/(rk0_w*rk0_w) !: for cool-skin parameterizations... (grav = 9.80665_wp)
   !                              => see eq.(14) in Fairall et al. 1996   (eq.(6) of Zeng aand Beljaars is WRONG! (typo?)


   INTERFACE gamma_moist
      MODULE PROCEDURE gamma_moist_vctr, gamma_moist_sclr
   END INTERFACE gamma_moist

   INTERFACE e_sat
      MODULE PROCEDURE e_sat_vctr, e_sat_sclr
   END INTERFACE e_sat

   INTERFACE L_vap
      MODULE PROCEDURE L_vap_vctr, L_vap_sclr
   END INTERFACE L_vap

   INTERFACE rho_air
      MODULE PROCEDURE rho_air_vctr, rho_air_sclr
   END INTERFACE rho_air

   INTERFACE cp_air
      MODULE PROCEDURE cp_air_vctr, cp_air_sclr
   END INTERFACE cp_air

   INTERFACE alpha_sw
      MODULE PROCEDURE alpha_sw_vctr, alpha_sw_sclr
   END INTERFACE alpha_sw

   INTERFACE bulk_formula
      MODULE PROCEDURE bulk_formula_vctr, bulk_formula_sclr
   END INTERFACE bulk_formula



   PUBLIC virt_temp
   PUBLIC rho_air
   PUBLIC visc_air
   PUBLIC L_vap
   PUBLIC cp_air
   PUBLIC gamma_moist
   PUBLIC One_on_L
   PUBLIC Ri_bulk
   PUBLIC q_sat
   PUBLIC q_air_rh
   !:
   PUBLIC update_qnsol_tau
   PUBLIC alpha_sw
   PUBLIC bulk_formula

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcblk.F90 10535 2019-01-16 17:36:47Z clem $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   FUNCTION virt_temp( pta, pqa )
      !!------------------------------------------------------------------------
      !!
      !! Compute the (absolute/potential) virtual temperature, knowing the
      !! (absolute/potential) temperature and specific humidity
      !!
      !! If input temperature is absolute then output vitual temperature is absolute
      !! If input temperature is potential then output vitual temperature is potential
      !!
      !! Author: L. Brodeau, June 2019 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: virt_temp         !: 1./(Monin Obukhov length) [m^-1]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pta,  &  !: absolute or potetntial air temperature [K]
         &                                        pqa      !: specific humidity of air   [kg/kg]
      !!-------------------------------------------------------------------
      !
      virt_temp(:,:) = pta(:,:) * (1._wp + rctv0*pqa(:,:))
      !!
      !! This is exactly the same sing that:
      !! virt_temp = pta * ( pwa + reps0) / (reps0*(1.+pwa))
      !! with wpa (mixing ration) defined as : pwa = pqa/(1.-pqa)
      !
   END FUNCTION virt_temp

   FUNCTION rho_air_vctr( ptak, pqa, pslp )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION rho_air_vctr  ***
      !!
      !! ** Purpose : compute density of (moist) air using the eq. of state of the atmosphere
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptak      ! air temperature             [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pqa       ! air specific humidity   [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pslp      ! pressure in                [Pa]
      REAL(wp), DIMENSION(jpi,jpj)             ::   rho_air_vctr   ! density of moist air   [kg/m^3]
      !!-------------------------------------------------------------------------------
      rho_air_vctr = MAX( pslp / (R_dry*ptak * ( 1._wp + rctv0*pqa )) , 0.8_wp )
   END FUNCTION rho_air_vctr

   FUNCTION rho_air_sclr( ptak, pqa, pslp )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION rho_air_sclr  ***
      !!
      !! ** Purpose : compute density of (moist) air using the eq. of state of the atmosphere
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!-------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: ptak           ! air temperature             [K]
      REAL(wp), INTENT(in) :: pqa            ! air specific humidity   [kg/kg]
      REAL(wp), INTENT(in) :: pslp           ! pressure in                [Pa]
      REAL(wp)             :: rho_air_sclr   ! density of moist air   [kg/m^3]
      !!-------------------------------------------------------------------------------
      rho_air_sclr = MAX( pslp / (R_dry*ptak * ( 1._wp + rctv0*pqa )) , 0.8_wp )
   END FUNCTION rho_air_sclr



   FUNCTION visc_air(ptak)
      !!----------------------------------------------------------------------------------
      !! Air kinetic viscosity (m^2/s) given from temperature in degrees...
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             ::   visc_air   !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptak       ! air temperature in (K)
      !
      INTEGER  ::   ji, jj      ! dummy loop indices
      REAL(wp) ::   ztc, ztc2   ! local scalar
      !!----------------------------------------------------------------------------------
      !
      DO_2D_11_11
         ztc  = ptak(ji,jj) - rt0   ! air temp, in deg. C
         ztc2 = ztc*ztc
         visc_air(ji,jj) = 1.326e-5*(1. + 6.542E-3*ztc + 8.301e-6*ztc2 - 4.84e-9*ztc2*ztc)
      END_2D
      !
   END FUNCTION visc_air

   FUNCTION L_vap_vctr( psst )
      !!---------------------------------------------------------------------------------
      !!                           ***  FUNCTION L_vap_vctr  ***
      !!
      !! ** Purpose : Compute the latent heat of vaporization of water from temperature
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             ::   L_vap_vctr   ! latent heat of vaporization   [J/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   psst   ! water temperature                [K]
      !!----------------------------------------------------------------------------------
      !
      L_vap_vctr = (  2.501_wp - 0.00237_wp * ( psst(:,:) - rt0)  ) * 1.e6_wp
      !
   END FUNCTION L_vap_vctr

   FUNCTION L_vap_sclr( psst )
      !!---------------------------------------------------------------------------------
      !!                           ***  FUNCTION L_vap_sclr  ***
      !!
      !! ** Purpose : Compute the latent heat of vaporization of water from temperature
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp)             ::   L_vap_sclr   ! latent heat of vaporization   [J/kg]
      REAL(wp), INTENT(in) ::   psst         ! water temperature                [K]
      !!----------------------------------------------------------------------------------
      !
      L_vap_sclr = (  2.501_wp - 0.00237_wp * ( psst - rt0)  ) * 1.e6_wp
      !
   END FUNCTION L_vap_sclr


   FUNCTION cp_air_vctr( pqa )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION cp_air_vctr  ***
      !!
      !! ** Purpose : Compute specific heat (Cp) of moist air
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!-------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pqa      ! air specific humidity         [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj)             ::   cp_air_vctr   ! specific heat of moist air   [J/K/kg]
      !!-------------------------------------------------------------------------------
      cp_air_vctr = rCp_dry + rCp_vap * pqa
   END FUNCTION cp_air_vctr

   FUNCTION cp_air_sclr( pqa )
      !!-------------------------------------------------------------------------------
      !!                           ***  FUNCTION cp_air_sclr  ***
      !!
      !! ** Purpose : Compute specific heat (Cp) of moist air
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!-------------------------------------------------------------------------------
      REAL(wp), INTENT(in) :: pqa           ! air specific humidity         [kg/kg]
      REAL(wp)             :: cp_air_sclr   ! specific heat of moist air   [J/K/kg]
      !!-------------------------------------------------------------------------------
      cp_air_sclr = rCp_dry + rCp_vap * pqa
   END FUNCTION cp_air_sclr





   FUNCTION gamma_moist_vctr( ptak, pqa )
      !!----------------------------------------------------------------------------------
      !!                           ***  FUNCTION gamma_moist_vctr  ***
      !!
      !! ** Purpose : Compute the moist adiabatic lapse-rate.
      !!     => http://glossary.ametsoc.org/wiki/Moist-adiabatic_lapse_rate
      !!     => http://www.geog.ucsb.edu/~joel/g266_s10/lecture_notes/chapt03/oh10_3_01/oh10_3_01.html
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptak          ! air temperature       [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pqa           ! specific humidity [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj)             ::   gamma_moist_vctr   ! moist adiabatic lapse-rate
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      !!----------------------------------------------------------------------------------
      DO_2D_11_11
         gamma_moist_vctr(ji,jj) = gamma_moist_sclr( ptak(ji,jj), pqa(ji,jj) )
      END_2D
   END FUNCTION gamma_moist_vctr

   FUNCTION gamma_moist_sclr( ptak, pqa )
      !!----------------------------------------------------------------------------------
      !! ** Purpose : Compute the moist adiabatic lapse-rate.
      !!     => http://glossary.ametsoc.org/wiki/Moist-adiabatic_lapse_rate
      !!     => http://www.geog.ucsb.edu/~joel/g266_s10/lecture_notes/chapt03/oh10_3_01/oh10_3_01.html
      !!
      !! ** Author: L. Brodeau, June 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp)             :: gamma_moist_sclr
      REAL(wp), INTENT(in) :: ptak, pqa ! air temperature (K) and specific humidity (kg/kg)
      !
      REAL(wp) :: zta, zqa, zwa, ziRT        ! local scalar
      !!----------------------------------------------------------------------------------
      zta = MAX( ptak,  180._wp) ! prevents screw-up over masked regions where field == 0.
      zqa = MAX( pqa,  1.E-6_wp) !    "                   "                     "
      !!
      zwa = zqa / (1._wp - zqa)   ! w is mixing ratio w = q/(1-q) | q = w/(1+w)
      ziRT = 1._wp / (R_dry*zta)    ! 1/RT
      gamma_moist_sclr = grav * ( 1._wp + rLevap*zwa*ziRT ) / ( rCp_dry + rLevap*rLevap*zwa*reps0*ziRT/zta )
      !!
   END FUNCTION gamma_moist_sclr

   FUNCTION One_on_L( ptha, pqa, pus, pts, pqs )
      !!------------------------------------------------------------------------
      !!
      !! Evaluates the 1./(Monin Obukhov length) from air temperature and
      !!  specific humidity, and frictional scales u*, t* and q*
      !!
      !! Author: L. Brodeau, June 2016 / AeroBulk
      !!         (https://github.com/brodeau/aerobulk/)
      !!------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: One_on_L         !: 1./(Monin Obukhov length) [m^-1]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ptha,  &  !: average potetntial air temperature [K]
         &                                        pqa,   &  !: average specific humidity of air   [kg/kg]
         &                                      pus, pts, pqs   !: frictional velocity, temperature and humidity
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::     zqa          ! local scalar
      !!-------------------------------------------------------------------
      !
      DO_2D_11_11
         !
         zqa = (1._wp + rctv0*pqa(ji,jj))
         !
         ! The main concern is to know whether, the vertical turbulent flux of virtual temperature, < u' theta_v' > is estimated with:
         !  a/  -u* [ theta* (1 + 0.61 q) + 0.61 theta q* ] => this is the one that seems correct! chose this one!
         !                      or
         !  b/  -u* [ theta*              + 0.61 theta q* ]
         !
         One_on_L(ji,jj) = grav*vkarmn*( pts(ji,jj)*zqa + rctv0*ptha(ji,jj)*pqs(ji,jj) ) &
            &               / MAX( pus(ji,jj)*pus(ji,jj)*ptha(ji,jj)*zqa , 1.E-9_wp )
         !
      END_2D
      !
      One_on_L = SIGN( MIN(ABS(One_on_L),200._wp), One_on_L ) ! (prevent FPE from stupid values over masked regions...)
      !
   END FUNCTION One_on_L

   FUNCTION Ri_bulk( pz, psst, ptha, pssq, pqa, pub )
      !!----------------------------------------------------------------------------------
      !! Bulk Richardson number according to "wide-spread equation"...
      !!
      !! ** Author: L. Brodeau, June 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: Ri_bulk
      REAL(wp)                    , INTENT(in) :: pz    ! height above the sea (aka "delta z")  [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psst  ! SST                                   [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ptha  ! pot. air temp. at height "pz"         [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pssq  ! 0.98*q_sat(SST)                   [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pqa   ! air spec. hum. at height "pz"     [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pub   ! bulk wind speed                     [m/s]
      !
      INTEGER  ::   ji, jj                                ! dummy loop indices
      REAL(wp) ::   zqa, zta, zgamma, zdth_v, ztv, zsstv  ! local scalars
      !!-------------------------------------------------------------------
      !
      DO_2D_11_11
         !
         zqa = 0.5_wp*(pqa(ji,jj)+pssq(ji,jj))                                        ! ~ mean q within the layer...
         zta = 0.5_wp*( psst(ji,jj) + ptha(ji,jj) - gamma_moist(ptha(ji,jj),zqa)*pz ) ! ~ mean absolute temperature of air within the layer
         zta = 0.5_wp*( psst(ji,jj) + ptha(ji,jj) - gamma_moist(zta,        zqa)*pz ) ! ~ mean absolute temperature of air within the layer
         zgamma =  gamma_moist(zta, zqa)                                              ! Adiabatic lapse-rate for moist air within the layer
         !
         zsstv = psst(ji,jj)*(1._wp + rctv0*pssq(ji,jj)) ! absolute==potential virtual SST (absolute==potential because z=0!)
         !
         zdth_v = ptha(ji,jj)*(1._wp + rctv0*pqa(ji,jj)) - zsstv ! air-sea delta of "virtual potential temperature"
         !
         ztv = 0.5_wp*( zsstv + (ptha(ji,jj) - zgamma*pz)*(1._wp + rctv0*pqa(ji,jj)) )  ! ~ mean absolute virtual temp. within the layer
         !
         Ri_bulk(ji,jj) = grav*zdth_v*pz / ( ztv*pub(ji,jj)*pub(ji,jj) )                            ! the usual definition of Ri_bulk
         !
      END_2D
   END FUNCTION Ri_bulk


   FUNCTION e_sat_vctr(ptak)
      !!**************************************************
      !! ptak:     air temperature [K]
      !! e_sat:  water vapor at saturation [Pa]
      !!
      !! Recommended by WMO
      !!
      !! Goff, J. A., 1957: Saturation pressure of water on the new kelvin
      !! temperature scale. Transactions of the American society of heating
      !! and ventilating engineers, 347–354.
      !!
      !! rt0 should be 273.16 (triple point of water) and not 273.15 like here
      !!**************************************************

      REAL(wp), DIMENSION(jpi,jpj)             :: e_sat_vctr !: vapour pressure at saturation  [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ptak    !: temperature (K)

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: ztmp

      ALLOCATE ( ztmp(jpi,jpj) )

      ztmp(:,:) = rtt0/ptak(:,:)

      e_sat_vctr = 100.*( 10.**(10.79574*(1. - ztmp) - 5.028*LOG10(ptak/rtt0)         &
         &       + 1.50475*10.**(-4)*(1. - 10.**(-8.2969*(ptak/rtt0 - 1.)) )   &
         &       + 0.42873*10.**(-3)*(10.**(4.76955*(1. - ztmp)) - 1.) + 0.78614) )

      DEALLOCATE ( ztmp )

   END FUNCTION e_sat_vctr


   FUNCTION e_sat_sclr( ptak )
      !!----------------------------------------------------------------------------------
      !!                   ***  FUNCTION e_sat_sclr  ***
      !!                  < SCALAR argument version >
      !! ** Purpose : water vapor at saturation in [Pa]
      !!              Based on accurate estimate by Goff, 1957
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!
      !!    Note: what rt0 should be here, is 273.16 (triple point of water) and not 273.15 like here
      !!----------------------------------------------------------------------------------
      REAL(wp), INTENT(in) ::   ptak    ! air temperature                  [K]
      REAL(wp)             ::   e_sat_sclr   ! water vapor at saturation   [kg/kg]
      !
      REAL(wp) ::   zta, ztmp   ! local scalar
      !!----------------------------------------------------------------------------------
      !
      zta = MAX( ptak , 180._wp )   ! air temp., prevents fpe0 errors dute to unrealistically low values over masked regions...
      ztmp = rt0 / zta
      !
      ! Vapour pressure at saturation [Pa] : WMO, (Goff, 1957)
      e_sat_sclr = 100.*( 10.**( 10.79574*(1. - ztmp) - 5.028*LOG10(zta/rt0)        &
         &    + 1.50475*10.**(-4)*(1. - 10.**(-8.2969*(zta/rt0 - 1.)) )  &
         &    + 0.42873*10.**(-3)*(10.**(4.76955*(1. - ztmp)) - 1.) + 0.78614) )
      !
   END FUNCTION e_sat_sclr


   FUNCTION q_sat( ptak, pslp )
      !!----------------------------------------------------------------------------------
      !!                           ***  FUNCTION q_sat  ***
      !!
      !! ** Purpose : Specific humidity at saturation in [kg/kg]
      !!              Based on accurate estimate of "e_sat"
      !!              aka saturation water vapor (Goff, 1957)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   ptak    ! air temperature                       [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   pslp    ! sea level atmospheric pressure       [Pa]
      REAL(wp), DIMENSION(jpi,jpj)             ::   q_sat   ! Specific humidity at saturation   [kg/kg]
      !
      INTEGER  ::   ji, jj         ! dummy loop indices
      REAL(wp) ::   ze_sat   ! local scalar
      !!----------------------------------------------------------------------------------
      !
      DO_2D_11_11
         !
         ze_sat =  e_sat_sclr( ptak(ji,jj) )
         !
         q_sat(ji,jj) = reps0 * ze_sat/( pslp(ji,jj) - (1._wp - reps0)*ze_sat )
         !
      END_2D
      !
   END FUNCTION q_sat

   FUNCTION q_air_rh(prha, ptak, pslp)
      !!----------------------------------------------------------------------------------
      !! Specific humidity of air out of Relative Humidity
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             :: q_air_rh
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: prha        !: relative humidity      [fraction, not %!!!]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: ptak        !: air temperature        [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pslp        !: atmospheric pressure   [Pa]
      !
      INTEGER  ::   ji, jj      ! dummy loop indices
      REAL(wp) ::   ze      ! local scalar
      !!----------------------------------------------------------------------------------
      !
      DO_2D_11_11
         ze = prha(ji,jj)*e_sat_sclr(ptak(ji,jj))
         q_air_rh(ji,jj) = ze*reps0/(pslp(ji,jj) - (1. - reps0)*ze)
      END_2D
      !
   END FUNCTION q_air_rh


   SUBROUTINE UPDATE_QNSOL_TAU( pzu, pTs, pqs, pTa, pqa, pust, ptst, pqst, pwnd, pUb, pslp, prlw, &
      &                         pQns, pTau,  &
      &                         Qlat)
      !!----------------------------------------------------------------------------------
      !! Purpose: returns the non-solar heat flux to the ocean aka "Qlat + Qsen + Qlw"
      !!          and the module of the wind stress => pTau = Tau
      !! ** Author: L. Brodeau, Sept. 2019 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp),                     INTENT(in)  :: pzu  ! height above the sea-level where all this takes place (normally 10m)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pTs  ! water temperature at the air-sea interface [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pqs  ! satur. spec. hum. at T=pTs   [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pTa  ! potential air temperature at z=pzu [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pqa  ! specific humidity at z=pzu [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pust ! u*
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: ptst ! t*
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pqst ! q*
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pwnd ! wind speed module at z=pzu [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pUb  ! bulk wind speed at z=pzu (inc. pot. effect of gustiness etc) [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pslp ! sea-level atmospheric pressure [Pa]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: prlw ! downwelling longwave radiative flux [W/m^2]
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pQns ! non-solar heat flux to the ocean aka "Qlat + Qsen + Qlw" [W/m^2]]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pTau ! module of the wind stress [N/m^2]
      !
      REAL(wp), DIMENSION(jpi,jpj), OPTIONAL, INTENT(out) :: Qlat
      !
      REAL(wp) :: zdt, zdq, zCd, zCh, zCe, zTs2, zz0, &
         &        zQlat, zQsen, zQlw
      INTEGER  ::   ji, jj     ! dummy loop indices
      !!----------------------------------------------------------------------------------
      DO_2D_11_11

         zdt = pTa(ji,jj) - pTs(ji,jj) ;  zdt = SIGN( MAX(ABS(zdt),1.E-6_wp), zdt )
         zdq = pqa(ji,jj) - pqs(ji,jj) ;  zdq = SIGN( MAX(ABS(zdq),1.E-9_wp), zdq )
         zz0 = pust(ji,jj)/pUb(ji,jj)
         zCd = zz0*zz0
         zCh = zz0*ptst(ji,jj)/zdt
         zCe = zz0*pqst(ji,jj)/zdq

         CALL BULK_FORMULA_SCLR( pzu, pTs(ji,jj), pqs(ji,jj), pTa(ji,jj), pqa(ji,jj), &
            &                    zCd, zCh, zCe,                                       &
            &                    pwnd(ji,jj), pUb(ji,jj), pslp(ji,jj),                &
            &                    pTau(ji,jj), zQsen, zQlat )
         
         zTs2  = pTs(ji,jj)*pTs(ji,jj)
         zQlw  = emiss_w*(prlw(ji,jj) - stefan*zTs2*zTs2) ! Net longwave flux

         pQns(ji,jj) = zQlat + zQsen + zQlw

         IF( PRESENT(Qlat) ) Qlat(ji,jj) = zQlat
      END_2D
   END SUBROUTINE UPDATE_QNSOL_TAU


   SUBROUTINE BULK_FORMULA_SCLR( pzu, pTs, pqs, pTa, pqa, &
      &                          pCd, pCh, pCe,           &
      &                          pwnd, pUb, pslp,         &
      &                          pTau, pQsen, pQlat,      &
      &                          pEvap, prhoa, pfact_evap )
      !!----------------------------------------------------------------------------------
      REAL(wp),                     INTENT(in)  :: pzu  ! height above the sea-level where all this takes place (normally 10m)
      REAL(wp), INTENT(in)  :: pTs  ! water temperature at the air-sea interface [K]
      REAL(wp), INTENT(in)  :: pqs  ! satur. spec. hum. at T=pTs   [kg/kg]
      REAL(wp), INTENT(in)  :: pTa  ! potential air temperature at z=pzu [K]
      REAL(wp), INTENT(in)  :: pqa  ! specific humidity at z=pzu [kg/kg]
      REAL(wp), INTENT(in)  :: pCd
      REAL(wp), INTENT(in)  :: pCh
      REAL(wp), INTENT(in)  :: pCe
      REAL(wp), INTENT(in)  :: pwnd ! wind speed module at z=pzu [m/s]
      REAL(wp), INTENT(in)  :: pUb  ! bulk wind speed at z=pzu (inc. pot. effect of gustiness etc) [m/s]
      REAL(wp), INTENT(in)  :: pslp ! sea-level atmospheric pressure [Pa]
      !!
      REAL(wp), INTENT(out) :: pTau  ! module of the wind stress [N/m^2]
      REAL(wp), INTENT(out) :: pQsen !  [W/m^2]
      REAL(wp), INTENT(out) :: pQlat !  [W/m^2]
      !!
      REAL(wp), INTENT(out), OPTIONAL :: pEvap ! Evaporation [kg/m^2/s]
      REAL(wp), INTENT(out), OPTIONAL :: prhoa ! Air density at z=pzu [kg/m^3]
      REAL(wp), INTENT(in) , OPTIONAL :: pfact_evap  ! ABOMINATION: corrective factor for evaporation (doing this against my will! /laurent)
      !!
      REAL(wp) :: ztaa, zgamma, zrho, zUrho, zevap, zfact_evap
      INTEGER  :: jq
      !!----------------------------------------------------------------------------------
      zfact_evap = 1._wp
      IF( PRESENT(pfact_evap) ) zfact_evap = pfact_evap
      
      !! Need ztaa, absolute temperature at pzu (formula to estimate rho_air needs absolute temperature, not the potential temperature "pTa")
      ztaa = pTa ! first guess...
      DO jq = 1, 4
         zgamma = gamma_moist( 0.5*(ztaa+pTs) , pqa )  !LOLO: why not "0.5*(pqs+pqa)" rather then "pqa" ???
         ztaa = pTa - zgamma*pzu   ! Absolute temp. is slightly colder...
      END DO
      zrho = rho_air(ztaa, pqa, pslp)
      zrho = rho_air(ztaa, pqa, pslp-zrho*grav*pzu) ! taking into account that we are pzu m above the sea level where SLP is given!

      zUrho = pUb*MAX(zrho, 1._wp)     ! rho*U10

      pTau = zUrho * pCd * pwnd ! Wind stress module

      zevap = zUrho * pCe * (pqa - pqs)
      pQsen = zUrho * pCh * (pTa - pTs) * cp_air(pqa)
      pQlat = L_vap(pTs) * zevap

      IF( PRESENT(pEvap) ) pEvap = - zfact_evap * zevap
      IF( PRESENT(prhoa) ) prhoa = zrho

   END SUBROUTINE BULK_FORMULA_SCLR

   SUBROUTINE BULK_FORMULA_VCTR( pzu, pTs, pqs, pTa, pqa, &
      &                          pCd, pCh, pCe,           &
      &                          pwnd, pUb, pslp,         &
      &                          pTau, pQsen, pQlat,      & 
      &                          pEvap, prhoa, pfact_evap )      
      !!----------------------------------------------------------------------------------
      REAL(wp),                     INTENT(in)  :: pzu  ! height above the sea-level where all this takes place (normally 10m)
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pTs  ! water temperature at the air-sea interface [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pqs  ! satur. spec. hum. at T=pTs   [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pTa  ! potential air temperature at z=pzu [K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pqa  ! specific humidity at z=pzu [kg/kg]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pCd
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pCh
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pCe
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pwnd ! wind speed module at z=pzu [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pUb  ! bulk wind speed at z=pzu (inc. pot. effect of gustiness etc) [m/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: pslp ! sea-level atmospheric pressure [Pa]
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pTau  ! module of the wind stress [N/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pQsen !  [W/m^2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: pQlat !  [W/m^2]
      !!
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL :: pEvap ! Evaporation [kg/m^2/s]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out), OPTIONAL :: prhoa ! Air density at z=pzu [kg/m^3]
      REAL(wp),                     INTENT(in) , OPTIONAL :: pfact_evap  ! ABOMINATION: corrective factor for evaporation (doing this against my will! /laurent)
      !!
      REAL(wp) :: ztaa, zgamma, zrho, zUrho, zevap, zfact_evap
      INTEGER  :: ji, jj
      !!----------------------------------------------------------------------------------
      zfact_evap = 1._wp
      IF( PRESENT(pfact_evap) ) zfact_evap = pfact_evap

      DO_2D_11_11

         CALL BULK_FORMULA_SCLR( pzu, pTs(ji,jj), pqs(ji,jj), pTa(ji,jj), pqa(ji,jj), &
            &                    pCd(ji,jj), pCh(ji,jj), pCe(ji,jj),                  &
            &                    pwnd(ji,jj), pUb(ji,jj), pslp(ji,jj),                &
            &                    pTau(ji,jj), pQsen(ji,jj), pQlat(ji,jj),             &
            &                    pEvap=zevap, prhoa=zrho, pfact_evap=zfact_evap       )

         IF( PRESENT(pEvap) ) pEvap(ji,jj) = zevap
         IF( PRESENT(prhoa) ) prhoa(ji,jj) = zrho
   
      END_2D
   END SUBROUTINE BULK_FORMULA_VCTR


   FUNCTION alpha_sw_vctr( psst )
      !!---------------------------------------------------------------------------------
      !!                           ***  FUNCTION alpha_sw_vctr  ***
      !!
      !! ** Purpose : ROUGH estimate of the thermal expansion coefficient of sea-water at the surface (P =~ 1010 hpa)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)             ::   alpha_sw_vctr   ! thermal expansion coefficient of sea-water [1/K]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) ::   psst   ! water temperature                [K]
      !!----------------------------------------------------------------------------------
      alpha_sw_vctr = 2.1e-5_wp * MAX(psst(:,:)-rt0 + 3.2_wp, 0._wp)**0.79
   END FUNCTION alpha_sw_vctr

   FUNCTION alpha_sw_sclr( psst )
      !!---------------------------------------------------------------------------------
      !!                           ***  FUNCTION alpha_sw_sclr  ***
      !!
      !! ** Purpose : ROUGH estimate of the thermal expansion coefficient of sea-water at the surface (P =~ 1010 hpa)
      !!
      !! ** Author: L. Brodeau, june 2016 / AeroBulk (https://github.com/brodeau/aerobulk/)
      !!----------------------------------------------------------------------------------
      REAL(wp)             ::   alpha_sw_sclr   ! thermal expansion coefficient of sea-water [1/K]
      REAL(wp), INTENT(in) ::   psst   ! sea-water temperature                   [K]
      !!----------------------------------------------------------------------------------
      alpha_sw_sclr = 2.1e-5_wp * MAX(psst-rt0 + 3.2_wp, 0._wp)**0.79
   END FUNCTION alpha_sw_sclr



   !!======================================================================
END MODULE sbcblk_phy
