MODULE usrdef_istate
   !!======================================================================
   !!                   ***  MODULE  usrdef_istate   ***
   !!
   !!                     ===  AM98 configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  4.0 ! 2016-03  (S. Flavoni) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE dom_oce , ONLY : glamt, gphit,gphiv, glamu, r1_e1u, r1_e2v

   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called in istate.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !!
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here AM98 configuration example : (double gyre with rotated domain)
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity   field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      INTEGER :: ji, jj, jk  ! dummy loop indices
      REAL(wp) :: zr, zr2, zr0, zr02, z0x, z0y, zgamma, zv_azimuth
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : analytical definition of initial state '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest, with an horizontally uniform T and S profiles'
      !
      pu  (:,:,:) = 0._wp        ! ocean at rest
      pv  (:,:,:) = 0._wp
      pssh(:,:)   = 0._wp
!!an
      zgamma = -5.e4  ! m2/s
      zr0 = 20000._wp  ! 20km
      ! zr0 = 200000._wp  ! 20km
      zr02 = zr0*zr0  !
      z0x = 150000._wp ; z0y = 150000._wp ! center of the vortex is (150km,150km)
      ! z0x = 250000._wp ; z0y = 250000._wp ! center of the vortex is (150km,150km)
      DO ji=1,jpi
        DO jj=1,jpj
          zr2 = ( glamt(ji,jj) - z0x )**2 + ( gphit(ji,jj) - z0y )**2
          zr  = SQRT(zr2)
          zv_azimuth =  (zgamma/(2._wp*rpi*zr)) * (1._wp - EXP(-zr2/(2._wp*zr02)) )
          IF ( gphiv(ji,jj) >= (250000._wp + 0.1_wp) .OR. glamu(ji,jj) >= (100000._wp + 0.1_wp) )  THEN
            !
            pu  (ji,jj,1) = - zv_azimuth * ( gphit(ji,jj)- z0y )/zr    !  v_azimuthal * sin(theta)
            pv  (ji,jj,1) =   zv_azimuth * ( glamt(ji,jj)- z0x )/zr    !  v_azimuthal * cos(theta)
            !
            ! pssh(ji,jj  ) = ( zv_azimuth*zv_azimuth )/ (grav*zr)
          END IF
        END DO
      END DO

      !
      ! vortex in a f plan (200m)
      ! DO ji=1,jpi
      !   DO jj=1,jpj
      !     zr2 = ( glamt(ji,jj)- z0x )**2 + ( gphit(ji,jj)- z0y )**2
      !     pssh(ji,jj) = 1._wp * EXP(-zr2 / (2._wp*zr02))
      !   END DO
      ! END DO

      !
      ! DO ji=1,jpim1
      !   DO jj=1,jpjm1
      !     pu(ji,jj,1) = r1_e1u * (pssh(ji+1,jj  ) - pssh(ji,jj))    !
      !     pv(ji,jj,1) = r1_e2v * (pssh(ji  ,jj+1) - pssh(ji,jj))    !
      !   END DO
      ! END DO
!!an
       IF(lwp) WRITE(numout,*) "end istate"
        call flush(numout)
      !
   END SUBROUTINE usr_def_istate

   !!======================================================================
END MODULE usrdef_istate
