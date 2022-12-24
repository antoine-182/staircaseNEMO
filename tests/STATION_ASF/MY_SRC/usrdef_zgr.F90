MODULE usrdef_zgr
   !!======================================================================
   !!                     ***  MODULE usrdef_zgr  ***
   !!
   !!                       ===  STATION_ASF case  ===
   !!
   !! user defined :  vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2019-10  (L. Brodeau)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system (required)
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   !USE dom_oce        ! ocean domain
   !USE depth_e3       ! depth <=> e3
   USE usrdef_nam     ! User defined : namelist variables
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr   ! called by domzgr.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 10072 2018-08-28 15:21:50Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_zgr( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
      &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,    &   ! 1D reference vertical coordinate
      &                    pdept , pdepw ,                             &   ! 3D t & w-points depth
      &                    pe3t  , pe3u  , pe3v   , pe3f ,             &   ! vertical scale factors
      &                    pe3w  , pe3uw , pe3vw         ,             &   !     -      -      -
      &                    k_top  , k_bot    )                             ! top & bottom ocean level
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(  out) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags ( read in namusr_def )
      LOGICAL                   , INTENT(  out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         ! i-scale factors
      INTEGER , DIMENSION(:,:)  , INTENT(  out) ::   k_top, k_bot                ! first & last ocean level
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : STATION_ASF configuration, setting first level properties.'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !

      ld_zco    = .TRUE.         ! z-coordinate without ocean cavities
      ld_zps    = .FALSE.
      ld_sco    = .FALSE.
      ld_isfcav = .FALSE.
      
      pdept_1d(1) = rn_dept1 ! depth (m) at which the SST is taken/measured == depth of first T point!
      pdepw_1d(1) = 0._wp
      pe3t_1d(1)  = 2._wp*rn_dept1
      pe3w_1d(1)  = rn_dept1 ! LB???

      pdept(:,:,:) = rn_dept1
      pdepw(:,:,:) = 0._wp
      pe3t(:,:,:) = 2._wp*rn_dept1
      pe3u(:,:,:) = 2._wp*rn_dept1
      pe3v(:,:,:) = 2._wp*rn_dept1
      pe3f(:,:,:) = 2._wp*rn_dept1
      pe3w(:,:,:)  = rn_dept1  ! LB???
      pe3uw(:,:,:) = rn_dept1  ! LB???
      pe3vw(:,:,:) = rn_dept1  ! LB???
      k_top = 1
      k_bot = 1
      !
   END SUBROUTINE usr_def_zgr
   !!======================================================================
END MODULE usrdef_zgr
