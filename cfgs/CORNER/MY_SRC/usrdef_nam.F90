MODULE usrdef_nam
   !!======================================================================
   !!                     ***  MODULE usrdef_nam   ***
   !!
   !!                     ===  AM98 configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp, njmpp       ! ocean space and time domain
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called in nemogcm.F90 module
   !                              !!* namusr_def namelist *!!
   INTEGER , PUBLIC ::   nn_AM98            ! 1/nn_AM98 = the resolution chosen in degrees and thus defining the horizontal domain size
   REAL(wp), PUBLIC ::   rn_theta           ! rotation angle (in degree) of the grid
   INTEGER , PUBLIC ::   nn_gc              ! number of ghostcells
   REAL(wp), PUBLIC ::   rn_domsiz          ! size of the domain (default 2000km)  [m]
   REAL(wp), PUBLIC ::   rn_dx              ! gridspacing (default 100km)          [m]
   REAL(wp), PUBLIC ::   rn_tau             ! wind stress on the surface           [N/m2]
   REAL(wp), PUBLIC ::   rn_f0              !    f-plan coriolis parameter         [1/s]
   REAL(wp), PUBLIC ::   rn_beta            ! beta-plan coriolis parameter         [1/m.s]
   REAL(wp), PUBLIC ::   rn_modified_grav   ! modified gravity                     [m/s2]
   REAL(wp), PUBLIC ::   rn_rfr             ! layer friction                       [1/s]
   !                              !!* temporary *!!
   INTEGER , PUBLIC ::   nn_dynldf_lap_typ            ! choose type of laplacian (ideally from namelist)
   !                                                  ! = 1   divrot    laplacian
   !                                                  ! = 2   symmetric laplacian (Griffies&Hallberg 2000)
   !                                                  ! = 3   symmetric laplacian (cartesian)
   LOGICAL , PUBLIC ::   ln_dynldf_lap_PM             ! if True - apply the P.Marchand boundary condition on the laplacian
   !
   LOGICAL , PUBLIC ::   ln_hdiv_AD          ! Layer width diffusion on height
   !
   LOGICAL , PUBLIC ::   ln_obstacle         ! Keep the 1°/4 staircase shape
   !
   LOGICAL , PUBLIC ::   ln_hldc             ! Layer width diffusion on height
   REAL(wp), PUBLIC ::   rn_hldc             ! friction                       [m2/s]
   !                              !!* penalisation *!!
   REAL(wp), PUBLIC ::   rn_abp             ! alpha boundary parameter                                       [-]
   REAL(wp), PUBLIC ::   rn_cnp             ! size of the penalised domain (in number of cells)
   INTEGER , PUBLIC ::   nn_smo             ! smoothing parameters (X x Y shapiro filters applied nn_smo times)
   REAL(wp), PUBLIC ::   rn_fsp             ! friction parameter 1/epsilon of the permeability               [1/s]
   !
   INTEGER , PUBLIC ::   nn_rfr             ! layer friction for penalisation                      [1/s]
   INTEGER , PUBLIC ::   nn_fsp             ! friction parameter type                              [1/s]
   !
   REAL(wp), PUBLIC ::   r1_abp             ! inverse alpha boundary parameter                            [-]
   !
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_nam( cd_cfg, kk_cfg, kpi, kpj, kpk, kperio )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!
      !! ** Purpose :   read user defined namelist and define the domain size
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here AM98 configuration
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(out) ::   cd_cfg          ! configuration name
      INTEGER         , INTENT(out) ::   kk_cfg          ! configuration resolution
      INTEGER         , INTENT(out) ::   kpi, kpj, kpk   ! global domain sizes
      INTEGER         , INTENT(out) ::   kperio          ! lateral global domain b.c.
      !
      INTEGER  ::   ios             ! Local integer
      REAL(wp) ::   ze1, zgl, zbl   ! gridspacing, length of the biggest square
      REAL(wp) ::   zfac, zek ,zei0, znI, znq, zel, zei
      REAL(wp) ::   ze0, zgl0, zbl0   ! ln_obstacle
      !!
      NAMELIST/namusr_def/ nn_AM98, rn_theta, jpkglo,           &   !
         &                 nn_gc ,rn_domsiz, rn_dx,             &   ! domain parameters
         &                 ln_obstacle,                         &   ! size of staircase
         &                 rn_f0 ,rn_beta, rn_tau,              &   ! coriolis parameter, wind
         &                 rn_modified_grav, rn_rfr, nn_rfr,    &   ! reduced gravity, friction
         &                 ln_hldc, rn_hldc, ln_hdiv_AD,       &   ! layer width diffusion
         &                 nn_dynldf_lap_typ, ln_dynldf_lap_PM, &   ! temporary parameter
         &                 rn_abp, rn_cnp, nn_smo, rn_fsp, nn_fsp           ! penalisation parameters
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )
      !
      cd_cfg = 'AM98'               ! name & resolution (not used)

#if defined key_agrif
      IF (.NOT.Agrif_root()) nn_AM98 = Agrif_parent(nn_AM98) * Agrif_irhox()
#endif

      kk_cfg = nn_AM98
      !
      ze1 =  rn_dx / REAL(nn_AM98, wp)                   ! [m] gridspacing used
      zgl =  rn_domsiz + (4._wp + 2._wp * REAL(nn_gc, wp) ) * ze1   ! [m] length of the square with ghostcells
      ! rotation
      zbl = zgl * ( COS( rn_theta * rad ) + SIN( rn_theta * rad ) )   ! length side bigger domain [m]
      !
        !
#if !defined key_refineE && defined key_refineW
        ! work only for theta = 0
        ! glam
        znI = 21._wp ; znq = 0.9_wp
        zel = 0._wp ; zei = -1._wp ; zei0 = 2._wp + REAL(nn_gc, wp) ; zfac = 1._wp
        DO WHILE( zel < zbl ) !
          zei = zei + 1
          zek = zei - zei0
          IF (zek <   0._wp                 ) zfac = 1._wp
          IF (zek >=  0._wp .AND. zek < znI ) zfac = znq**(znI - zek)
          IF (zek >=     znI                 ) zfac = 1._wp
          zel = zel + ze1 * zfac
        END DO
        kpi = zei
#elif !defined key_refineW && defined key_refineE
        ! work only for theta = 0
        ! kpi = ceiling(zbl / ze1 )    ! Global Domain size + ghost cells
        znI = 21._wp ; znq = 0.9_wp  ! znI set with znq so the grispace 10 times smaller
        zel = 0._wp ; zei = -1._wp ; zei0 = nint(zbl/(ze1*2)) ; zfac = 1._wp
        DO WHILE( zel < zbl ) !
          zei = zei + 1
          zek = zei - zei0
          IF (zek <   0._wp                 ) zfac = 1._wp
          IF (zek >=  0._wp .AND. zek < znI ) zfac = znq**(zek)
          IF (zek >=    znI                 ) zfac = znq**(znI)
          zel = zel + ze1 * zfac
        END DO
        kpi = zei
#elif defined key_refineWlin
      ! work only for theta = 0
      ! kpi = ceiling(zbl / ze1 )    ! Global Domain size + ghost cells
      zel = - (2._wp + REAL(nn_gc, wp)) * ze1 ; zei = -1._wp ; zfac = 1._wp
      DO WHILE( zel < zbl ) ! tant que la longueur cumulée < 2000km projeté
        zei = zei + 1
        ! 25km est retiré car le cumul de la zone pénalisée est alpha/2 fois plus grande
        zfac = (1._wp - rn_abp)*(zel - (0._wp + (rn_cnp*ze1/2._wp) ))/(rn_cnp * ze1) + 0.5_wp*(1._wp + rn_abp)
        zfac = MAX(MIN(zfac,1._wp), rn_abp)
        !
        zel = zel + ze1 * zfac
      END DO
      kpi = zei
#else
      kpi = ceiling(zbl / ze1 )    ! Global Domain size + ghost cells
#endif

#if defined key_refineN
      ! work only for theta = 0
      ! kpi = ceiling(zbl / ze1 )    ! Global Domain size + ghost cells
      znI = 21._wp ; znq = 0.9_wp  ! znI set with znq so the grispace 10 times smaller
      zel = 0._wp ; zei = -1._wp ; zei0 = nint(3._wp*zbl/(ze1*4._wp)) ; zfac = 1._wp
      DO WHILE( zel < zbl ) !
        zei = zei + 1
        zek = zei - zei0
        IF (zek <   0._wp                 ) zfac = 1._wp
        IF (zek >=  0._wp .AND. zek < znI ) zfac = znq**(zek)
        IF (zek >=    znI                 ) zfac = znq**(znI)
        zel = zel + ze1 * zfac
      END DO
      kpj = zei
#elif defined key_refineNlin
      ! work only for theta = 0
      ! kpi = ceiling(zbl / ze1 )    ! Global Domain size + ghost cells
      zel = - (2._wp + REAL(nn_gc, wp)) * ze1 ; zei = -1._wp ; zfac = 1._wp
      DO WHILE( zel < zbl ) ! tant que la longueur cumulée < 2000km projeté
        zei = zei + 1
        ! 25km est retiré car le cumul de la zone pénalisée est alpha/2 fois plus grande
        zfac = -(1._wp - rn_abp)*(zel - (2000000._wp - (rn_cnp*ze1/4._wp) ))/(rn_cnp * ze1) + 0.5_wp*(1._wp + rn_abp)
        zfac = MAX(MIN(zfac,1._wp), rn_abp)
        !
        zel = zel + ze1 * zfac
      END DO
      kpj = zei
#elif defined key_refineR
      ! adapt dx so dx = R/2
      ! work only for theta = 0
      ! kpi = ceiling(zbl / ze1 )    ! Global Domain size + ghost cells
      zel = - (2._wp + REAL(nn_gc, wp)) * ze1 ; zei = -1._wp ; zfac = 1._wp
      DO WHILE( zel < zbl ) ! tant que la longueur cumulée < 2000km projeté
        zei = zei + 1
        ze1 = 0.5_wp * SQRT(rn_modified_grav*500._wp) / (rn_f0 + rn_beta*zel)
        !
        zel = zel + ze1
      END DO
      kpj = zei
#elif defined key_onecell
      kpj = ceiling(zbl / ze1 ) + 1 ! add a northern cell
#else
      kpj = ceiling(zbl / ze1 )    ! Global Domain size + ghost cells
#endif

IF (ln_obstacle ) THEN
  IF(nn_AM98<4) WRITE(numout,*) 'nn_AM98 must be >= 4'
  ze0  =  rn_dx / 4._wp                                   ! [m] gridspacing BIGGER cells
  zgl0 =  rn_domsiz + (4._wp + 2._wp * REAL(nn_gc, wp) ) * ze0       ! [m] length of the square with ghostcells
  zbl0 = zgl0 * (  COS( rn_theta * rad ) + SIN( rn_theta * rad ) )   ! length side bigger domain [m]
  !
  kpi = ceiling(zbl0 / ze0 )
  kpj = ceiling(zbl0 / ze0 )
  IF(lwp) WRITE(numout,*) 'before kpi=',kpi,' kpj',kpj
  kpi = kpi * nn_AM98 / 4; kpj = kpj * nn_AM98 / 4
  IF(lwp) WRITE(numout,*) 'after  kpi=',kpi,' kpj',kpj
ENDIF
      !
      IF( rn_modified_grav /= 0._wp) grav = rn_modified_grav   ! update the gravity
      !
      kpk = jpkglo
      !                             ! Set the lateral boundary condition of the global domain
      kperio = 0                    ! AM98 configuration : closed domain
      !
# if defined key_bvp
      r1_abp = 1._wp / rn_abp
#endif
      !                             ! control print
      IF(lwp) THEN
         WRITE(numout,*) '   '
         WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namusr_def : AM98 case'
         WRITE(numout,*) '                                   domain size       rn_domsiz   = ', rn_domsiz, 'm'
         WRITE(numout,*) '                                   gridspacing           rn_dx   = ', rn_dx, 'm'
         WRITE(numout,*) '      inverse resolution & implied domain size         nn_AM98   = ', nn_AM98
         WRITE(numout,*) '                           implied gridspacing           rn_dx   = ', rn_dx, 'm'
         WRITE(numout,*) '                          number of ghostcells           nn_gc   = ', nn_gc
         WRITE(numout,*) '   '
         WRITE(numout,*) '                           minimum porosity             rn_abp   = ', rn_abp
         WRITE(numout,*) '                           friction parameter           rn_fsp   = ', rn_fsp, '1/s'
         WRITE(numout,*) '                           size of the penz             rn_cnp   = ', rn_cnp
         WRITE(numout,*) '                           smoothing parameter          nn_smo   = ', nn_smo
#if defined key_bvp
         WRITE(numout,*) '                           friction type                nn_rfr   = ', nn_rfr
         WRITE(numout,*) '                           bvp friction type            nn_fsp   = ', nn_fsp
#endif
         WRITE(numout,*) '   '
         WRITE(numout,*) '                              wind stress               rn_tau   = ', rn_tau, 'N/m2'
         WRITE(numout,*) '                   linear friction chosen               rn_rfr   = ', rn_rfr, '1/s'
         WRITE(numout,*) '                   1°/4 staircase shape             ln_obstacle  = ', ln_obstacle
         WRITE(numout,*) '                   New viscous boundary cond.  ln_dynldf_lap_PM  = ', ln_dynldf_lap_PM
         WRITE(numout,*) '                   type of stress tensor      nn_dynldf_lap_typ  = ', nn_dynldf_lap_typ
         WRITE(numout,*) '                   Alternative Direction hdiv        ln_hdiv_AD  = ', ln_hdiv_AD
         WRITE(numout,*) '                   Layer with diffusion                 ln_hldc  = ', ln_hldc
         WRITE(numout,*) '                   diffusion parameter                  rn_hldc  = ', rn_hldc, '[m2/s]'
         WRITE(numout,*) '                    rotation angle chosen              rn_theta  = ', rn_theta, 'deg'
         WRITE(numout,*) '                    modified gravity           rn_modified_grav  = ', rn_modified_grav, 'm2/s'
         WRITE(numout,*) '      number of model levels                              jpkglo = ', kpk
         WRITE(numout,*) '   '
         WRITE(numout,*) '   Lateral b.c. of the global domain set to closed        jperio = ', kperio
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
