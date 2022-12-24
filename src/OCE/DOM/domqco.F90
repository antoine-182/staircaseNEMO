MODULE domqco
   !!======================================================================
   !!                       ***  MODULE domqco   ***
   !! Ocean :
   !!======================================================================
   !! History :  2.0  !  2006-06  (B. Levier, L. Marie)  original code
   !!            3.1  !  2009-02  (G. Madec, M. Leclair, R. Benshila)  pure z* coordinate
   !!            3.3  !  2011-10  (M. Leclair) totally rewrote domvvl: vvl option includes z_star and z_tilde coordinates
   !!            3.6  !  2014-11  (P. Mathiot) add ice shelf capability
   !!            4.1  !  2019-08  (A. Coward, D. Storkey) rename dom_vvl_sf_swp -> dom_vvl_sf_update for new timestepping
   !!            4.x  !  2020-02  (G. Madec, S. Techene) pure z* (quasi-eulerian) coordinate
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_qe_init   : define initial vertical scale factors, depths and column thickness
   !!   dom_qe_r3c    : Compute ssh/h_0 ratioat t-, u-, v-, and optionally f-points
   !!       qe_rst_read : read/write restart file
   !!   dom_qe_ctl    : Check the vvl options
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE phycst         ! physical constant
   USE dom_oce        ! ocean space and time domain
   USE dynadv  , ONLY : ln_dynadv_vec
   USE isf_oce        ! iceshelf cavities
   USE sbc_oce        ! ocean surface boundary condition
   USE wet_dry        ! wetting and drying
   USE usrdef_istate  ! user defined initial state (wad only)
   USE restart        ! ocean restart
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! distributed memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC  dom_qco_init       ! called by domain.F90
   PUBLIC  dom_qco_zgr        ! called by isfcpl.F90
   PUBLIC  dom_qco_r3c        ! called by steplf.F90

   !                                                      !!* Namelist nam_vvl
   LOGICAL , PUBLIC :: ln_vvl_zstar           = .FALSE.    ! zstar  vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_ztilde          = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_layer           = .FALSE.    ! level  vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_ztilde_as_zstar = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_zstar_at_eqtor  = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_kepe            = .FALSE.    ! kinetic/potential energy transfer
   !                                                       ! conservation: not used yet
   REAL(wp)         :: rn_ahe3                             ! thickness diffusion coefficient
   REAL(wp)         :: rn_rst_e3t                          ! ztilde to zstar restoration timescale [days]
   REAL(wp)         :: rn_lf_cutoff                        ! cutoff frequency for low-pass filter  [days]
   REAL(wp)         :: rn_zdef_max                         ! maximum fractional e3t deformation
   LOGICAL , PUBLIC :: ln_vvl_dbg = .FALSE.                ! debug control prints

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: un_td, vn_td                ! thickness diffusion transport

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: domvvl.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_qco_init( Kbb, Kmm, Kaa )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_qco_init  ***
      !!
      !! ** Purpose :  Initialization of all ssh. to h._0 ratio
      !!
      !! ** Method  :  - use restart file and/or initialize
      !!               - compute ssh. to h._0 ratio
      !!
      !! ** Action  : - r3(t/u/v)_b
      !!              - r3(t/u/v/f)_n
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kbb, Kmm, Kaa
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dom_qco_init : Variable volume activated'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
      !
      CALL dom_qco_ctl     ! choose vertical coordinate (z_star, z_tilde or layer)
      !
      !                    ! Read or initialize e3t_(b/n), tilde_e3t_(b/n) and hdiv_lf
      CALL qe_rst_read( nit000, Kbb, Kmm )
      !
      CALL dom_qco_zgr(Kbb, Kmm, Kaa) ! interpolation scale factor, depth and water column
      !
      ! IF(lwxios) THEN   ! define variables in restart file when writing with XIOS
      !    CALL iom_set_rstw_var_active('e3t_b')
      !    CALL iom_set_rstw_var_active('e3t_n')
      ! ENDIF
      !
   END SUBROUTINE dom_qco_init


   SUBROUTINE dom_qco_zgr(Kbb, Kmm, Kaa)
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dom_qco_init  ***
      !!
      !! ** Purpose :  Initialization of all ssh. to h._0 ratio
      !!
      !! ** Method  :  - interpolate scale factors
      !!
      !! ** Action  : - r3(t/u/v)_b
      !!              - r3(t/u/v/f)_n
      !!
      !! Reference  : Leclair, M., and G. Madec, 2011, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kbb, Kmm, Kaa
      !!----------------------------------------------------------------------
      !
      !                    !== Set of all other vertical scale factors  ==!  (now and before)
      !                                ! Horizontal interpolation of e3t
      CALL dom_qco_r3c( ssh(:,:,Kbb), r3t(:,:,Kbb), r3u(:,:,Kbb), r3v(:,:,Kbb) )
      CALL dom_qco_r3c( ssh(:,:,Kmm), r3t(:,:,Kmm), r3u(:,:,Kmm), r3v(:,:,Kmm), r3f(:,:) )
      !
   END SUBROUTINE dom_qco_zgr


   SUBROUTINE dom_qco_r3c( pssh, pr3t, pr3u, pr3v, pr3f )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE r3c  ***
      !!
      !! ** Purpose :   compute the filtered ratio ssh/h_0 at t-,u-,v-,f-points
      !!
      !! ** Method  : - compute the ssh at u- and v-points (f-point optional)
      !!                   Vector Form : surface weighted averaging
      !!                   Flux   Form : simple           averaging
      !!              - compute the ratio ssh/h_0 at t-,u-,v-pts, (f-pt optional)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)          , INTENT(in   )  ::   pssh               ! sea surface height   [m]
      REAL(wp), DIMENSION(:,:)          , INTENT(  out)  ::   pr3t, pr3u, pr3v   ! ssh/h0 ratio at t-, u-, v-,points  [-]
      REAL(wp), DIMENSION(:,:), OPTIONAL, INTENT(  out)  ::   pr3f               ! ssh/h0 ratio at f-point   [-]
      !
      INTEGER ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      !
      pr3t(:,:) = pssh(:,:) * r1_ht_0(:,:)   !==  ratio at t-point  ==!
      !
      !
      !                                      !==  ratio at u-,v-point  ==!
      !
      IF( ln_dynadv_vec ) THEN                     !- Vector Form   (thickness weighted averaging)
         DO_2D_00_00
            pr3u(ji,jj) = 0.5_wp * (  e1e2t(ji  ,jj) * pssh(ji  ,jj)  &
               &                    + e1e2t(ji+1,jj) * pssh(ji+1,jj)  ) * r1_hu_0(ji,jj) * r1_e1e2u(ji,jj)
            pr3v(ji,jj) = 0.5_wp * (  e1e2t(ji,jj  ) * pssh(ji,jj  )  &
               &                    + e1e2t(ji,jj+1) * pssh(ji,jj+1)  ) * r1_hv_0(ji,jj) * r1_e1e2v(ji,jj)
         END_2D
      ELSE                                         !- Flux Form   (simple averaging)
         DO_2D_00_00
            pr3u(ji,jj) = 0.5_wp * (  pssh(ji  ,jj) + pssh(ji+1,jj)  ) * r1_hu_0(ji,jj)
            pr3v(ji,jj) = 0.5_wp * (  pssh(ji,jj  ) + pssh(ji,jj+1)  ) * r1_hv_0(ji,jj)
         END_2D
      ENDIF
      !
      IF( .NOT.PRESENT( pr3f ) ) THEN              !- lbc on ratio at u-, v-points only
         CALL lbc_lnk_multi( 'dom_qco_r3c', pr3u, 'U', 1._wp, pr3v, 'V', 1._wp )
         !
         !
      ELSE                                   !==  ratio at f-point  ==!
         !
         IF( ln_dynadv_vec )   THEN                !- Vector Form   (thickness weighted averaging)
            DO_2D_10_10                               ! start from 1 since lbc_lnk('F') doesn't update the 1st row/line
               pr3f(ji,jj) = 0.25_wp * (  e1e2t(ji  ,jj  ) * pssh(ji  ,jj  )  &
                  &                     + e1e2t(ji+1,jj  ) * pssh(ji+1,jj  )  &
                  &                     + e1e2t(ji  ,jj+1) * pssh(ji  ,jj+1)  &
                  &                     + e1e2t(ji+1,jj+1) * pssh(ji+1,jj+1)  ) * r1_hf_0(ji,jj) * r1_e1e2f(ji,jj)
            END_2D
         ELSE                                      !- Flux Form   (simple averaging)
            DO_2D_10_10                               ! start from 1 since lbc_lnk('F') doesn't update the 1st row/line
               pr3f(ji,jj) = 0.25_wp * (  pssh(ji  ,jj  ) + pssh(ji+1,jj  )  &
                  &                     + pssh(ji  ,jj+1) + pssh(ji+1,jj+1)  ) * r1_hf_0(ji,jj)
            END_2D
         ENDIF
         !                                                 ! lbc on ratio at u-,v-,f-points
         CALL lbc_lnk_multi( 'dom_qco_r3c', pr3u, 'U', 1._wp, pr3v, 'V', 1._wp, pr3f, 'F', 1._wp )
         !
      ENDIF
      !
   END SUBROUTINE dom_qco_r3c


   SUBROUTINE qe_rst_read( kt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE qe_rst_read  ***
      !!
      !! ** Purpose :   Read ssh in restart file
      !!
      !! ** Method  :   use of IOM library
      !!                if the restart does not contain ssh,
      !!                it is set to the _0 values.
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt        ! ocean time-step
      INTEGER         , INTENT(in) ::   Kbb, Kmm  ! ocean time level indices
      !
      INTEGER ::   ji, jj, jk
      INTEGER ::   id1, id2     ! local integers
      !!----------------------------------------------------------------------
      !
         IF( ln_rstart ) THEN                   !* Read the restart file
            CALL rst_read_open                  !  open the restart file if necessary
            !
            id1 = iom_varid( numror, 'sshb', ldstop = .FALSE. )
            id2 = iom_varid( numror, 'sshn', ldstop = .FALSE. )
            !
            !                             ! --------- !
            !                             ! all cases !
            !                             ! --------- !
            !
            IF( MIN( id1, id2 ) > 0 ) THEN       ! all required arrays exist
               CALL iom_get( numror, jpdom_autoglo, 'sshb'   , ssh(:,:,Kbb), ldxios = lrxios    )
               CALL iom_get( numror, jpdom_autoglo, 'sshn'   , ssh(:,:,Kmm), ldxios = lrxios    )
               ! needed to restart if land processor not computed
               IF(lwp) write(numout,*) 'qe_rst_read : ssh(:,:,Kbb) and ssh(:,:,Kmm) found in restart files'
               WHERE ( ssmask(:,:) == 0.0_wp )   !!gm/st ==> sm should not be necessary on ssh when it was required on e3
                  ssh(:,:,Kmm) = 0._wp
                  ssh(:,:,Kbb) = 0._wp
               END WHERE
               IF( l_1st_euler ) THEN
                  ssh(:,:,Kbb) = ssh(:,:,Kmm)
               ENDIF
            ELSE IF( id1 > 0 ) THEN
               IF(lwp) write(numout,*) 'qe_rst_read WARNING : ssh(:,:,Kmm) not found in restart files'
               IF(lwp) write(numout,*) 'sshn set equal to sshb.'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               CALL iom_get( numror, jpdom_autoglo, 'sshb', ssh(:,:,Kbb), ldxios = lrxios )
               ssh(:,:,Kmm) = ssh(:,:,Kbb)
               l_1st_euler = .TRUE.
            ELSE IF( id2 > 0 ) THEN
               IF(lwp) write(numout,*) 'qe_rst_read WARNING : ssh(:,:,Kbb) not found in restart files'
               IF(lwp) write(numout,*) 'sshb set equal to sshn.'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               CALL iom_get( numror, jpdom_autoglo, 'sshn', ssh(:,:,Kmm), ldxios = lrxios )
               ssh(:,:,Kbb) = ssh(:,:,Kmm)
               l_1st_euler = .TRUE.
            ELSE
               IF(lwp) write(numout,*) 'qe_rst_read WARNING : ssh(:,:,Kmm) not found in restart file'
               IF(lwp) write(numout,*) 'ssh_b and ssh_n set to zero'
               IF(lwp) write(numout,*) 'neuler is forced to 0'
               ssh(:,:,:) = 0._wp
               l_1st_euler = .TRUE.
            ENDIF
            !
         ELSE                                   !* Initialize at "rest"
            !
            IF( ll_wd ) THEN   ! MJB ll_wd edits start here - these are essential
               !
               IF( cn_cfg == 'wad' ) THEN            ! Wetting and drying test case
                  CALL usr_def_istate( gdept(:,:,:,Kbb), tmask, ts(:,:,:,:,Kbb), uu(:,:,:,Kbb), vv(:,:,:,Kbb), ssh(:,:,Kbb)  )
                  ts (:,:,:,:,Kmm) = ts (:,:,:,:,Kbb)       ! set now values from to before ones
                  ssh(:,:    ,Kmm) = ssh(:,:    ,Kbb)
                  uu (:,:,:  ,Kmm) = uu (:,:,:  ,Kbb)
                  vv (:,:,:  ,Kmm) = vv (:,:,:  ,Kbb)
               ELSE                                  ! if not test case
                  ssh(:,:,Kmm) = -ssh_ref
                  ssh(:,:,Kbb) = -ssh_ref
                  !
                  DO_2D_11_11
                     IF( ht_0(ji,jj)-ssh_ref <  rn_wdmin1 ) THEN ! if total depth is less than min depth
                        ssh(ji,jj,Kbb) = rn_wdmin1 - (ht_0(ji,jj) )
                        ssh(ji,jj,Kmm) = rn_wdmin1 - (ht_0(ji,jj) )
                     ENDIF
                  END_2D
               ENDIF

               DO ji = 1, jpi
                  DO jj = 1, jpj
                     IF ( ht_0(ji,jj) .LE. 0.0 .AND. NINT( ssmask(ji,jj) ) .EQ. 1) THEN
                       CALL ctl_stop( 'qe_rst_read: ht_0 must be positive at potentially wet points' )
                     ENDIF
                  END DO
               END DO
               !
            ELSE
               !
               ! Just to read set ssh in fact, called latter once vertical grid
               ! is set up:
!               CALL usr_def_istate( gdept_0, tmask, ts(:,:,:,:,Kbb), uu(:,:,:,Kbb), vv(:,:,:,Kbb), ssh(:,:,Kbb)  )
!               !
                ssh(:,:,:) = 0._wp
               !
            ENDIF           ! end of ll_wd edits
            !
         ENDIF
      !
   END SUBROUTINE qe_rst_read


   SUBROUTINE dom_qco_ctl
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_qco_ctl  ***
      !!
      !! ** Purpose :   Control the consistency between namelist options
      !!                for vertical coordinate
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ios
      !!
      NAMELIST/nam_vvl/ ln_vvl_zstar, ln_vvl_ztilde, ln_vvl_layer, ln_vvl_ztilde_as_zstar, &
         &              ln_vvl_zstar_at_eqtor      , rn_ahe3     , rn_rst_e3t            , &
         &              rn_lf_cutoff               , rn_zdef_max , ln_vvl_dbg                ! not yet implemented: ln_vvl_kepe
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_ref, nam_vvl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nam_vvl in reference namelist' )
      READ  ( numnam_cfg, nam_vvl, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'nam_vvl in configuration namelist' )
      IF(lwm) WRITE ( numond, nam_vvl )
      !
      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_qco_ctl : choice/control of the variable vertical coordinate'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist nam_vvl : chose a vertical coordinate'
         WRITE(numout,*) '      zstar                      ln_vvl_zstar   = ', ln_vvl_zstar
         WRITE(numout,*) '      ztilde                     ln_vvl_ztilde  = ', ln_vvl_ztilde
         WRITE(numout,*) '      layer                      ln_vvl_layer   = ', ln_vvl_layer
         WRITE(numout,*) '      ztilde as zstar   ln_vvl_ztilde_as_zstar  = ', ln_vvl_ztilde_as_zstar
         WRITE(numout,*) '      ztilde near the equator    ln_vvl_zstar_at_eqtor  = ', ln_vvl_zstar_at_eqtor
         WRITE(numout,*) '      !'
         WRITE(numout,*) '      thickness diffusion coefficient                      rn_ahe3      = ', rn_ahe3
         WRITE(numout,*) '      maximum e3t deformation fractional change            rn_zdef_max  = ', rn_zdef_max
         IF( ln_vvl_ztilde_as_zstar ) THEN
            WRITE(numout,*) '      ztilde running in zstar emulation mode (ln_vvl_ztilde_as_zstar=T) '
            WRITE(numout,*) '         ignoring namelist timescale parameters and using:'
            WRITE(numout,*) '            hard-wired : z-tilde to zstar restoration timescale (days)'
            WRITE(numout,*) '                         rn_rst_e3t     = 0.e0'
            WRITE(numout,*) '            hard-wired : z-tilde cutoff frequency of low-pass filter (days)'
            WRITE(numout,*) '                         rn_lf_cutoff   = 1.0/rn_Dt'
         ELSE
            WRITE(numout,*) '      z-tilde to zstar restoration timescale (days)        rn_rst_e3t   = ', rn_rst_e3t
            WRITE(numout,*) '      z-tilde cutoff frequency of low-pass filter (days)   rn_lf_cutoff = ', rn_lf_cutoff
         ENDIF
         WRITE(numout,*) '         debug prints flag                                 ln_vvl_dbg   = ', ln_vvl_dbg
      ENDIF
      !
      ioptio = 0                      ! Parameter control
      IF( ln_vvl_ztilde_as_zstar )   ln_vvl_ztilde = .true.
      IF( ln_vvl_zstar           )   ioptio = ioptio + 1
      IF( ln_vvl_ztilde          )   ioptio = ioptio + 1
      IF( ln_vvl_layer           )   ioptio = ioptio + 1
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE vertical coordinate in namelist nam_vvl' )
      !
      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         IF( ln_vvl_zstar           ) WRITE(numout,*) '      ==>>>   zstar vertical coordinate is used'
         IF( ln_vvl_ztilde          ) WRITE(numout,*) '      ==>>>   ztilde vertical coordinate is used'
         IF( ln_vvl_layer           ) WRITE(numout,*) '      ==>>>   layer vertical coordinate is used'
         IF( ln_vvl_ztilde_as_zstar ) WRITE(numout,*) '      ==>>>   to emulate a zstar coordinate'
      ENDIF
      !
#if defined key_agrif
      IF( (.NOT.Agrif_Root()).AND.(.NOT.ln_vvl_zstar) )   CALL ctl_stop( 'AGRIF is implemented with zstar coordinate only' )
#endif
      !
   END SUBROUTINE dom_qco_ctl

   !!======================================================================
END MODULE domqco
