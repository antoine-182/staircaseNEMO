MODULE stpctl
   !!======================================================================
   !!                       ***  MODULE  stpctl  ***
   !! Ocean run control :  gross check of the ocean time stepping
   !!                      version for STATION_ASF test-case
   !!======================================================================
   !! History :  OPA  ! 1991-03  (G. Madec) Original code
   !!            6.0  ! 1992-06  (M. Imbard)
   !!            8.0  ! 1997-06  (A.M. Treguier)
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            2.0  ! 2009-07  (G. Madec)  Add statistic for time-spliting
   !!            3.7  ! 2016-09  (G. Madec)  Remove solver
   !!            4.0  ! 2017-04  (G. Madec)  regroup global communications
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp_ctl      : Control the run
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain variables
   USE sbc_oce         ! surface fluxes and stuff
   USE diawri          ! Standard run outputs       (dia_wri_state routine)
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing

   USE netcdf          ! NetCDF library
   IMPLICIT NONE
   PRIVATE

   PUBLIC stp_ctl           ! routine called by step.F90

   INTEGER  ::   idrun, idtime, idtau, idqns, idemp, istatus
   LOGICAL  ::   lsomeoce
   !!----------------------------------------------------------------------
   !! NEMO/SAS 4.0 , NEMO Consortium (2018)
   !! $Id: stpctl.F90 10603 2019-01-29 11:18:09Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp_ctl( kt, Kbb, Kmm, kindic )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl  ***
      !!
      !! ** Purpose :   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Print it each 50 time steps
      !!              - Stop the run IF problem encountered by setting indic=-3
      !!
      !! ** Actions :   "time.step" file = last ocean time-step
      !!                "run.stat"  file = run statistics
      !!                nstop indicator sheared among all local domain (lk_mpp=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER, INTENT(in   ) ::   Kbb, Kmm      ! ocean time level index
      INTEGER, INTENT(inout) ::   kindic   ! error indicator
      !!
      REAL(wp), DIMENSION(3) ::   zmax
      LOGICAL                ::   ll_wrtstp, ll_colruns, ll_wrtruns
      CHARACTER(len=20) :: clname
      !!----------------------------------------------------------------------
      !
      ll_wrtstp  = ( MOD( kt, sn_cfctl%ptimincr ) == 0 ) .OR. ( kt == nitend )
      ll_colruns = ll_wrtstp .AND. ( sn_cfctl%l_runstat )
      ll_wrtruns = ll_colruns .AND. lwm
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'stp_ctl : time-stepping control'
         WRITE(numout,*) '~~~~~~~'
         !                                ! open time.step file
         IF( lwm ) CALL ctl_opn( numstp, 'time.step', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
         !                                ! open run.stat file(s) at start whatever
         !                                ! the value of sn_cfctl%ptimincr
         IF( lwm .AND. ( sn_cfctl%l_runstat ) ) THEN
            CALL ctl_opn( numrun, 'run.stat', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
            clname = 'run.stat.nc'
            IF( .NOT. Agrif_Root() )   clname = TRIM(Agrif_CFixed())//"_"//TRIM(clname)
            istatus = NF90_CREATE( TRIM(clname), NF90_CLOBBER, idrun )
            istatus = NF90_DEF_DIM( idrun, 'time', NF90_UNLIMITED, idtime )
            istatus = NF90_DEF_VAR( idrun, 'tau_max', NF90_DOUBLE, (/ idtime /), idtau )
            istatus = NF90_DEF_VAR( idrun, 'qns_max', NF90_DOUBLE, (/ idtime /), idqns   )
            istatus = NF90_DEF_VAR( idrun, 'emp_max', NF90_DOUBLE, (/ idtime /), idemp   )
            istatus = NF90_ENDDEF(idrun)
         ENDIF
      ENDIF
      IF( kt == nit000 )   lsomeoce = COUNT( ssmask(:,:) == 1._wp ) > 0
      !
      IF(lwm .AND. ll_wrtstp) THEN        !==  current time step  ==!   ("time.step" file)
         WRITE ( numstp, '(1x, i8)' )   kt
         REWIND( numstp )
      ENDIF
      !
      !                                   !==  test of extrema  ==!
      zmax(1) = MAXVAL(     taum(:,:)   , mask = tmask(:,:,1) == 1._wp )                                         ! max wind stress module
      zmax(2) = MAXVAL( ABS( qns(:,:) ) , mask = tmask(:,:,1) == 1._wp )                                         ! max non-solar heat flux
      zmax(3) = MAXVAL( ABS( emp(:,:) ) , mask = tmask(:,:,1) == 1._wp )                                         ! max E-P
      !
      IF( ll_colruns ) THEN
         CALL mpp_max( "stpctl", zmax )          ! max over the global domain
         nstop = NINT( zmax(3) )                 ! nstop indicator sheared among all local domains
      ENDIF
      !                                   !==  run statistics  ==!   ("run.stat" files)
      IF( ll_wrtruns ) THEN
         WRITE(numrun,9500) kt, zmax(1), zmax(2), zmax(3)
         istatus = NF90_PUT_VAR( idrun, idtau, (/ zmax(1)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun, idqns, (/ zmax(2)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun, idemp, (/ zmax(3)/), (/kt/), (/1/) )
         IF( MOD( kt , 100 ) == 0 ) istatus = NF90_SYNC(idrun)
         IF( kt == nitend         ) istatus = NF90_CLOSE(idrun)
      END IF
      !                                   !==  error handling  ==!
      IF( ( sn_cfctl%l_glochk .OR. lsomeoce ) .AND. (   &  ! domain contains some ocean points, check for sensible ranges
         &  zmax(1) >    5._wp .OR.   &             ! too large wind stress ( > 5 N/m^2 )
         &  zmax(2) > 2000._wp .OR.   &             ! too large non-solar heat flux ( > 2000 W/m^2)
         &  zmax(3) > 1.E-3_wp .OR.   &             ! too large net freshwater flux ( kg/m^2/s)
         &  ISNAN( zmax(1) + zmax(2) + zmax(3) ) ) ) THEN   ! NaN encounter in the tests

         !! We are 1D so no need to find a spatial location of the rogue point.

         WRITE(ctmp1,*) ' stp_ctl: |tau_mod| > 5 N/m2  or  |qns| > 2000 W/m2  or |emp| > 1.E-3 or  NaN encounter in the tests'
         WRITE(ctmp2,9500) kt,   zmax(1), zmax(2), zmax(3)
         WRITE(ctmp6,*) '      ===> output of last computed fields in output.abort.nc file'

         CALL dia_wri_state( Kmm, 'output.abort' )     ! create an output.abort file

         IF( .NOT. sn_cfctl%l_glochk ) THEN
            WRITE(ctmp8,*) 'E R R O R message from sub-domain: ', narea
            CALL ctl_stop( 'STOP', ctmp1, ' ', ctmp2, ' ', ctmp6, ' ' )
         ELSE
            CALL ctl_stop( ctmp1, ' ', ctmp2, ' ', ctmp6, ' ' )
         ENDIF

         kindic = -3
         !
      ENDIF
      !
9500  FORMAT(' it :', i8, '    tau_max: ', D23.16, ' |qns|_max: ', D23.16,' |emp|_max: ', D23.16)
      !
   END SUBROUTINE stp_ctl

   !!======================================================================
END MODULE stpctl
