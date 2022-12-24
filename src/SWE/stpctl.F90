MODULE stpctl
   !!======================================================================
   !!                       ***  MODULE  stpctl  ***
   !! Ocean run control :  gross check of the ocean time stepping
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
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE c1d             ! 1D vertical configuration
   USE diawri          ! Standard run outputs       (dia_wri_state routine)
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing
   USE zdf_oce ,  ONLY : ln_zad_Aimp       ! ocean vertical physics variables
   USE wet_dry,   ONLY : ll_wd, ssh_ref    ! reference depth for negative bathy

   USE netcdf          ! NetCDF library
   IMPLICIT NONE
   PRIVATE

   PUBLIC stp_ctl           ! routine called by step.F90

   INTEGER  ::   idrun, idtime, idssh, idu, ids1, ids2, idt1, idt2, idc1, idw1, istatus
   LOGICAL  ::   lsomeoce
!!stoops
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: stpctl.F90 12983 2020-05-27 15:46:51Z techene $
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
      !!                Problems checked: |ssh| maximum larger than 10 m
      !!                                  |U|   maximum larger than 10 m/s
      !!                                  negative sea surface salinity
      !!
      !! ** Actions :   "time.step" file = last ocean time-step
      !!                "run.stat"  file = run statistics
      !!                nstop indicator sheared among all local domain (lk_mpp=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER, INTENT(in   ) ::   Kbb, Kmm      ! ocean time level index
      INTEGER, INTENT(inout) ::   kindic   ! error indicator
      !!
      INTEGER                ::   ji, jj, jk          ! dummy loop indices
      INTEGER, DIMENSION(2)  ::   ih                  ! min/max loc indices
      INTEGER, DIMENSION(3)  ::   iu, is1, is2        ! min/max loc indices
      REAL(wp)               ::   zzz                 ! local real
      REAL(wp), DIMENSION(3) ::   zmax
      LOGICAL                ::   ll_wrtstp, ll_colruns, ll_wrtruns
      CHARACTER(len=20) :: clname
      ! REAL, DIMENSION(:,:,:) ::z3d
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
            istatus = NF90_DEF_VAR( idrun, 'abs_ssh_max', NF90_DOUBLE, (/ idtime /), idssh )
            istatus = NF90_DEF_VAR( idrun,   'abs_u_max', NF90_DOUBLE, (/ idtime /), idu   )
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
      IF( ll_wd ) THEN
         zmax(1) = MAXVAL(  ABS( ssh(:,:,Kmm) + ssh_ref*tmask(:,:,1) )  )        ! ssh max
      ELSE
         zmax(1) = MINVAL( e3t(:,:,1,Kmm)  )                                         ! ssh min
      ENDIF
!!an   ne pas utiliser
      ! z3d = uu(:,:,:,Kmm)*rpo(:,:,:)
      ! zmax(2) = MAXVAL( ABS( z3d )  )                                  ! velocity max (zonal only)
!!an
      zmax(2) = MAXVAL( ABS( uu(:,:,:,Kmm) )  )                                  ! velocity max (zonal only)
      zmax(3) = REAL( nstop , wp )                                            ! stop indicator
      !
      IF( ll_colruns ) THEN
         CALL mpp_max( "stpctl", zmax )          ! max over the global domain
         nstop = NINT( zmax(3) )                 ! nstop indicator sheared among all local domains
      ENDIF
      !                                   !==  run statistics  ==!   ("run.stat" files)
      IF( ll_wrtruns ) THEN
         WRITE(numrun,9500) kt, zmax(1), zmax(2)
         istatus = NF90_PUT_VAR( idrun, idssh, (/ zmax(1)/), (/kt/), (/1/) )
         istatus = NF90_PUT_VAR( idrun,   idu, (/ zmax(2)/), (/kt/), (/1/) )
         IF( MOD( kt , 100 ) == 0 ) istatus = NF90_SYNC(idrun)
         IF( kt == nitend         ) istatus = NF90_CLOSE(idrun)
      END IF
      !                                   !==  error handling  ==!
      IF( ( sn_cfctl%l_glochk .OR. lsomeoce ) .AND. (   &  ! domain contains some ocean points, check for sensible ranges
         &  zmax(1) <    0._wp .OR.   &                    ! negative sea surface height
         ! &  zmax(2) >   10._wp .OR.   &                    ! too large velocity ( > 10 m/s)
         &  zmax(2) >   10._wp .OR.   &                    ! too large velocity ( > 10 m/s)
         &  ISNAN( zmax(1) + zmax(2) ) ) ) THEN            ! NaN encounter in the tests
         IF( lk_mpp .AND. sn_cfctl%l_glochk ) THEN
            ! have use mpp_max (because sn_cfctl%l_glochk=.T. and distributed)
            CALL mpp_maxloc( 'stpctl', ABS(ssh(:,:,Kmm))        , ssmask(:,:)  , zzz, ih  )
            CALL mpp_maxloc( 'stpctl', ABS(uu(:,:,:,Kmm))          , umask (:,:,:), zzz, iu  )
         ELSE
            ! find local min and max locations
            ih(:)  = MAXLOC( ABS( ssh(:,:,Kmm)   )                              ) + (/ nimpp - 1, njmpp - 1    /)
            iu(:)  = MAXLOC( ABS( uu  (:,:,:,Kmm) )                              ) + (/ nimpp - 1, njmpp - 1, 0 /)
         ENDIF

         WRITE(ctmp1,*) ' stp_ctl: (e3t0) ssh < 0 m  or  |U| > 10 m/s  or  NaN encounter in the tests'
         WRITE(ctmp2,9100) kt,   zmax(1), ih(1) , ih(2)
         WRITE(ctmp3,9200) kt,   zmax(2), iu(1) , iu(2) , iu(3)
         WRITE(ctmp4,*) '      ===> output of last computed fields in output.abort.nc file'

         CALL dia_wri_state( Kmm, 'output.abort' )     ! create an output.abort file

         IF( .NOT. sn_cfctl%l_glochk ) THEN
            WRITE(ctmp8,*) 'E R R O R message from sub-domain: ', narea
            CALL ctl_stop( 'STOP', ctmp1, ' ', ctmp8, ' ', ctmp2, ctmp3, ctmp4 )
         ELSE
            CALL ctl_stop( ctmp1, ' ', ctmp2, ctmp3, ctmp4 )
         ENDIF

         kindic = -3
         !
      ENDIF
      !
9100  FORMAT (' kt=',i8,'   |ssh| min: ',1pg11.4,', at  i j  : ',2i5)
9200  FORMAT (' kt=',i8,'   |U|   max: ',1pg11.4,', at  i j k: ',3i5)
9500  FORMAT(' it :', i8, '    |ssh|_max: ', D23.16, ' |U|_max: ', D23.16)
      !
   END SUBROUTINE stp_ctl

   !!======================================================================
END MODULE stpctl
