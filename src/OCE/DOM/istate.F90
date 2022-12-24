MODULE istate
   !!======================================================================
   !!                     ***  MODULE  istate  ***
   !! Ocean state   :  initial state setting
   !!=====================================================================
   !! History :  OPA  !  1989-12  (P. Andrich)  Original code
   !!            5.0  !  1991-11  (G. Madec)  rewritting
   !!            6.0  !  1996-01  (G. Madec)  terrain following coordinates
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_eel
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_uvg
   !!   NEMO     1.0  !  2003-08  (G. Madec, C. Talandier)  F90: Free form, modules + EEL R5
   !!             -   !  2004-05  (A. Koch-Larrouy)  istate_gyre 
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!            3.3  !  2010-10  (C. Ethe) merge TRC-TRA
   !!            3.4  !  2011-04  (G. Madec) Merge of dtatem and dtasal & suppression of tb,tn/sb,sn 
   !!            3.7  !  2016-04  (S. Flavoni) introduce user defined initial state 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   istate_init   : initial state setting
   !!   istate_uvg    : initial velocity in geostropic balance
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers 
   USE dom_oce        ! ocean space and time domain 
   USE daymod         ! calendar
   USE divhor         ! horizontal divergence            (div_hor routine)
   USE dtatsd         ! data temperature and salinity   (dta_tsd routine)
   USE dtauvd         ! data: U & V current             (dta_uvd routine)
   USE domvvl          ! varying vertical mesh
   USE wet_dry         ! wetting and drying (needed for wad_istate)
   USE usrdef_istate   ! User defined initial state
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lib_mpp         ! MPP library
   USE restart         ! restart

   IMPLICIT NONE
   PRIVATE

   PUBLIC   istate_init   ! routine called by step.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: istate.F90 12489 2020-02-28 15:55:11Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE istate_init( Kbb, Kmm, Kaa )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_init  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracer fields.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::  Kbb, Kmm, Kaa   ! ocean time level indices
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zgdept     ! 3D table  !!st patch to use gdept subtitute
!!gm see comment further down
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   zuvd    ! U & V data workspace
!!gm end
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'istate_init : Initialization of the dynamics and tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

!!gm  Why not include in the first call of dta_tsd ?  
!!gm  probably associated with the use of internal damping...
                     CALL dta_tsd_init        ! Initialisation of T & S input data
!!gm to be moved in usrdef of C1D case
!      IF( lk_c1d )   CALL dta_uvd_init        ! Initialization of U & V input data
!!gm

      rhd  (:,:,:  ) = 0._wp   ;   rhop (:,:,:  ) = 0._wp      ! set one for all to 0 at level jpk
      rn2b (:,:,:  ) = 0._wp   ;   rn2  (:,:,:  ) = 0._wp      ! set one for all to 0 at levels 1 and jpk
      ts  (:,:,:,:,Kaa) = 0._wp                                   ! set one for all to 0 at level jpk
      rab_b(:,:,:,:) = 0._wp   ;   rab_n(:,:,:,:) = 0._wp      ! set one for all to 0 at level jpk
#if defined key_agrif
      uu   (:,:,:  ,Kaa) = 0._wp   ! used in agrif_oce_sponge at initialization
      vv   (:,:,:  ,Kaa) = 0._wp   ! used in agrif_oce_sponge at initialization    
#endif

      IF( ln_rstart ) THEN                    ! Restart from a file
         !                                    ! -------------------
         CALL rst_read( Kbb, Kmm )            ! Read the restart file
         CALL day_init                        ! model calendar (using both namelist and restart infos)
         !
      ELSE                                    ! Start from rest
         !                                    ! ---------------
         numror = 0                           ! define numror = 0 -> no restart file to read
         l_1st_euler = .true.                 ! Set time-step indicator at nit000 (euler forward)
         CALL day_init                        ! model calendar (using both namelist and restart infos)
         !                                    ! Initialization of ocean to zero
         !
         IF( ln_tsd_init ) THEN               
            CALL dta_tsd( nit000, ts(:,:,:,:,Kbb) )       ! read 3D T and S data at nit000
            !
            ssh(:,:,Kbb)   = 0._wp               ! set the ocean at rest
            IF( ll_wd ) THEN
               ssh(:,:,Kbb) =  -ssh_ref  ! Added in 30 here for bathy that adds 30 as Iterative test CEOD 
               !
               ! Apply minimum wetdepth criterion
               !
               DO_2D_11_11
                  IF( ht_0(ji,jj) + ssh(ji,jj,Kbb)  < rn_wdmin1 ) THEN
                     ssh(ji,jj,Kbb) = tmask(ji,jj,1)*( rn_wdmin1 - (ht_0(ji,jj)) )
                  ENDIF
               END_2D
            ENDIF 
            uu  (:,:,:,Kbb) = 0._wp
            vv  (:,:,:,Kbb) = 0._wp  
            !
         ELSE                                 ! user defined initial T and S
            DO jk = 1, jpk
               zgdept(:,:,jk) = gdept(:,:,jk,Kbb)
            END DO
            CALL usr_def_istate( zgdept, tmask, ts(:,:,:,:,Kbb), uu(:,:,:,Kbb), vv(:,:,:,Kbb), ssh(:,:,Kbb)  )         
         ENDIF
         ts  (:,:,:,:,Kmm) = ts (:,:,:,:,Kbb)       ! set now values from to before ones
         ssh (:,:,Kmm)     = ssh(:,:,Kbb)   
         uu   (:,:,:,Kmm)   = uu  (:,:,:,Kbb)
         vv   (:,:,:,Kmm)   = vv  (:,:,:,Kbb)
         hdiv(:,:,jpk) = 0._wp               ! bottom divergence set one for 0 to zero at jpk level
         CALL div_hor( 0, Kbb, Kmm )         ! compute interior hdiv value  
!!gm                                    hdiv(:,:,:) = 0._wp

!!gm POTENTIAL BUG :
!!gm  ISSUE :  if ssh(:,:,Kbb) /= 0  then, in non linear free surface, the e3._n, e3._b should be recomputed
!!             as well as gdept_ and gdepw_....   !!!!! 
!!      ===>>>>   probably a call to domvvl initialisation here....


         !
!!gm to be moved in usrdef of C1D case
!         IF ( ln_uvd_init .AND. lk_c1d ) THEN ! read 3D U and V data at nit000
!            ALLOCATE( zuvd(jpi,jpj,jpk,2) )
!            CALL dta_uvd( nit000, zuvd )
!            uu(:,:,:,Kbb) = zuvd(:,:,:,1) ;  uu(:,:,:,Kmm) = uu(:,:,:,Kbb)
!            vv(:,:,:,Kbb) = zuvd(:,:,:,2) ;  vv(:,:,:,Kmm) = vv(:,:,:,Kbb)
!            DEALLOCATE( zuvd )
!         ENDIF
         !
!!gm This is to be changed !!!!
!         ! - ML - ssh(:,:,Kmm) could be modified by istate_eel, so that initialization of e3t(:,:,:,Kbb) is done here
!         IF( .NOT.ln_linssh ) THEN
!            DO jk = 1, jpk
!               e3t(:,:,jk,Kbb) = e3t(:,:,jk,Kmm)
!            END DO
!         ENDIF
!!gm 
         ! 
      ENDIF 
      ! 
      ! Initialize "now" and "before" barotropic velocities:
      ! Do it whatever the free surface method, these arrays being eventually used
      !
      uu_b(:,:,Kmm) = 0._wp   ;   vv_b(:,:,Kmm) = 0._wp
      uu_b(:,:,Kbb) = 0._wp   ;   vv_b(:,:,Kbb) = 0._wp
      !
!!gm  the use of umsak & vmask is not necessary below as uu(:,:,:,Kmm), vv(:,:,:,Kmm), uu(:,:,:,Kbb), vv(:,:,:,Kbb) are always masked
      DO_3D_11_11( 1, jpkm1 )
         uu_b(ji,jj,Kmm) = uu_b(ji,jj,Kmm) + e3u(ji,jj,jk,Kmm) * uu(ji,jj,jk,Kmm) * umask(ji,jj,jk)
         vv_b(ji,jj,Kmm) = vv_b(ji,jj,Kmm) + e3v(ji,jj,jk,Kmm) * vv(ji,jj,jk,Kmm) * vmask(ji,jj,jk)
         !
         uu_b(ji,jj,Kbb) = uu_b(ji,jj,Kbb) + e3u(ji,jj,jk,Kbb) * uu(ji,jj,jk,Kbb) * umask(ji,jj,jk)
         vv_b(ji,jj,Kbb) = vv_b(ji,jj,Kbb) + e3v(ji,jj,jk,Kbb) * vv(ji,jj,jk,Kbb) * vmask(ji,jj,jk)
      END_3D
      !
      uu_b(:,:,Kmm) = uu_b(:,:,Kmm) * r1_hu(:,:,Kmm)
      vv_b(:,:,Kmm) = vv_b(:,:,Kmm) * r1_hv(:,:,Kmm)
      !
      uu_b(:,:,Kbb) = uu_b(:,:,Kbb) * r1_hu(:,:,Kbb)
      vv_b(:,:,Kbb) = vv_b(:,:,Kbb) * r1_hv(:,:,Kbb)
      !
   END SUBROUTINE istate_init

   !!======================================================================
END MODULE istate
