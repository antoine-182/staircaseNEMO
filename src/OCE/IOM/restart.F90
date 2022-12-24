MODULE restart
   !!======================================================================
   !!                     ***  MODULE  restart  ***
   !! Ocean restart :  write the ocean restart file
   !!======================================================================
   !! History :  OPA  !  1999-11  (M. Imbard)  Original code
   !!   NEMO     1.0  !  2002-08  (G. Madec)  F90: Free form
   !!            2.0  !  2006-07  (S. Masson)  use IOM for restart
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  modified LF-RA
   !!            - -  !  2010-10  (C. Ethe, G. Madec) TRC-TRA merge (T-S in 4D)
   !!            3.7  !  2014-01  (G. Madec) suppression of curl and hdiv from the restart
   !!             -   !  2014-12  (G. Madec) remove KPP scheme
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   rst_opn    : open the ocean restart file
   !!   rst_write  : write the ocean restart file
   !!   rst_read   : read the ocean restart file
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE sbc_ice         ! only lk_si3 
   USE phycst          ! physical constants
   USE eosbn2          ! equation of state            (eos bn2 routine)
   USE trdmxl_oce      ! ocean active mixed layer tracers trends variables
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module
   USE diu_bulk
   USE lib_mpp         ! distribued memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   rst_opn         ! routine called by step module
   PUBLIC   rst_write       ! routine called by step module
   PUBLIC   rst_read        ! routine called by istate module
   PUBLIC   rst_read_open   ! routine called in rst_read and (possibly) in dom_vvl_init

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: restart.F90 12489 2020-02-28 15:55:11Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE rst_opn( kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rst_opn  ***
      !!                     
      !! ** Purpose : + initialization (should be read in the namelist) of nitrst 
      !!              + open the restart when we are one time step before nitrst
      !!                   - restart header is defined when kt = nitrst-1
      !!                   - restart data  are written when kt = nitrst
      !!              + define lrst_oce to .TRUE. when we need to define or write the restart
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      !!
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step deine as a character
      CHARACTER(LEN=50)   ::   clname   ! ocean output restart file name
      CHARACTER(lc)       ::   clpath   ! full path to ocean output restart file
      CHARACTER(LEN=52)   ::   clpname   ! ocean output restart file name including prefix for AGRIF
      CHARACTER(LEN=256)  ::   clinfo    ! info character
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN   ! default definitions
         lrst_oce = .FALSE.   
         IF( ln_rst_list ) THEN
            nrst_lst = 1
            nitrst = nn_stocklist( nrst_lst )
         ELSE
            nitrst = nitend
         ENDIF
      ENDIF
      
      IF( .NOT. ln_rst_list .AND. nn_stock == -1 )   RETURN   ! we will never do any restart

      ! frequency-based restart dumping (nn_stock)
      IF( .NOT. ln_rst_list .AND. MOD( kt - 1, nn_stock ) == 0 ) THEN   
         ! we use kt - 1 and not kt - nit000 to keep the same periodicity from the beginning of the experiment
         nitrst = kt + nn_stock - 1                  ! define the next value of nitrst for restart writing
         IF( nitrst > nitend )   nitrst = nitend   ! make sure we write a restart at the end of the run
      ENDIF
      ! to get better performances with NetCDF format:
      ! we open and define the ocean restart file one time step before writing the data (-> at nitrst - 1)
      ! except if we write ocean restart files every time step or if an ocean restart file was writen at nitend - 1
      IF( kt == nitrst - 1 .OR. nn_stock == 1 .OR. ( kt == nitend .AND. .NOT. lrst_oce ) ) THEN
         IF( nitrst <= nitend .AND. nitrst > 0 ) THEN 
            ! beware of the format used to write kt (default is i8.8, that should be large enough...)
            IF( nitrst > 999999999 ) THEN   ;   WRITE(clkt, *       ) nitrst
            ELSE                            ;   WRITE(clkt, '(i8.8)') nitrst
            ENDIF
            ! create the file
            clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_"//TRIM(cn_ocerst_out)
            clpath = TRIM(cn_ocerst_outdir)
            IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
            IF(lwp) THEN
               WRITE(numout,*)
               IF(.NOT.lwxios) THEN
                  WRITE(numout,*) '             open ocean restart NetCDF file: ',TRIM(clpath)//TRIM(clname)
                  IF ( snc4set%luse )      WRITE(numout,*) '             opened for NetCDF4 chunking and compression'
                  IF( kt == nitrst - 1 ) THEN   ;   WRITE(numout,*) '             kt = nitrst - 1 = ', kt
                  ELSE                          ;   WRITE(numout,*) '             kt = '             , kt
                  ENDIF
               ENDIF
            ENDIF
            !
            IF(.NOT.lwxios) THEN
               CALL iom_open( TRIM(clpath)//TRIM(clname), numrow, ldwrt = .TRUE. )
            ELSE
#if defined key_iomput
               cwxios_context = "rstw_"//TRIM(ADJUSTL(clkt))
               IF( TRIM(Agrif_CFixed()) == '0' ) THEN
                  clpname = clname
               ELSE
                  clpname = TRIM(Agrif_CFixed())//"_"//clname   
               ENDIF
               CALL iom_init( cwxios_context, TRIM(clpath)//TRIM(clpname), .false. )
               CALL xios_update_calendar(nitrst)
               CALL iom_swap(      cxios_context          )
#else
               clinfo = 'Can not use XIOS in rst_opn'
               CALL ctl_stop(TRIM(clinfo))
#endif
            ENDIF
            lrst_oce = .TRUE.
         ENDIF
      ENDIF
      !
   END SUBROUTINE rst_opn


   SUBROUTINE rst_write( kt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rstwrite  ***
      !!                     
      !! ** Purpose :   Write restart fields in NetCDF format
      !!
      !! ** Method  :   Write in numrow when kt == nitrst in NetCDF
      !!              file, save fields which are necessary for restart
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt         ! ocean time-step
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      !!----------------------------------------------------------------------
                     IF(lwxios) CALL iom_swap(      cwxios_context          )
                     CALL iom_rstput( kt, nitrst, numrow, 'rdt'    , rn_Dt       , ldxios = lwxios)   ! dynamics time step
                     CALL iom_delay_rst( 'WRITE', 'OCE', numrow )   ! save only ocean delayed global communication variables

      IF ( .NOT. ln_diurnal_only ) THEN
                     CALL iom_rstput( kt, nitrst, numrow, 'ub'     , uu(:,:,:       ,Kbb), ldxios = lwxios        )     ! before fields
                     CALL iom_rstput( kt, nitrst, numrow, 'vb'     , vv(:,:,:       ,Kbb), ldxios = lwxios        )
                     CALL iom_rstput( kt, nitrst, numrow, 'tb'     , ts(:,:,:,jp_tem,Kbb), ldxios = lwxios )
                     CALL iom_rstput( kt, nitrst, numrow, 'sb'     , ts(:,:,:,jp_sal,Kbb), ldxios = lwxios )
                     CALL iom_rstput( kt, nitrst, numrow, 'sshb'   ,ssh(:,:         ,Kbb), ldxios = lwxios      )
                     !
                     CALL iom_rstput( kt, nitrst, numrow, 'un'     , uu(:,:,:       ,Kmm), ldxios = lwxios        )     ! now fields
                     CALL iom_rstput( kt, nitrst, numrow, 'vn'     , vv(:,:,:       ,Kmm), ldxios = lwxios        )
                     CALL iom_rstput( kt, nitrst, numrow, 'tn'     , ts(:,:,:,jp_tem,Kmm), ldxios = lwxios )
                     CALL iom_rstput( kt, nitrst, numrow, 'sn'     , ts(:,:,:,jp_sal,Kmm), ldxios = lwxios )
                     CALL iom_rstput( kt, nitrst, numrow, 'sshn'   ,ssh(:,:         ,Kmm), ldxios = lwxios      )
                     CALL iom_rstput( kt, nitrst, numrow, 'rhop'   , rhop, ldxios = lwxios      )
      ENDIF
      
      IF (ln_diurnal) CALL iom_rstput( kt, nitrst, numrow, 'Dsst', x_dsst, ldxios = lwxios )  
      IF(lwxios) CALL iom_swap(      cxios_context          )
      IF( kt == nitrst ) THEN
         IF(.NOT.lwxios) THEN
            CALL iom_close( numrow )     ! close the restart file (only at last time step)
         ELSE
            CALL iom_context_finalize(      cwxios_context          )
         ENDIF
!!gm         IF( .NOT. lk_trdmld )   lrst_oce = .FALSE.
!!gm  not sure what to do here   ===>>>  ask to Sebastian
         lrst_oce = .FALSE.
            IF( ln_rst_list ) THEN
               nrst_lst = MIN(nrst_lst + 1, SIZE(nn_stocklist,1))
               nitrst = nn_stocklist( nrst_lst )
            ENDIF
      ENDIF
      !
   END SUBROUTINE rst_write


   SUBROUTINE rst_read_open
      !!---------------------------------------------------------------------- 
      !!                   ***  ROUTINE rst_read_open  ***
      !! 
      !! ** Purpose :   Open read files for NetCDF restart
      !! 
      !! ** Method  :   Use a non-zero, positive value of numror to assess whether or not
      !!                the file has already been opened
      !!----------------------------------------------------------------------
      LOGICAL        ::   llok
      CHARACTER(lc)  ::   clpath   ! full path to ocean output restart file
      !!----------------------------------------------------------------------
      !
      IF( numror <= 0 ) THEN
         IF(lwp) THEN                                             ! Contol prints
            WRITE(numout,*)
            WRITE(numout,*) 'rst_read : read oce NetCDF restart file'
            IF ( snc4set%luse )      WRITE(numout,*) 'rst_read : configured with NetCDF4 support'
            WRITE(numout,*) '~~~~~~~~'
         ENDIF
         lxios_sini = .FALSE.
         clpath = TRIM(cn_ocerst_indir)
         IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
         CALL iom_open( TRIM(clpath)//cn_ocerst_in, numror )
! are we using XIOS to read the data? Part above will have to modified once XIOS
! can handle checking if variable is in the restart file (there will be no need to open
! restart)
         IF(.NOT.lxios_set) lrxios = lrxios.AND.lxios_sini
         IF( lrxios) THEN
             crxios_context = 'nemo_rst'
             IF( .NOT.lxios_set ) THEN
                 IF(lwp) WRITE(numout,*) 'Enable restart reading by XIOS'
                 CALL iom_init( crxios_context, ld_tmppatch = .false. )
                 lxios_set = .TRUE.
             ENDIF
         ENDIF
         IF( TRIM(Agrif_CFixed()) /= '0' .AND. lrxios) THEN
             CALL iom_init( crxios_context, ld_tmppatch = .false. )
             IF(lwp) WRITE(numout,*) 'Enable restart reading by XIOS for AGRIF'
             lxios_set = .TRUE.
         ENDIF 
      ENDIF

   END SUBROUTINE rst_read_open


   SUBROUTINE rst_read( Kbb, Kmm )
      !!---------------------------------------------------------------------- 
      !!                   ***  ROUTINE rst_read  ***
      !! 
      !! ** Purpose :   Read files for NetCDF restart
      !! 
      !! ** Method  :   Read in restart.nc file fields which are necessary for restart
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Kmm   ! ocean time level indices
      REAL(wp) ::   zrdt
      INTEGER  ::   jk
      REAL(wp), DIMENSION(jpi, jpj, jpk) :: w3d
      !!----------------------------------------------------------------------

      CALL rst_read_open           ! open restart for reading (if not already opened)

      ! Check dynamics and tracer time-step consistency and force Euler restart if changed
      IF( iom_varid( numror, 'rdt', ldstop = .FALSE. ) > 0 )   THEN
         CALL iom_get( numror, 'rdt', zrdt, ldxios = lrxios )
         IF( zrdt /= rn_Dt ) THEN
            IF(lwp) WRITE( numout,*)
            IF(lwp) WRITE( numout,*) 'rst_read:  rdt not equal to the read one'
            IF(lwp) WRITE( numout,*)
            IF(lwp) WRITE( numout,*) '      ==>>>   forced euler first time-step'
            l_1st_euler =  .TRUE.
         ENDIF
      ENDIF

      CALL iom_delay_rst( 'READ', 'OCE', numror )   ! read only ocean delayed global communication variables
      
      ! Diurnal DSST 
      IF( ln_diurnal ) CALL iom_get( numror, jpdom_autoglo, 'Dsst' , x_dsst, ldxios = lrxios ) 
      IF ( ln_diurnal_only ) THEN 
         IF(lwp) WRITE( numout, * ) &
         &   "rst_read:- ln_diurnal_only set, setting rhop=rho0" 
         rhop = rho0
         CALL iom_get( numror, jpdom_autoglo, 'tn'     , w3d, ldxios = lrxios ) 
         ts(:,:,1,jp_tem,Kmm) = w3d(:,:,1)
         RETURN 
      ENDIF  
      
      IF( iom_varid( numror, 'ub', ldstop = .FALSE. ) > 0 ) THEN
         CALL iom_get( numror, jpdom_autoglo, 'ub'     , uu(:,:,:       ,Kbb), ldxios = lrxios )   ! before fields
         CALL iom_get( numror, jpdom_autoglo, 'vb'     , vv(:,:,:       ,Kbb), ldxios = lrxios )
         CALL iom_get( numror, jpdom_autoglo, 'tb'     , ts(:,:,:,jp_tem,Kbb), ldxios = lrxios )
         CALL iom_get( numror, jpdom_autoglo, 'sb'     , ts(:,:,:,jp_sal,Kbb), ldxios = lrxios )
         CALL iom_get( numror, jpdom_autoglo, 'sshb'   ,ssh(:,:         ,Kbb), ldxios = lrxios )
      ELSE
         l_1st_euler =  .TRUE.      ! before field not found, forced euler 1st time-step
      ENDIF
      !
      CALL iom_get( numror, jpdom_autoglo, 'un'     , uu(:,:,:       ,Kmm), ldxios = lrxios )       ! now    fields
      CALL iom_get( numror, jpdom_autoglo, 'vn'     , vv(:,:,:       ,Kmm), ldxios = lrxios )
      CALL iom_get( numror, jpdom_autoglo, 'tn'     , ts(:,:,:,jp_tem,Kmm), ldxios = lrxios )
      CALL iom_get( numror, jpdom_autoglo, 'sn'     , ts(:,:,:,jp_sal,Kmm), ldxios = lrxios )
      CALL iom_get( numror, jpdom_autoglo, 'sshn'   ,ssh(:,:         ,Kmm), ldxios = lrxios )
      IF( iom_varid( numror, 'rhop', ldstop = .FALSE. ) > 0 ) THEN
         CALL iom_get( numror, jpdom_autoglo, 'rhop'   , rhop, ldxios = lrxios )   ! now    potential density
      ELSE
         CALL eos( ts(:,:,:,:,Kmm), rhd, rhop, gdept(:,:,:,Kmm) )   
      ENDIF
      !
      IF( l_1st_euler ) THEN                                  ! Euler restart 
         ts (:,:,:,:,Kbb) = ts (:,:,:,:,Kmm)                  ! all before fields set to now values
         uu (:,:,:  ,Kbb) = uu (:,:,:  ,Kmm)
         vv (:,:,:  ,Kbb) = vv (:,:,:  ,Kmm)
         ssh(:,:    ,Kbb) = ssh(:,:    ,Kmm)
      ENDIF
      !
   END SUBROUTINE rst_read

   !!=====================================================================
END MODULE restart
