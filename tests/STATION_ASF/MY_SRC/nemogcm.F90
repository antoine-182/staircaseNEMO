MODULE nemogcm
   !!======================================================================
   !!                       ***  MODULE nemogcm   ***
   !! StandAlone Surface module : surface fluxes
   !!======================================================================
   !! History :  3.6  ! 2011-11  (S. Alderson, G. Madec) original code
   !!             -   ! 2013-06  (I. Epicoco, S. Mocavero, CMCC) nemo_northcomms: setup avoiding MPI communication
   !!             -   ! 2014-12  (G. Madec) remove KPP scheme and cross-land advection (cla)
   !!            4.0  ! 2016-10  (G. Madec, S. Flavoni)  domain configuration / user defined interface
   !!            4.x  ! 2019-12  (L. Brodeau)  STATION_ASF test-case
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   nemo_gcm      : solve ocean dynamics, tracer, biogeochemistry and/or sea-ice
   !!   nemo_init     : initialization of the NEMO system
   !!   nemo_ctl      : initialisation of the contol print
   !!   nemo_closefile: close remaining open files
   !!   nemo_alloc    : dynamical allocation
   !!----------------------------------------------------------------------
   USE step_oce       ! module used in the ocean time stepping module (step.F90)
   USE sbc_oce        ! surface boundary condition: ocean #LB: rm?
   USE phycst         ! physical constant                  (par_cst routine)
   USE domain         ! domain initialization   (dom_init & dom_cfg routines)
   USE closea         ! treatment of closed seas (for ln_closea)
   USE usrdef_nam     ! user defined configuration
   USE step, ONLY : Nbb, Nnn, Naa, Nrhs ! time level indices
   USE daymod         ! calendar
   USE restart        ! open  restart file
   !LB:USE step           ! NEMO time-stepping                 (stp     routine)
   USE c1d            ! 1D configuration
   USE step_c1d       ! Time stepping loop for the 1D configuration
   USE sbcssm         !
   !
   USE lib_mpp        ! distributed memory computing
   USE mppini         ! shared/distributed memory setting (mpp_init routine)
   USE lib_fortran    ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined)
#if defined key_iomput
   USE xios           ! xIOserver
#endif

   IMPLICIT NONE
   PRIVATE

   PUBLIC   nemo_gcm    ! called by model.F90
   PUBLIC   nemo_init   ! needed by AGRIF

   CHARACTER(lc) ::   cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: nemogcm.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE nemo_gcm
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_gcm  ***
      !!
      !! ** Purpose :   NEMO solves the primitive equations on an orthogonal
      !!              curvilinear mesh on the sphere.
      !!
      !! ** Method  : - model general initialization
      !!              - launch the time-stepping (stp routine)
      !!              - finalize the run by closing files and communications
      !!
      !! References : Madec, Delecluse, Imbard, and Levy, 1997:  internal report, IPSL.
      !!              Madec, 2008, internal report, IPSL.
      !!----------------------------------------------------------------------
      INTEGER ::   istp   ! time step index
      REAL(wp)::   zstptiming   ! elapsed time for 1 time step
      !!----------------------------------------------------------------------
      !
      !                            !-----------------------!
      CALL nemo_init               !==  Initialisations  ==!
      !                            !-----------------------!
      ! check that all process are still there... If some process have an error,
      ! they will never enter in step and other processes will wait until the end of the cpu time!
      CALL mpp_max( 'nemogcm', nstop )

      IF(lwp) WRITE(numout,cform_aaa)   ! Flag AAAAAAA

      !                            !-----------------------!
      !                            !==   time stepping   ==!
      !                            !-----------------------!
      istp = nit000
      !
      DO WHILE ( istp <= nitend .AND. nstop == 0 )    !==  C1D time-stepping  ==!
         CALL stp_c1d( istp )
         istp = istp + 1
      END DO
      !
      !                            !------------------------!
      !                            !==  finalize the run  ==!
      !                            !------------------------!
      IF(lwp) WRITE(numout,cform_aaa)        ! Flag AAAAAAA
      !
      IF( nstop /= 0 .AND. lwp ) THEN        ! error print
         WRITE(ctmp1,*) '   ==>>>   nemo_gcm: a total of ', nstop, ' errors have been found'
         CALL ctl_stop( ctmp1 )
      ENDIF
      !
      IF( ln_timing )   CALL timing_finalize
      !
      CALL nemo_closefile
      !
#if defined key_iomput
      CALL xios_finalize  ! end mpp communications with xios
#else
      IF( lk_mpp   ) THEN   ;   CALL mppstop      ! end mpp communications
      ENDIF
#endif
      !
      IF(lwm) THEN
         IF( nstop == 0 ) THEN   ;   STOP 0
         ELSE                    ;   STOP 123
         ENDIF
      ENDIF
      !
   END SUBROUTINE nemo_gcm


   SUBROUTINE nemo_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_init  ***
      !!
      !! ** Purpose :   initialization of the NEMO GCM
      !!----------------------------------------------------------------------
      INTEGER ::   ios, ilocal_comm   ! local integers
      !!
      NAMELIST/namctl/ sn_cfctl,  nn_print, nn_ictls, nn_ictle,             &
         &             nn_isplt , nn_jsplt, nn_jctls, nn_jctle,             &
         &             ln_timing, ln_diacfl
      NAMELIST/namcfg/ ln_read_cfg, cn_domcfg, ln_closea, ln_write_cfg, cn_domcfg_out, ln_use_jattr
      !!----------------------------------------------------------------------
      !
      cxios_context = 'nemo'
      !
      !                             !-------------------------------------------------!
      !                             !     set communicator & select the local rank    !
      !                             !  must be done as soon as possible to get narea  !
      !                             !-------------------------------------------------!
      !
#if defined key_iomput
      IF( Agrif_Root() ) THEN
            CALL xios_initialize( "for_xios_mpi_id", return_comm=ilocal_comm )   ! nemo local communicator given by xios
      ENDIF
      CALL mpp_start( ilocal_comm )
#else
         CALL mpp_start( )
#endif
      !
      narea = mpprank + 1               ! mpprank: the rank of proc (0 --> mppsize -1 )
      lwm = (narea == 1)                ! control of output namelists
      !
      !                             !---------------------------------------------------------------!
      !                             ! Open output files, reference and configuration namelist files !
      !                             !---------------------------------------------------------------!
      !
      ! open ocean.output as soon as possible to get all output prints (including errors messages)
      IF( lwm )   CALL ctl_opn(     numout,        'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
      ! open reference and configuration namelist files
      CALL load_nml( numnam_ref,        'namelist_ref',                                           -1, lwm )
      CALL load_nml( numnam_cfg,        'namelist_cfg',                                           -1, lwm )
      IF( lwm )   CALL ctl_opn(     numond, 'output.namelist.dyn', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
      ! open /dev/null file to be able to supress output write easily
      CALL ctl_opn(     numnul,           '/dev/null', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE. )
      !
      !                             !--------------------!
      !                             ! Open listing units !  -> need sn_cfctl from namctl to define lwp
      !                             !--------------------!
      !
      READ  ( numnam_ref, namctl, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namctl in reference namelist' )
      READ  ( numnam_cfg, namctl, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namctl in configuration namelist' )
      !
      ! finalize the definition of namctl variables
      IF( sn_cfctl%l_allon ) THEN
         ! Turn on all options.
         CALL nemo_set_cfctl( sn_cfctl, .TRUE., .TRUE. )
         ! Ensure all processors are active
         sn_cfctl%procmin = 0 ; sn_cfctl%procmax = 1000000 ; sn_cfctl%procincr = 1
      ELSEIF( sn_cfctl%l_config ) THEN
         ! Activate finer control of report outputs
         ! optionally switch off output from selected areas (note this only
         ! applies to output which does not involve global communications)
         IF( ( narea < sn_cfctl%procmin .OR. narea > sn_cfctl%procmax  ) .OR. &
           & ( MOD( narea - sn_cfctl%procmin, sn_cfctl%procincr ) /= 0 ) )    &
           &   CALL nemo_set_cfctl( sn_cfctl, .FALSE., .FALSE. )
      ELSE
         ! turn off all options.
         CALL nemo_set_cfctl( sn_cfctl, .FALSE., .TRUE. )
      ENDIF
      !
      lwp = (narea == 1) .OR. sn_cfctl%l_oceout    ! control of all listing output print
      !
      IF(lwp) THEN                      ! open listing units
         !
         IF( .NOT. lwm )   &            ! alreay opened for narea == 1
            &            CALL ctl_opn( numout, 'ocean.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, -1, .FALSE., narea )
         !
         WRITE(numout,*)
         WRITE(numout,*) '   CNRS - NERC - Met OFFICE - MERCATOR-ocean - CMCC'
         WRITE(numout,*) '                       NEMO team'
         WRITE(numout,*) '            Ocean General Circulation Model'
         WRITE(numout,*) '                NEMO version 4.0  (2019) '
         WRITE(numout,*) '                      STATION_ASF '
         WRITE(numout,*) "           ._      ._      ._      ._      ._    "
         WRITE(numout,*) "       _.-._)`\_.-._)`\_.-._)`\_.-._)`\_.-._)`\_ "
         WRITE(numout,*)
         WRITE(numout,*) "           o         _,           _,             "
         WRITE(numout,*) "            o      .' (        .-' /             "
         WRITE(numout,*) "           o     _/..._'.    .'   /              "
         WRITE(numout,*) "      (    o .-'`      ` '-./  _.'               "
         WRITE(numout,*) "       )    ( o)           ;= <_         (       "
         WRITE(numout,*) "      (      '-.,\\__ __.-;`\   '.        )      "
         WRITE(numout,*) "       )  )       \) |`\ \)  '.   \      (   (   "
         WRITE(numout,*) "      (  (           \_/       '-._\      )   )  "
         WRITE(numout,*) "       )  ) jgs                     `    (   (   "
         WRITE(numout,*) "     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ "
         WRITE(numout,*)
         !
         WRITE(numout,cform_aaa)                                        ! Flag AAAAAAA
         !
      ENDIF
      !
      IF(lwm) WRITE( numond, namctl )
      !
      !                             !------------------------------------!
      !                             !  Set global domain size parameters !
      !                             !------------------------------------!
      !
      READ  ( numnam_ref, namcfg, IOSTAT = ios, ERR = 903 )
903   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namcfg in reference namelist' )
      READ  ( numnam_cfg, namcfg, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 )   CALL ctl_nam ( ios , 'namcfg in configuration namelist' )
      !
      IF( ln_read_cfg ) THEN            ! Read sizes in domain configuration file
         CALL domain_cfg ( cn_cfg, nn_cfg, jpiglo, jpjglo, jpkglo, jperio )
      ELSE                              ! user-defined namelist
         CALL usr_def_nam( cn_cfg, nn_cfg, jpiglo, jpjglo, jpkglo, jperio )
      ENDIF
      !
      IF(lwm)   WRITE( numond, namcfg )
      !
      !                             !-----------------------------------------!
      !                             ! mpp parameters and domain decomposition !
      !                             !-----------------------------------------!
      CALL mpp_init

      ! Now we know the dimensions of the grid and numout has been set: we can allocate arrays
      CALL nemo_alloc()

      ! Initialise time level indices
      Nbb = 1; Nnn = 2; Naa = 3; Nrhs = Naa

      !                             !-------------------------------!
      !                             !  NEMO general initialization  !
      !                             !-------------------------------!

      CALL nemo_ctl                          ! Control prints
      !
      !                                      ! General initialization
      IF( ln_timing    )   CALL timing_init     ! timing
      IF( ln_timing    )   CALL timing_start( 'nemo_init')
      !
      CALL     phy_cst         ! Physical constants
      CALL     eos_init        ! Equation of state
      IF( lk_c1d       )   CALL     c1d_init        ! 1D column configuration
      CALL     dom_init( Nbb, Nnn, Naa, "OPA") ! Domain
      IF( sn_cfctl%l_prtctl )   &
         &                 CALL prt_ctl_init        ! Print control

      IF( ln_rstart ) THEN                    ! Restart from a file                                                                                 
         !                                    ! -------------------                                                                                 
         CALL rst_read( Nbb, Nnn )            ! Read the restart file                                                                               
         CALL day_init                        ! model calendar (using both namelist and restart infos)                                              
         !                                                                                                                                          
      ELSE                                    ! Start from rest                                                                                     
         !                                    ! ---------------                                                                                     
         numror = 0                           ! define numror = 0 -> no restart file to read                                                        
         neuler = 0                           ! Set time-step indicator at nit000 (euler forward)                                                   
         CALL day_init                        ! model calendar (using both namelist and restart infos)                                              
      ENDIF
      !

      !                                      ! external forcing
      CALL     sbc_init( Nbb, Nnn, Naa )    ! surface boundary conditions (including sea-ice)

      !
      IF(lwp) WRITE(numout,cform_aaa)           ! Flag AAAAAAA
      !
      IF( ln_timing    )   CALL timing_stop( 'nemo_init')
      !
   END SUBROUTINE nemo_init


   SUBROUTINE nemo_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_ctl  ***
      !!
      !! ** Purpose :   control print setting
      !!
      !! ** Method  : - print namctl and namcfg information and check some consistencies
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'nemo_ctl: Control prints'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '   Namelist namctl'
         WRITE(numout,*) '                              sn_cfctl%l_glochk  = ', sn_cfctl%l_glochk
         WRITE(numout,*) '                              sn_cfctl%l_allon   = ', sn_cfctl%l_allon
         WRITE(numout,*) '       finer control over o/p sn_cfctl%l_config  = ', sn_cfctl%l_config
         WRITE(numout,*) '                              sn_cfctl%l_runstat = ', sn_cfctl%l_runstat
         WRITE(numout,*) '                              sn_cfctl%l_trcstat = ', sn_cfctl%l_trcstat
         WRITE(numout,*) '                              sn_cfctl%l_oceout  = ', sn_cfctl%l_oceout
         WRITE(numout,*) '                              sn_cfctl%l_layout  = ', sn_cfctl%l_layout
         WRITE(numout,*) '                              sn_cfctl%l_prtctl  = ', sn_cfctl%l_prtctl
         WRITE(numout,*) '                              sn_cfctl%l_prttrc  = ', sn_cfctl%l_prttrc
         WRITE(numout,*) '                              sn_cfctl%l_oasout  = ', sn_cfctl%l_oasout
         WRITE(numout,*) '                              sn_cfctl%procmin   = ', sn_cfctl%procmin
         WRITE(numout,*) '                              sn_cfctl%procmax   = ', sn_cfctl%procmax
         WRITE(numout,*) '                              sn_cfctl%procincr  = ', sn_cfctl%procincr
         WRITE(numout,*) '                              sn_cfctl%ptimincr  = ', sn_cfctl%ptimincr
         WRITE(numout,*) '      level of print                  nn_print   = ', nn_print
         WRITE(numout,*) '      Start i indice for SUM control  nn_ictls   = ', nn_ictls
         WRITE(numout,*) '      End i indice for SUM control    nn_ictle   = ', nn_ictle
         WRITE(numout,*) '      Start j indice for SUM control  nn_jctls   = ', nn_jctls
         WRITE(numout,*) '      End j indice for SUM control    nn_jctle   = ', nn_jctle
         WRITE(numout,*) '      number of proc. following i     nn_isplt   = ', nn_isplt
         WRITE(numout,*) '      number of proc. following j     nn_jsplt   = ', nn_jsplt
         WRITE(numout,*) '      timing by routine               ln_timing  = ', ln_timing
         WRITE(numout,*) '      CFL diagnostics                 ln_diacfl  = ', ln_diacfl
      ENDIF
      !
      nprint    = nn_print          ! convert DOCTOR namelist names into OLD names
      nictls    = nn_ictls
      nictle    = nn_ictle
      njctls    = nn_jctls
      njctle    = nn_jctle
      isplt     = nn_isplt
      jsplt     = nn_jsplt

      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) '   Namelist namcfg'
         WRITE(numout,*) '      read domain configuration file                ln_read_cfg      = ', ln_read_cfg
         WRITE(numout,*) '         filename to be read                           cn_domcfg     = ', TRIM(cn_domcfg)
         WRITE(numout,*) '         keep closed seas in the domain (if exist)     ln_closea     = ', ln_closea
         WRITE(numout,*) '      create a configuration definition file        ln_write_cfg     = ', ln_write_cfg
         WRITE(numout,*) '         filename to be written                        cn_domcfg_out = ', TRIM(cn_domcfg_out)
         WRITE(numout,*) '      use file attribute if exists as i/p j-start   ln_use_jattr     = ', ln_use_jattr
      ENDIF
      IF( .NOT.ln_read_cfg )   ln_closea = .false.   ! dealing possible only with a domcfg file
      !
      !                             ! Parameter control
      !
      IF( sn_cfctl%l_prtctl .OR. sn_cfctl%l_prttrc ) THEN              ! sub-domain area indices for the control prints
         IF( lk_mpp .AND. jpnij > 1 ) THEN
            isplt = jpni   ;   jsplt = jpnj   ;   ijsplt = jpni*jpnj   ! the domain is forced to the real split domain
         ELSE
            IF( isplt == 1 .AND. jsplt == 1  ) THEN
               CALL ctl_warn( ' - isplt & jsplt are equal to 1',   &
                  &           ' - the print control will be done over the whole domain' )
            ENDIF
            ijsplt = isplt * jsplt            ! total number of processors ijsplt
         ENDIF
         IF(lwp) WRITE(numout,*)'          - The total number of processors over which the'
         IF(lwp) WRITE(numout,*)'            print control will be done is ijsplt : ', ijsplt
         !
         !                              ! indices used for the SUM control
         IF( nictls+nictle+njctls+njctle == 0 )   THEN    ! print control done over the default area
            lsp_area = .FALSE.
         ELSE                                             ! print control done over a specific  area
            lsp_area = .TRUE.
            IF( nictls < 1 .OR. nictls > jpiglo )   THEN
               CALL ctl_warn( '          - nictls must be 1<=nictls>=jpiglo, it is forced to 1' )
               nictls = 1
            ENDIF
            IF( nictle < 1 .OR. nictle > jpiglo )   THEN
               CALL ctl_warn( '          - nictle must be 1<=nictle>=jpiglo, it is forced to jpiglo' )
               nictle = jpiglo
            ENDIF
            IF( njctls < 1 .OR. njctls > jpjglo )   THEN
               CALL ctl_warn( '          - njctls must be 1<=njctls>=jpjglo, it is forced to 1' )
               njctls = 1
            ENDIF
            IF( njctle < 1 .OR. njctle > jpjglo )   THEN
               CALL ctl_warn( '          - njctle must be 1<=njctle>=jpjglo, it is forced to jpjglo' )
               njctle = jpjglo
            ENDIF
         ENDIF
      ENDIF
      !
      IF( 1._wp /= SIGN(1._wp,-0._wp)  )   CALL ctl_stop( 'nemo_ctl: The intrinsec SIGN function follows f2003 standard.',  &
         &                                                'Compile with key_nosignedzero enabled:',   &
         &                                                '--> add -Dkey_nosignedzero to the definition of %CPP in your arch file' )
      !
   END SUBROUTINE nemo_ctl


   SUBROUTINE nemo_closefile
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_closefile  ***
      !!
      !! ** Purpose :   Close the files
      !!----------------------------------------------------------------------
      !
      IF( lk_mpp )   CALL mppsync
      !
      CALL iom_close                                 ! close all input/output files managed by iom_*
      !
      IF( numstp          /= -1 )   CLOSE( numstp          )   ! time-step file
      IF( numrun          /= -1 )   CLOSE( numrun          )   ! run statistics file
      IF( lwm.AND.numond  /= -1 )   CLOSE( numond          )   ! oce output namelist
      IF( lwm.AND.numoni  /= -1 )   CLOSE( numoni          )   ! ice output namelist
      IF( numevo_ice      /= -1 )   CLOSE( numevo_ice      )   ! ice variables (temp. evolution)
      IF( numout          /=  6 )   CLOSE( numout          )   ! standard model output file
      !
      numout = 6                                     ! redefine numout in case it is used after this point...
      !
   END SUBROUTINE nemo_closefile


   SUBROUTINE nemo_alloc
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_alloc  ***
      !!
      !! ** Purpose :   Allocate all the dynamic arrays of the OPA modules
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      USE diawri    , ONLY : dia_wri_alloc
      USE dom_oce   , ONLY : dom_oce_alloc
      !
      INTEGER :: ierr
      !!----------------------------------------------------------------------
      !
      ierr =        oce_alloc    ()    ! ocean
      ierr = ierr + dia_wri_alloc()
      ierr = ierr + dom_oce_alloc()    ! ocean domain
      !
      CALL mpp_sum( 'nemogcm', ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'nemo_alloc: unable to allocate standard ocean arrays' )
      !
   END SUBROUTINE nemo_alloc


   SUBROUTINE nemo_set_cfctl(sn_cfctl, setto, for_all )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE nemo_set_cfctl  ***
      !!
      !! ** Purpose :   Set elements of the output control structure to setto.
      !!                for_all should be .false. unless all areas are to be
      !!                treated identically.
      !!
      !! ** Method  :   Note this routine can be used to switch on/off some
      !!                types of output for selected areas but any output types
      !!                that involve global communications (e.g. mpp_max, glob_sum)
      !!                should be protected from selective switching by the
      !!                for_all argument
      !!----------------------------------------------------------------------
      LOGICAL :: setto, for_all
      TYPE(sn_ctl) :: sn_cfctl
      !!----------------------------------------------------------------------
      IF( for_all ) THEN
         sn_cfctl%l_runstat = setto
         sn_cfctl%l_trcstat = setto
      ENDIF
      sn_cfctl%l_oceout  = setto
      sn_cfctl%l_layout  = setto
      sn_cfctl%l_prtctl  = setto
      sn_cfctl%l_prttrc  = setto
      sn_cfctl%l_oasout  = setto
   END SUBROUTINE nemo_set_cfctl

   !!======================================================================
END MODULE nemogcm
