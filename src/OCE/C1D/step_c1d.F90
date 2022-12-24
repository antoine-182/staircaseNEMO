MODULE step_c1d
   !!======================================================================
   !!                       ***  MODULE step_c1d  ***
   !! Time-stepping    : manager of the ocean, tracer and ice time stepping - c1d case
   !!======================================================================
   !! History :   2.0  !  2004-04  (C. Ethe)  adapted from step.F90 for C1D
   !!             3.0  !  2008-04  (G. Madec)  redo the adaptation to include SBC
   !!             4.1  !  2019-08  (A. Coward, D. Storkey) rewrite in preparation for new timestepping scheme
   !!----------------------------------------------------------------------
#if defined key_c1d
   !!----------------------------------------------------------------------
   !!   'key_c1d'                                       1D Configuration
   !!----------------------------------------------------------------------  
   !!   stp_c1d        : NEMO system time-stepping in c1d case
   !!----------------------------------------------------------------------
   USE step_oce        ! time stepping definition modules 
   USE step, ONLY : Nbb, Nnn, Naa, Nrhs ! time level indices
#if defined key_top
   USE trcstp          ! passive tracer time-stepping      (trc_stp routine)
#endif
   USE dyncor_c1d      ! Coriolis term (c1d case)         (dyn_cor_1d     )
   USE dynatf          ! time filtering                   (dyn_atf routine)
   USE dyndmp          ! U & V momentum damping           (dyn_dmp routine)
   USE restart         ! restart 

   IMPLICIT NONE
   PRIVATE

   PUBLIC stp_c1d      ! called by opa.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: step_c1d.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE stp_c1d( kstp )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp_c1d  ***
      !!                      
      !! ** Purpose :  - Time stepping of SBC including sea ice (dynamic and thermodynamic eqs.)
      !!               - Time stepping of OPA (momentum and active tracer eqs.)
      !!               - Time stepping of TOP (passive tracer eqs.)
      !! 
      !! ** Method  : -1- Update forcings and data  
      !!              -2- Update vertical ocean physics 
      !!              -3- Compute the t and s trends 
      !!              -4- Update t and s 
      !!              -5- Compute the momentum trends
      !!              -6- Update the horizontal velocity
      !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w)
      !!              -8- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
      !
      INTEGER ::   jk       ! dummy loop indice
      INTEGER ::   indic    ! error indicator if < 0
      !! ---------------------------------------------------------------------

                             indic = 0                ! reset to no error condition
      IF( kstp == nit000 )   CALL iom_init( "nemo")   ! iom_put initialization (must be done after nemo_init for AGRIF+XIOS+OASIS)
      IF( kstp /= nit000 )   CALL day( kstp )         ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp - nit000 + 1, "nemo" )   ! say to iom that we are at time step kstp

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update data, open boundaries, surface boundary condition (including sea-ice)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL sbc    ( kstp, Nbb, Nnn )  ! Sea Boundary Condition (including sea-ice)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Ocean physics update        
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL eos_rab( ts(:,:,:,:,Nbb), rab_b, Nnn )  ! before local thermal/haline expension ratio at T-points
                         CALL eos_rab( ts(:,:,:,:,Nnn), rab_n, Nnn )  ! now    local thermal/haline expension ratio at T-points
                         CALL bn2( ts(:,:,:,:,Nbb), rab_b, rn2b, Nnn ) ! before Brunt-Vaisala frequency
                         CALL bn2( ts(:,:,:,:,Nnn), rab_n, rn2 , Nnn ) ! now    Brunt-Vaisala frequency
      
      !  VERTICAL PHYSICS
                         CALL zdf_phy( kstp, Nbb, Nnn, Nrhs  )    ! vertical physics update (bfr, avt, avs, avm + MLD)

      IF(.NOT.ln_linssh )   CALL ssh_nxt       ( kstp, Nbb, Nnn, ssh, Naa )  ! after ssh (includes call to div_hor)
      IF(.NOT.ln_linssh )   CALL dom_vvl_sf_nxt( kstp, Nbb, Nnn,      Naa )  ! after vertical scale factors 

      IF(.NOT.ln_linssh )   CALL wzv           ( kstp, Nbb, Nnn, Naa, ww )  ! now cross-level velocity 
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! diagnostics and outputs       
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL dia_wri( kstp, Nnn )  ! ocean model: outputs
      IF( lk_diahth  )   CALL dia_hth( kstp, Nnn )  ! Thermocline depth (20°C)


#if defined key_top
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Passive Tracer Model
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        CALL trc_stp( kstp, Nbb, Nnn, Nrhs, Naa  )   ! time-stepping
#endif

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Active tracers                              (uu(:,:,:,Nrhs), vv(:,:,:,Nrhs) used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        ts(:,:,:,:,Nrhs) = 0._wp       ! set tracer trends to zero

                        CALL tra_sbc( kstp,      Nnn, ts, Nrhs  )  ! surface boundary condition
      IF( ln_traqsr )   CALL tra_qsr( kstp,      Nnn, ts, Nrhs  )  ! penetrative solar radiation qsr
      IF( ln_tradmp )   CALL tra_dmp( kstp, Nbb, Nnn, ts, Nrhs  )  ! internal damping trends- tracers
      IF(.NOT.ln_linssh)CALL tra_adv( kstp, Nbb, Nnn, ts, Nrhs  )  ! horizontal & vertical advection
      IF( ln_zdfosm  )  CALL tra_osm( kstp, Nnn     , ts, Nrhs  )  ! OSMOSIS non-local tracer fluxes
                        CALL tra_zdf( kstp, Nbb, Nnn, Nrhs, ts, Naa   )         ! vertical mixing
                        CALL eos( ts(:,:,:,:,Nnn), rhd, rhop, gdept_0(:,:,:) )  ! now potential density for zdfmxl
      IF( ln_zdfnpc )   CALL tra_npc( kstp,      Nnn, Nrhs, ts, Naa   )         ! applied non penetrative convective adjustment on (t,s)
                        CALL tra_atf( kstp, Nbb, Nnn, Nrhs,     Naa, ts   )     ! time filtering of "now" tracer fields



      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Dynamics                                    (ta, sa used as workspace)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        uu(:,:,:,Nrhs) = 0._wp          ! set dynamics trends to zero
                        vv(:,:,:,Nrhs) = 0._wp

      IF( ln_dyndmp )   CALL dyn_dmp    ( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! internal damping trends- momentum
                        CALL dyn_cor_c1d( kstp,      Nnn      , uu, vv, Nrhs )  ! vorticity term including Coriolis
      IF( ln_zdfosm  )  CALL dyn_osm    ( kstp,      Nnn      , uu, vv, Nrhs )  ! OSMOSIS non-local velocity fluxes
                        CALL dyn_zdf    ( kstp, Nbb, Nnn, Nrhs, uu, vv, Naa  )  ! vertical diffusion
                        CALL dyn_atf    ( kstp, Nbb, Nnn, Naa , uu, vv, e3t, e3u, e3v )  ! time filtering of "now" fields
      IF(.NOT.ln_linssh)CALL ssh_atf    ( kstp, Nbb, Nnn, Naa , ssh )                    ! time filtering of "now" sea surface height
      !
      ! Swap time levels
      Nrhs = Nbb
      Nbb = Nnn
      Nnn = Naa
      Naa = Nrhs
      !
      IF(.NOT.ln_linssh)CALL dom_vvl_sf_update( kstp, Nbb, Nnn, Naa )                    ! update of vertical scale factors

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control and restarts
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                             CALL stp_ctl( kstp, Nnn, indic )
      IF( kstp == nit000 )   CALL iom_close( numror )          ! close input  ocean restart file
      IF( lrst_oce       )   CALL rst_write( kstp, Nbb, Nnn )  ! write output ocean restart file
      !
#if defined key_iomput
      IF( kstp == nitend .OR. indic < 0 )   CALL xios_context_finalize()   ! needed for XIOS
      !
#endif
   END SUBROUTINE stp_c1d

#else
   !!----------------------------------------------------------------------
   !!   Default key                                            NO 1D Config
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE stp_c1d ( kt )      ! dummy routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'stp_c1d: You should not have seen this print! error?', kt
   END SUBROUTINE stp_c1d
#endif

   !!======================================================================
END MODULE step_c1d
