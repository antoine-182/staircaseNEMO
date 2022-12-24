MODULE stpRK3
   !!======================================================================
   !!                       ***  MODULE step  ***
   !! Time-stepping   : manager of the shallow water equation time stepping
   !!======================================================================
   !! History :  NEMO !  2020-03  (A. Nasser, G. Madec)  Original code from  4.0.2
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stpRK3             : Shallow Water time-stepping
   !!----------------------------------------------------------------------
   USE step_oce         ! time stepping definition modules
   USE phycst           ! physical constants
   USE usrdef_nam
   !
   USE iom              ! xIOs server 
   USE domqco

   IMPLICIT NONE
   PRIVATE

   PUBLIC   stp_RK3   ! called by nemogcm.F90
   
   !!----------------------------------------------------------------------
   !! time level indices
   !!----------------------------------------------------------------------
   INTEGER, PUBLIC ::   Nbb, Nnn, Naa, Nrhs   !! used by nemo_init
      
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: step.F90 12614 2020-03-26 14:59:52Z gm $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_agrif
   RECURSIVE SUBROUTINE stp_RK3( )
      INTEGER             ::   kstp   ! ocean time-step index
#else
   SUBROUTINE stp_RK3( kstp )
      INTEGER, INTENT(in) ::   kstp   ! ocean time-step index
#endif
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp_RK3  ***
      !!
      !! ** Purpose : - Time stepping of shallow water (SHW) (momentum and ssh eqs.)
      !!
      !! ** Method  : -1- Update forcings
      !!              -2- Update the ssh at Naa
      !!              -3- Compute the momentum trends (Nrhs)
      !!              -4- Update the horizontal velocity
      !!              -5- Apply Asselin time filter to uu,vv,ssh
      !!              -6- Outputs and diagnostics
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj, jk   ! dummy loop indice
      INTEGER ::   indic        ! error indicator if < 0
!!gm kcall can be removed, I guess
      INTEGER ::   kcall        ! optional integer argument (dom_vvl_sf_nxt)
      REAL(wp)::   z1_2rho0,  z5_6,  z3_4  ! local scalars
      
      REAL(wp) ::   zue3a, zue3n, zue3b    ! local scalars
      REAL(wp) ::   zve3a, zve3n, zve3b    !   -      -
      REAL(wp) ::   ze3t_tf, ze3u_tf, ze3v_tf, zua, zva
      !! ---------------------------------------------------------------------
#if defined key_agrif
      kstp = nit000 + Agrif_Nb_Step()
      Kbb_a = Nbb; Kmm_a = Nnn; Krhs_a = Nrhs   ! agrif_oce module copies of time level indices
      IF( lk_agrif_debug ) THEN
         IF( Agrif_Root() .and. lwp)   WRITE(*,*) '---'
         IF(lwp)   WRITE(*,*) 'Grid Number', Agrif_Fixed(),' time step ', kstp, 'int tstep', Agrif_NbStepint()
      ENDIF
      IF( kstp == nit000 + 1 )   lk_agrif_fstep = .FALSE.
# if defined key_iomput
      IF( Agrif_Nbstepint() == 0 )   CALL iom_swap( cxios_context )
# endif
#endif
      !
      IF( ln_timing )   CALL timing_start('stp_RK3')
      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! model timestep
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      IF ( kstp == nit000 )   ww(:,:,:) = 0._wp   ! initialize vertical velocity one for all to zero

      !
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! update I/O and calendar 
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                             indic = 0                ! reset to no error condition
                             
      IF( kstp == nit000 ) THEN                       ! initialize IOM context (must be done after nemo_init for AGRIF+XIOS+OASIS)
                             CALL iom_init( cxios_context, ld_closedef=.FALSE. )   ! for model grid (including passible AGRIF zoom)
         IF( lk_diamlr   )   CALL dia_mlr_iom_init    ! with additional setup for multiple-linear-regression analysis
                             CALL iom_init_closedef
         IF( ln_crs      )   CALL iom_init( TRIM(cxios_context)//"_crs" )  ! for coarse grid
      ENDIF
      IF( kstp /= nit000 )   CALL day( kstp )         ! Calendar (day was already called at nit000 in day_init)
                             CALL iom_setkt( kstp - nit000 + 1,      cxios_context          )   ! tell IOM we are at time step kstp
      IF( ln_crs         )   CALL iom_setkt( kstp - nit000 + 1, TRIM(cxios_context)//"_crs" )   ! tell IOM we are at time step kstp

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update external forcing (tides, open boundaries, ice shelf interaction and surface boundary condition (including sea-ice)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_tide    )   CALL tide_update( kstp )                     ! update tide potential
      IF( ln_apr_dyn )   CALL sbc_apr ( kstp )                        ! atmospheric pressure (NB: call before bdy_dta which needs ssh_ib) 
      IF( ln_bdy     )   CALL bdy_dta ( kstp, Nnn )                   ! update dynamic & tracer data at open boundaries
                         CALL sbc     ( kstp, Nbb, Nnn )              ! Sea Boundary Condition

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Ocean physics update
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  LATERAL  PHYSICS
      !                                                                        ! eddy diffusivity coeff.
      IF( l_ldfdyn_time )   CALL ldf_dyn( kstp, Nbb )                          ! eddy viscosity coeff. 


      !======================================================================
      !======================================================================
      !                     =====       RK3       =====
      !======================================================================
      !======================================================================

      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  RK3 1st stage Ocean dynamics : hdiv, ssh, e3, u, v, w
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      rDt   = rn_Dt / 3._wp  
      r1_Dt = 1._wp / rDt
      
                            CALL ssh_nxt       ( kstp, Nbb, Nbb, ssh, Naa )    ! after ssh (includes call to div_hor)

                         uu(:,:,:,Nrhs) = 0._wp            ! set dynamics trends to zero
                         vv(:,:,:,Nrhs) = 0._wp

!!                            CALL dyn_adv( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! advection (VF or FF)	==> RHS
 
                            CALL dyn_vor( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! vorticity           	==> RHS
#if defined key_RK3all 
                            CALL dyn_ldf( kstp, Nbb, Nbb      , uu, vv, Nrhs )  ! lateral mixing
#endif
      !
!!an - calcul du gradient de pression horizontal (explicit)
      DO_3D_00_00( 1, jpkm1 )
         uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - grav * ( ssh(ji+1,jj,Nbb) - ssh(ji,jj,Nbb) ) * r1_e1u(ji,jj)
         vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - grav * ( ssh(ji,jj+1,Nbb) - ssh(ji,jj,Nbb) ) * r1_e2v(ji,jj)
      END_3D
      !
#if defined key_RK3all 
      ! add wind stress forcing and layer linear friction to the RHS 
      z5_6 = 5._wp/6._wp
      DO_3D_00_00(1,jpkm1)
         uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) + r1_rho0 * ( z5_6*utau_b(ji,jj) + (1._wp - z5_6)*utau(ji,jj) ) / e3u(ji,jj,jk,Nbb)   &
            &                                  - rn_rfr * uu(ji,jj,jk,Nbb)
         vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) + r1_rho0 * ( z5_6*vtau_b(ji,jj) + (1._wp - z5_6)*vtau(ji,jj) ) / e3v(ji,jj,jk,Nbb)   &
            &                                  - rn_rfr * vv(ji,jj,jk,Nbb)
      END_3D
#endif
!!an
                            CALL dom_qco_r3c   ( ssh(:,:,Naa), r3t(:,:,Naa), r3u(:,:,Naa), r3v(:,:,Naa), r3f(:,:) )   ! "after" ssh./h._0 ratio explicit
      IF( ln_dynadv_vec ) THEN      ! vector invariant form : applied on velocity
         DO_3D_00_00(1,jpkm1)
            uu(ji,jj,jk,Naa) = uu(ji,jj,jk,Nbb) + rDt * uu(ji,jj,jk,Nrhs) * umask(ji,jj,jk)
            vv(ji,jj,jk,Naa) = vv(ji,jj,jk,Nbb) + rDt * vv(ji,jj,jk,Nrhs) * vmask(ji,jj,jk)
         END_3D          
      ELSE
         DO_3D_00_00(1,jpkm1)       ! flux form : applied on thickness weighted velocity
            uu(ji,jj,jk,Naa) = (         uu(ji,jj,jk,Nbb )*e3u(ji,jj,jk,Nbb)                              &
               &                 + rDt * uu(ji,jj,jk,Nrhs)*e3t(ji,jj,jk,Nbb) * umask(ji,jj,jk)        )   &
               &               /                           e3t(ji,jj,jk,Naa)
            vv(ji,jj,jk,Naa) = (         vv(ji,jj,jk,Nbb )*e3v(ji,jj,jk,Nbb)                              &
               &                 + rDt * vv(ji,jj,jk,Nrhs)*e3t(ji,jj,jk,Nbb) * vmask(ji,jj,jk)        )   &
               &               /                           e3t(ji,jj,jk,Naa)
         END_3D
      ENDIF
      
      CALL lbc_lnk_multi( 'stp_RK3', uu(:,:,:,Nnn), 'U', -1., vv(:,:,:,Nnn), 'V', -1.,   &   !* local domain boundaries
             &                       uu(:,:,:,Naa), 'U', -1., vv(:,:,:,Naa), 'V', -1.    )     


      ! Swap time levels
      Nrhs= Nnn
      Nnn = Naa
      Naa = Nrhs

      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  RK3 2nd stage Ocean dynamics : hdiv, ssh, e3, u, v, w
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      rDt   = rn_Dt / 2._wp  
      r1_Dt = 1._wp / rDt
      
                            CALL ssh_nxt       ( kstp, Nbb, Nnn, ssh, Naa )    ! after ssh (includes call to div_hor)

                         uu(:,:,:,Nrhs) = 0._wp            ! set dynamics trends to zero
                         vv(:,:,:,Nrhs) = 0._wp
!!st TBC for dyn_adv
!!                            CALL dyn_adv( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! advection (VF or FF)	==> RHS
 
                            CALL dyn_vor( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! vorticity           	==> RHS
#if defined key_RK3all  
                            CALL dyn_ldf( kstp, Nbb, Nbb      , uu, vv, Nrhs )  ! lateral mixing
#endif
                            
      !
!!an - calcul du gradient de pression horizontal (explicit)
      DO_3D_00_00( 1, jpkm1 )
         uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - grav * ( ssh(ji+1,jj,Nnn) - ssh(ji,jj,Nnn) ) * r1_e1u(ji,jj)
         vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - grav * ( ssh(ji,jj+1,Nnn) - ssh(ji,jj,Nnn) ) * r1_e2v(ji,jj)
      END_3D
      !
      ! add wind stress forcing and layer linear friction to the RHS 
#if defined key_RK3all
      z3_4 = 3._wp/4._wp
      DO_3D_00_00(1,jpkm1)
         uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) + r1_rho0 * ( z3_4*utau_b(ji,jj) + (1._wp - z3_4)*utau(ji,jj) ) / e3u(ji,jj,jk,Nbb)   &
            &                                  - rn_rfr * uu(ji,jj,jk,Nbb)
         vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) + r1_rho0 * ( z3_4*vtau_b(ji,jj) + (1._wp - z3_4)*vtau(ji,jj) ) / e3v(ji,jj,jk,Nbb)   &
            &                                  - rn_rfr * vv(ji,jj,jk,Nbb)
      END_3D
#endif
!!an
                           CALL dom_qco_r3c   ( ssh(:,:,Naa), r3t(:,:,Naa), r3u(:,:,Naa), r3v(:,:,Naa), r3f(:,:) )   ! "after" ssh./h._0 ratio explicit
      IF( ln_dynadv_vec ) THEN      ! vector invariant form : applied on velocity
         DO_3D_00_00(1,jpkm1)
            uu(ji,jj,jk,Naa) = uu(ji,jj,jk,Nbb) + rDt * uu(ji,jj,jk,Nrhs) * umask(ji,jj,jk)
            vv(ji,jj,jk,Naa) = vv(ji,jj,jk,Nbb) + rDt * vv(ji,jj,jk,Nrhs) * vmask(ji,jj,jk)
         END_3D          
      ELSE
         DO_3D_00_00(1,jpkm1)       ! flux form : applied on thickness weighted velocity
            uu(ji,jj,jk,Naa) = (         uu(ji,jj,jk,Nbb )*e3u(ji,jj,jk,Nbb)                              &
               &                 + rDt * uu(ji,jj,jk,Nrhs)*e3t(ji,jj,jk,Nnn) * umask(ji,jj,jk)        )   &
               &               /                           e3t(ji,jj,jk,Naa)
            vv(ji,jj,jk,Naa) = (         vv(ji,jj,jk,Nbb )*e3v(ji,jj,jk,Nbb)                              &
               &                 + rDt * vv(ji,jj,jk,Nrhs)*e3t(ji,jj,jk,Nnn) * vmask(ji,jj,jk)        )   &
               &               /                           e3t(ji,jj,jk,Naa)
         END_3D
      ENDIF
      
      CALL lbc_lnk_multi( 'stp_RK3', uu(:,:,:,Nnn), 'U', -1., vv(:,:,:,Nnn), 'V', -1.,   &   !* local domain boundaries
             &                       uu(:,:,:,Naa), 'U', -1., vv(:,:,:,Naa), 'V', -1.    )     

      ! Swap time levels
      Nrhs= Nnn
      Nnn = Naa
      Naa = Nrhs
       
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  RK3 3rd stage Ocean dynamics : hdiv, ssh, e3, u, v, w
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      rDt   = rn_Dt
      r1_Dt = 1._wp / rDt

                            CALL ssh_nxt       ( kstp, Nbb, Nnn, ssh, Naa )    ! after ssh (includes call to div_hor)
      
                         uu(:,:,:,Nrhs) = 0._wp            ! set dynamics trends to zero
                         vv(:,:,:,Nrhs) = 0._wp

      IF( ln_bdy     )      CALL bdy_dyn3d_dmp ( kstp, Nbb,      uu, vv, Nrhs )  ! bdy damping trends

#if defined key_agrif
      IF(.NOT. Agrif_Root())  & 
               &            CALL Agrif_Sponge_dyn        ! momentum sponge
#endif
!!                            CALL dyn_adv( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! advection (VF or FF)	==> RHS
 
                            CALL dyn_vor( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! vorticity           	==> RHS
 
                            CALL dyn_ldf( kstp, Nbb, Nnn      , uu, vv, Nrhs )  ! lateral mixing

!!an - calcul du gradient de pression horizontal (explicit)
      DO_3D_00_00( 1, jpkm1 )
         uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) - grav * ( ssh(ji+1,jj,Nnn) - ssh(ji,jj,Nnn) ) * r1_e1u(ji,jj)
         vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) - grav * ( ssh(ji,jj+1,Nnn) - ssh(ji,jj,Nnn) ) * r1_e2v(ji,jj)
      END_3D
      !
      ! add wind stress forcing and layer linear friction to the RHS 
      z1_2rho0 = 0.5_wp * r1_rho0
      DO_3D_00_00(1,jpkm1)
         uu(ji,jj,jk,Nrhs) = uu(ji,jj,jk,Nrhs) + z1_2rho0 * ( utau_b(ji,jj) + utau(ji,jj) ) / e3u(ji,jj,jk,Nnn)   &
            &                                  - rn_rfr * uu(ji,jj,jk,Nbb)
         vv(ji,jj,jk,Nrhs) = vv(ji,jj,jk,Nrhs) + z1_2rho0 * ( vtau_b(ji,jj) + vtau(ji,jj) ) / e3v(ji,jj,jk,Nnn)   &
            &                                  - rn_rfr * vv(ji,jj,jk,Nbb)
      END_3D   
!!an         
                            CALL dom_qco_r3c   ( ssh(:,:,Naa), r3t(:,:,Naa), r3u(:,:,Naa), r3v(:,:,Naa), r3f(:,:) )   ! "after" ssh./h._0 ratio explicit      
      IF( ln_dynadv_vec ) THEN      ! vector invariant form : applied on velocity
         DO_3D_11_11(1,jpkm1)
            zua = uu(ji,jj,jk,Nbb) + rDt * uu(ji,jj,jk,Nrhs) * umask(ji,jj,jk)
            zva = vv(ji,jj,jk,Nbb) + rDt * vv(ji,jj,jk,Nrhs) * vmask(ji,jj,jk)
            !                                                                  ! Asselin time filter on u,v (Nnn)
            uu(ji,jj,jk,Nnn) = uu(ji,jj,jk,Nnn) + rn_atfp * (uu(ji,jj,jk,Nbb) - 2._wp * uu(ji,jj,jk,Nnn) + zua)
            vv(ji,jj,jk,Nnn) = vv(ji,jj,jk,Nnn) + rn_atfp * (vv(ji,jj,jk,Nbb) - 2._wp * vv(ji,jj,jk,Nnn) + zva)
            !              
            uu(ji,jj,jk,Naa) = zua
            vv(ji,jj,jk,Naa) = zva
         END_3D
         !
      ELSE                          ! flux form : applied on thickness weighted velocity
         DO_3D_11_11(1,jpkm1)
            zue3n = e3u(ji,jj,jk,Nnn) * uu(ji,jj,jk,Nnn)
            zve3n = e3v(ji,jj,jk,Nnn) * vv(ji,jj,jk,Nnn)
            zue3b = e3u(ji,jj,jk,Nbb) * uu(ji,jj,jk,Nbb)
            zve3b = e3v(ji,jj,jk,Nbb) * vv(ji,jj,jk,Nbb)
            !                                                ! LF time stepping
            zue3a = zue3b + rDt * e3t(ji,jj,jk,Nbb) * uu(ji,jj,jk,Nrhs) * umask(ji,jj,jk)
            zve3a = zve3b + rDt * e3t(ji,jj,jk,Nbb) * vv(ji,jj,jk,Nrhs) * vmask(ji,jj,jk)
            !
            uu(ji,jj,jk,Naa) = zue3a / e3t(ji,jj,jk,Naa)    
            vv(ji,jj,jk,Naa) = zve3a / e3t(ji,jj,jk,Naa)
         END_3D
!!st je ne comprends pas l'histoire des e3t et du Nbb et pas du Nnn pour le rhs ?        
      ENDIF


      CALL lbc_lnk_multi( 'stp_RK3', uu(:,:,:,Nnn), 'U', -1., vv(:,:,:,Nnn), 'V', -1.,   &   !* local domain boundaries
             &                       uu(:,:,:,Naa), 'U', -1., vv(:,:,:,Naa), 'V', -1.    )     
      ! Swap time levels
      Nrhs = Nbb
      Nbb = Naa
      Naa = Nrhs
      !
!                         CALL dom_vvl_sf_update_st( kstp, Nbb, Nnn, Naa )  ! recompute vertical scale factors
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! diagnostics and outputs
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_floats  )   CALL flo_stp   ( kstp, Nbb, Nnn )      ! drifting Floats
      IF( ln_diacfl  )   CALL dia_cfl   ( kstp,      Nnn )      ! Courant number diagnostics
    
                         CALL dia_wri   ( kstp,      Nnn )      ! ocean model: outputs

      !
      IF( lrst_oce   )   CALL rst_write    ( kstp, Nbb, Nnn )   ! write output ocean restart file
          

#if defined key_agrif
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! AGRIF
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
                         Kbb_a = Nbb; Kmm_a = Nnn; Krhs_a = Nrhs      ! agrif_oce module copies of time level indices
                         CALL Agrif_Integrate_ChildGrids( stp_RK3 )       ! allows to finish all the Child Grids before updating

                         IF( Agrif_NbStepint() == 0 ) THEN
                            CALL Agrif_update_all( )                  ! Update all components
                         ENDIF
#endif
      IF( ln_diaobs  )   CALL dia_obs      ( kstp, Nnn )      ! obs-minus-model (assimilation) diagnostics (call after dynamics update)

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                         CALL stp_ctl      ( kstp, Nbb, Nnn, indic )
               
                         
      IF( kstp == nit000 ) THEN                          ! 1st time step only
                                        CALL iom_close( numror )   ! close input  ocean restart file
         IF(lwm)                        CALL FLUSH    ( numond )   ! flush output namelist oce
         IF(lwm .AND. numoni /= -1 )    CALL FLUSH    ( numoni )   ! flush output namelist ice (if exist)
      ENDIF

      !
#if defined key_iomput
      IF( kstp == nitend .OR. indic < 0 ) THEN 
                      CALL iom_context_finalize(      cxios_context          ) ! needed for XIOS+AGRIF
                      IF(lrxios) CALL iom_context_finalize(      crxios_context          )
      ENDIF
#endif
      !
      IF( l_1st_euler ) THEN         ! recover Leap-frog timestep
         rDt = 2._wp * rn_Dt   
         r1_Dt = 1._wp / rDt
         l_1st_euler = .FALSE.      
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('stp_RK3')
      !
   END SUBROUTINE stp_RK3
   !
   !!======================================================================
END MODULE stpRK3
