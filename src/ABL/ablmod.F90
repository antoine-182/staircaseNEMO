MODULE ablmod
   !!======================================================================
   !!                       ***  MODULE  ablmod  ***
   !! Surface module :  ABL computation to provide atmospheric data 
   !!                   for surface fluxes computation
   !!======================================================================
   !! History :  3.6  ! 2019-03  (F. Lemarié & G. Samson)  Original code
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   abl_stp       : ABL single column model
   !!   abl_zdf_tke   : atmospheric vertical closure
   !!----------------------------------------------------------------------
   USE abl            ! ABL dynamics and tracers
   USE par_abl        ! ABL constants

   USE phycst         ! physical constants
   USE dom_oce, ONLY  : tmask  
   USE sbc_oce, ONLY  : ght_abl, ghw_abl, e3t_abl, e3w_abl, jpka, jpkam1, rhoa
   USE sbcblk         ! use rn_?fac
   USE sbcblk_phy     ! use some physical constants for flux computation
   !
   USE prtctl         ! Print control                    (prt_ctl routine)
   USE iom            ! IOM library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE

   PUBLIC   abl_stp   ! called by sbcabl.F90
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ustar2
   !! * Substitutions
#  include "do_loop_substitute.h90"

CONTAINS


!===================================================================================================
   SUBROUTINE abl_stp( kt, psst, pssu, pssv, pssq, &                            ! in
              &            pu_dta, pv_dta, pt_dta, pq_dta,    & 
              &            pslp_dta, pgu_dta, pgv_dta,        &
              &            pcd_du, psen, pevp,                &                            ! in/out
              &            pwndm, ptaui, ptauj, ptaum         & 
#if defined key_si3			  
              &          , ptm_su,pssu_ice,pssv_ice,pssq_ice,pcd_du_ice  & 
              &          , psen_ice, pevp_ice, pwndm_ice, pfrac_oce      & 
              &          , ptaui_ice, ptauj_ice                          &
#endif			  
			  &   )
!---------------------------------------------------------------------------------------------------

      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE abl_stp ***
      !!
      !! ** Purpose :   Time-integration of the ABL model 
      !!
      !! ** Method  :   Compute atmospheric variables : vertical turbulence 
      !!                             + Coriolis term + newtonian relaxation
      !!                
      !! ** Action  : - Advance TKE to time n+1 and compute Avm_abl, Avt_abl, PBLh
      !!              - Advance tracers to time n+1 (Euler backward scheme)
      !!              - Compute Coriolis term with forward-backward scheme (possibly with geostrophic guide)
      !!              - Advance u,v to time n+1 (Euler backward scheme)
      !!              - Apply newtonian relaxation on the dynamics and the tracers
      !!              - Finalize flux computation in psen, pevp, pwndm, ptaui, ptauj, ptaum
      !!
      !!----------------------------------------------------------------------
      INTEGER  , INTENT(in   )                   ::   kt         ! time step index
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   psst       ! sea-surface temperature [Celsius]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pssu       ! sea-surface u (U-point)
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pssv       ! sea-surface v (V-point) 
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pssq       ! sea-surface humidity
      REAL(wp) , INTENT(in   ), DIMENSION(:,:,:) ::   pu_dta     ! large-scale windi
      REAL(wp) , INTENT(in   ), DIMENSION(:,:,:) ::   pv_dta     ! large-scale windj
      REAL(wp) , INTENT(in   ), DIMENSION(:,:,:) ::   pgu_dta    ! large-scale hpgi
      REAL(wp) , INTENT(in   ), DIMENSION(:,:,:) ::   pgv_dta    ! large-scale hpgj
      REAL(wp) , INTENT(in   ), DIMENSION(:,:,:) ::   pt_dta     ! large-scale pot. temp.
      REAL(wp) , INTENT(in   ), DIMENSION(:,:,:) ::   pq_dta     ! large-scale humidity
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pslp_dta   ! sea-level pressure
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pcd_du     ! Cd x Du (T-point)
      REAL(wp) , INTENT(inout), DIMENSION(:,:  ) ::   psen       ! Ch x Du
      REAL(wp) , INTENT(inout), DIMENSION(:,:  ) ::   pevp       ! Ce x Du
      REAL(wp) , INTENT(inout), DIMENSION(:,:  ) ::   pwndm      ! ||uwnd||      
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ) ::   ptaui      ! taux
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ) ::   ptauj      ! tauy
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ) ::   ptaum      ! ||tau|| 
      !
#if defined key_si3	  
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   ptm_su       ! ice-surface temperature [K]
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pssu_ice     ! ice-surface u (U-point)
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pssv_ice     ! ice-surface v (V-point) 	  
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pssq_ice     ! ice-surface humidity 
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pcd_du_ice   ! Cd x Du over ice (T-point)
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   psen_ice     ! Ch x Du over ice (T-point)
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pevp_ice     ! Ce x Du over ice (T-point)
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pwndm_ice    ! ||uwnd - uice||  
      !REAL(wp) , INTENT(inout), DIMENSION(:,:  ) ::   pfrac_oce	   !!GS: out useless ?
      REAL(wp) , INTENT(in   ), DIMENSION(:,:  ) ::   pfrac_oce	   !
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ) ::   ptaui_ice    ! ice-surface taux stress (U-point)
      REAL(wp) , INTENT(  out), DIMENSION(:,:  ) ::   ptauj_ice    ! ice-surface tauy stress (V-point) 	  
#endif	  
      !
      REAL(wp), DIMENSION(1:jpi,1:jpj   )        ::   zwnd_i, zwnd_j
      REAL(wp), DIMENSION(1:jpi,2:jpka  )        ::   zCF    
      REAL(wp), DIMENSION(1:jpi,1:jpj,1:jpka)    ::   z_cft      !--FL--to be removed after the test phase   
      !
      REAL(wp), DIMENSION(1:jpi,1:jpka  )        ::   z_elem_a
      REAL(wp), DIMENSION(1:jpi,1:jpka  )        ::   z_elem_b
      REAL(wp), DIMENSION(1:jpi,1:jpka  )        ::   z_elem_c
      !
      INTEGER             ::   ji, jj, jk, jtra, jbak               ! dummy loop indices
      REAL(wp)            ::   zztmp, zcff, ztemp, zhumi, zcff1, zztmp1, zztmp2
      REAL(wp)            ::   zcff2, zfcor, zmsk, zsig, zcffu, zcffv, zzice,zzoce
      !
      !!---------------------------------------------------------------------      
      !
      IF(lwp .AND. kt == nit000) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'abl_stp : ABL time stepping'
         WRITE(numout,*) '~~~~~~'
      ENDIF 
      !
      IF( kt == nit000 ) ALLOCATE ( ustar2( 1:jpi, 1:jpj ) )
      !! Compute ustar squared as Cd || Uatm-Uoce ||^2 
      !! needed for surface boundary condition of TKE 
      !! pwndm contains | U10m - U_oce | (see blk_oce_1 in sbcblk)
      DO_2D_11_11
         zzoce         = pCd_du    (ji,jj) * pwndm    (ji,jj)
#if defined key_si3
         zzice         = pCd_du_ice(ji,jj) * pwndm_ice(ji,jj) 
         ustar2(ji,jj) = zzoce * pfrac_oce(ji,jj) + (1._wp - pfrac_oce(ji,jj)) * zzice 
#else
         ustar2(ji,jj) = zzoce	
#endif
      END_2D
      !
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  1 *** Advance TKE to time n+1 and compute Avm_abl, Avt_abl, PBLh
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      CALL abl_zdf_tke( )                       !--> Avm_abl, Avt_abl, pblh defined on (1,jpi) x (1,jpj) 
         
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  2 *** Advance tracers to time n+1
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !-------------
      DO jj = 1, jpj    ! outer loop             !--> tq_abl computed on (1:jpi) x (1:jpj)
      !-------------    
         ! Compute matrix elements for interior points 
         DO jk = 3, jpkam1
            DO ji = 1, jpi   ! vector opt.
               z_elem_a( ji,     jk              ) = - rDt_abl * Avt_abl( ji, jj, jk-1 ) / e3w_abl( jk-1 )   ! lower-diagonal
               z_elem_c( ji,     jk              ) = - rDt_abl * Avt_abl( ji, jj, jk   ) / e3w_abl( jk   )   ! upper-diagonal       
               z_elem_b( ji,     jk              ) = e3t_abl(jk) - z_elem_a( ji, jk ) - z_elem_c( ji, jk )   !       diagonal           
            END DO
         END DO  
         ! Boundary conditions  
         DO ji = 1, jpi   ! vector opt. 
            ! Neumann at the bottom          
            z_elem_a( ji,     2              ) = 0._wp
            z_elem_c( ji,     2              ) = - rDt_abl * Avt_abl( ji, jj, 2   ) / e3w_abl( 2   )                                         
            ! Homogeneous Neumann at the top
            z_elem_a( ji,     jpka           ) = - rDt_abl * Avt_abl( ji, jj, jpka ) / e3w_abl( jpka ) 
            z_elem_c( ji,     jpka           ) = 0._wp
            z_elem_b( ji,     jpka           ) = e3t_abl( jpka ) - z_elem_a( ji,     jpka )
         END DO            

         DO jtra = 1,jptq  ! loop on active tracers
               
            DO jk = 3, jpkam1
               DO ji = 1,jpi
                  tq_abl  ( ji, jj, jk, nt_a, jtra ) = e3t_abl(jk) * tq_abl  ( ji, jj, jk, nt_n, jtra )   ! initialize right-hand-side
               END DO
            END DO

            IF(jtra == jp_ta) THEN
               DO ji = 1,jpi  ! boundary conditions for temperature              
                  zztmp1 = psen(ji, jj) 
                  zztmp2 = psen(ji, jj) * ( psst(ji, jj) + rt0 )	
#if defined key_si3				  
                  zztmp1 = zztmp1 * pfrac_oce(ji,jj) + (1._wp - pfrac_oce(ji,jj)) * psen_ice(ji,jj) 
                  zztmp2 = zztmp2 * pfrac_oce(ji,jj) + (1._wp - pfrac_oce(ji,jj)) * psen_ice(ji,jj) * ptm_su(ji,jj)
#endif				  
                  z_elem_b( ji,     2                ) = e3t_abl( 2 ) - z_elem_c( ji, 2 ) + rDt_abl * zztmp1				  
                  tq_abl  ( ji, jj, 2   , nt_a, jtra ) = e3t_abl( 2    ) * tq_abl  ( ji, jj, 2   , nt_n, jtra ) + rDt_abl * zztmp2               
                  tq_abl  ( ji, jj, jpka, nt_a, jtra ) = e3t_abl( jpka ) * tq_abl  ( ji, jj, jpka, nt_n, jtra )
               END DO 
            ELSE   
               DO ji = 1,jpi  ! boundary conditions for humidity              
                  zztmp1 = pevp(ji, jj) 
                  zztmp2 = pevp(ji, jj) * pssq(ji, jj)   
#if defined key_si3				  
                  zztmp1 = zztmp1 * pfrac_oce(ji,jj) + (1._wp - pfrac_oce(ji,jj)) * pevp_ice(ji,jj) 
                  zztmp2 = zztmp2 * pfrac_oce(ji,jj) + (1._wp - pfrac_oce(ji,jj)) * pevp_ice(ji, jj) * pssq_ice(ji, jj)    
#endif	
                  z_elem_b( ji,     2                ) = e3t_abl( 2 ) - z_elem_c( ji, 2 ) + rDt_abl * zztmp1
                  tq_abl  ( ji, jj, 2   , nt_a, jtra ) = e3t_abl( 2    ) * tq_abl  ( ji, jj, 2   , nt_n, jtra ) + rDt_abl * zztmp2               
                  tq_abl  ( ji, jj, jpka, nt_a, jtra ) = e3t_abl( jpka ) * tq_abl  ( ji, jj, jpka, nt_n, jtra )
               END DO                                     
            END IF
            !!
            !! Matrix inversion
            !! ----------------------------------------------------------
            DO ji = 1,jpi            
               zcff                       =  1._wp / z_elem_b( ji, 2 )
               zCF   (ji,   2           ) = - zcff * z_elem_c( ji, 2 )
               tq_abl(ji,jj,2,nt_a,jtra) =    zcff * tq_abl(ji,jj,2,nt_a,jtra)
            END DO             

            DO jk = 3, jpka            
               DO ji = 1,jpi            
                  zcff = 1._wp / ( z_elem_b( ji, jk ) + z_elem_a( ji, jk ) * zCF   (ji, jk-1 ) )
                  zCF(ji,jk) = - zcff * z_elem_c( ji, jk )
                  tq_abl(ji,jj,jk,nt_a,jtra) = zcff * ( tq_abl(ji,jj,jk  ,nt_a,jtra)   &
                  &                - z_elem_a(ji, jk) * tq_abl(ji,jj,jk-1,nt_a,jtra) )
               END DO
            END DO
!!FL at this point we could check positivity of tq_abl(:,:,:,nt_a,jp_qa) ... test to do ...         
            DO jk = jpkam1,2,-1            
               DO ji = 1,jpi  
                  tq_abl(ji,jj,jk,nt_a,jtra) = tq_abl(ji,jj,jk,nt_a,jtra) +    &
                     &                        zCF(ji,jk) * tq_abl(ji,jj,jk+1,nt_a,jtra)
               END DO
            END DO
         
         END DO   !<-- loop on tracers                    
         !! 
      !-------------
      END DO             ! end outer loop 
      !-------------   
      
      
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  3 *** Compute Coriolis term with geostrophic guide
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
      !-------------        
      DO jk = 2, jpka    ! outer loop
      !-------------
         !             
         ! Advance u_abl & v_abl to time n+1
         DO_2D_11_11
            zcff = ( fft_abl(ji,jj) * rDt_abl )*( fft_abl(ji,jj) * rDt_abl )  ! (f dt)**2
   
            u_abl( ji, jj, jk, nt_a ) = e3t_abl(jk) *(  &
               &        (1._wp-gamma_Cor*(1._wp-gamma_Cor)*zcff)*u_abl( ji, jj, jk, nt_n )    &
               &                 +  rDt_abl * fft_abl(ji, jj) * v_abl ( ji , jj  , jk, nt_n ) )  &
               &                               / (1._wp + gamma_Cor*gamma_Cor*zcff)
               
            v_abl( ji, jj, jk, nt_a ) =  e3t_abl(jk) *(  &
               &        (1._wp-gamma_Cor*(1._wp-gamma_Cor)*zcff)*v_abl( ji, jj, jk, nt_n )   &
               &                 -  rDt_abl * fft_abl(ji, jj) * u_abl ( ji   , jj, jk, nt_n )  ) &
               &                                / (1._wp + gamma_Cor*gamma_Cor*zcff)                
         END_2D
         !                                   
      !-------------
      END DO             ! end outer loop  !<-- u_abl and v_abl are properly updated on (1:jpi) x (1:jpj)
      !-------------                
      !
      IF( ln_geos_winds ) THEN
         DO jj = 1, jpj    ! outer loop       
            DO jk = 1, jpka
               DO ji = 1, jpi 
                  u_abl( ji, jj, jk, nt_a ) = u_abl( ji, jj, jk, nt_a )   &
                     &                      - rDt_abl * e3t_abl(jk) * fft_abl(ji  , jj) * pgv_dta(ji  ,jj  ,jk)
                  v_abl( ji, jj, jk, nt_a ) = v_abl( ji, jj, jk, nt_a )   &
                     &                      + rDt_abl * e3t_abl(jk) * fft_abl(ji, jj  ) * pgu_dta(ji  ,jj  ,jk)
               END DO
            END DO
         END DO      
      END IF 
      !-------------            
      !
      IF( ln_hpgls_frc ) THEN
         DO jj = 1, jpj    ! outer loop       
            DO jk = 1, jpka
               DO ji = 1, jpi 
                  u_abl( ji, jj, jk, nt_a ) = u_abl( ji, jj, jk, nt_a ) - rDt_abl * e3t_abl(jk) * pgu_dta(ji,jj,jk)
                  v_abl( ji, jj, jk, nt_a ) = v_abl( ji, jj, jk, nt_a ) - rDt_abl * e3t_abl(jk) * pgv_dta(ji,jj,jk)  
               ENDDO
            ENDDO
         ENDDO 
      END IF
	  
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  4 *** Advance u,v to time n+1 
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
      !
      !  Vertical diffusion for u_abl     
      !-------------
      DO jj = 1, jpj    ! outer loop
      !-------------    

         DO jk = 3, jpkam1
            DO ji = 1, jpi  
               z_elem_a( ji,     jk ) = - rDt_abl * Avm_abl( ji, jj, jk-1 ) / e3w_abl( jk-1 )  ! lower-diagonal
               z_elem_c( ji,     jk ) = - rDt_abl * Avm_abl( ji, jj, jk   ) / e3w_abl( jk   )  ! upper-diagonal                
               z_elem_b( ji,     jk ) = e3t_abl(jk) - z_elem_a( ji, jk ) - z_elem_c( ji, jk )                             !       diagonal
            END DO
         END DO  
            
         DO ji = 2, jpi   ! boundary conditions   (Avm_abl and pcd_du must be available at ji=jpi)            
            !++ Surface boundary condition
            z_elem_a( ji,     2    ) = 0._wp
            z_elem_c( ji,     2    ) = - rDt_abl * Avm_abl( ji, jj, 2   ) / e3w_abl( 2   )                                       
            !
			zztmp1  = pcd_du(ji, jj)
			zztmp2  = 0.5_wp * pcd_du(ji, jj) * ( pssu(ji-1, jj) + pssu(ji,jj) )			
#if defined key_si3				  
            zztmp1 = zztmp1 * pfrac_oce(ji,jj) + (1._wp - pfrac_oce(ji,jj)) * pcd_du_ice(ji, jj)
			zzice  = 0.5_wp * ( pssu_ice(ji-1, jj) + pssu_ice(ji,jj) )
			zztmp2 = zztmp2 * pfrac_oce(ji,jj) + (1._wp - pfrac_oce(ji,jj)) * pcd_du_ice(ji, jj) * zzice
#endif           
			z_elem_b( ji,     2       ) = e3t_abl( 2 ) - z_elem_c( ji, 2 ) + rDt_abl * zztmp1 			
            u_abl( ji, jj,    2, nt_a ) =      u_abl( ji, jj,    2, nt_a ) + rDt_abl * zztmp2
			
            !++ Top Neumann B.C.
            !z_elem_a( ji,     jpka ) = - 0.5_wp * rDt_abl * ( Avm_abl( ji, jj, jpka )+ Avm_abl( ji+1, jj, jpka ) ) / e3w_abl( jpka ) 
            !z_elem_c( ji,     jpka ) = 0._wp
            !z_elem_b( ji,     jpka ) = e3t_abl( jpka ) - z_elem_a( ji,     jpka )                                               
            !++ Top Dirichlet B.C.
            z_elem_a( ji,     jpka )       = 0._wp
            z_elem_c( ji,     jpka )       = 0._wp
            z_elem_b( ji,     jpka )       = e3t_abl( jpka )
            u_abl   ( ji, jj, jpka, nt_a ) = e3t_abl( jpka ) * pu_dta(ji,jj,jk)              
         END DO 
         !!
         !! Matrix inversion
         !! ----------------------------------------------------------
         DO ji = 2, jpi          
            zcff                 =   1._wp / z_elem_b( ji, 2 )
            zCF   (ji,   2     ) =  - zcff * z_elem_c( ji, 2 ) 
            u_abl (ji,jj,2,nt_a) =    zcff * u_abl(ji,jj,2,nt_a)
         END DO             

         DO jk = 3, jpka            
            DO ji = 2, jpi              
               zcff = 1._wp / ( z_elem_b( ji, jk ) + z_elem_a( ji, jk ) * zCF   (ji, jk-1 ) )
               zCF(ji,jk) = - zcff * z_elem_c( ji, jk )
               u_abl(ji,jj,jk,nt_a) = zcff * ( u_abl(ji,jj,jk  ,nt_a)   &
               &          - z_elem_a(ji, jk) * u_abl(ji,jj,jk-1,nt_a) )
            END DO
         END DO
            
         DO jk = jpkam1,2,-1            
            DO ji = 2, jpi
               u_abl(ji,jj,jk,nt_a) = u_abl(ji,jj,jk,nt_a) + zCF(ji,jk) * u_abl(ji,jj,jk+1,nt_a)
            END DO
         END DO
                                          
      !-------------
      END DO             ! end outer loop 
      !-------------   

      !
      !  Vertical diffusion for v_abl     
      !-------------
      DO jj = 2, jpj   ! outer loop
      !-------------    
         !
         DO jk = 3, jpkam1
            DO ji = 1, jpi   
               z_elem_a( ji,     jk ) = -rDt_abl * Avm_abl( ji, jj, jk-1 ) / e3w_abl( jk-1 )   ! lower-diagonal
               z_elem_c( ji,     jk ) = -rDt_abl * Avm_abl( ji, jj, jk   ) / e3w_abl( jk   )   ! upper-diagonal              
               z_elem_b( ji,     jk ) = e3t_abl(jk) - z_elem_a( ji, jk ) - z_elem_c( ji, jk )                              !       diagonal
            END DO
         END DO

         DO ji = 1, jpi   ! boundary conditions (Avm_abl and pcd_du must be available at jj=jpj)          
            !++ Surface boundary condition
            z_elem_a( ji,     2    ) = 0._wp
            z_elem_c( ji,     2    ) = - rDt_abl * Avm_abl( ji, jj, 2   ) / e3w_abl( 2   )        
            !
			zztmp1 = pcd_du(ji, jj)
			zztmp2 = 0.5_wp * pcd_du(ji, jj) * ( pssv(ji, jj) + pssv(ji, jj-1) )  
#if defined key_si3				  
            zztmp1 = zztmp1 * pfrac_oce(ji,jj) + (1._wp - pfrac_oce(ji,jj)) * pcd_du_ice(ji, jj)
			zzice  = 0.5_wp * ( pssv_ice(ji, jj) + pssv_ice(ji,jj-1) )
			zztmp2 = zztmp2 * pfrac_oce(ji,jj) + (1._wp - pfrac_oce(ji,jj)) * pcd_du_ice(ji, jj) * zzice
#endif 			
            z_elem_b( ji,     2       ) = e3t_abl( 2 ) - z_elem_c( ji, 2 ) + rDt_abl * zztmp1 
            v_abl( ji, jj,    2, nt_a ) =         v_abl( ji, jj, 2, nt_a ) + rDt_abl * zztmp2
            !++ Top Neumann B.C.            
            !z_elem_a( ji,     jpka ) = -rDt_abl * Avm_abl( ji, jj, jpka ) / e3w_abl( jpka ) 
            !z_elem_c( ji,     jpka ) = 0._wp
            !z_elem_b( ji,     jpka ) = e3t_abl( jpka ) - z_elem_a( ji,     jpka )                                               
            !++ Top Dirichlet B.C.
            z_elem_a( ji,     jpka )       = 0._wp
            z_elem_c( ji,     jpka )       = 0._wp
            z_elem_b( ji,     jpka )       = e3t_abl( jpka )
            v_abl   ( ji, jj, jpka, nt_a ) = e3t_abl( jpka ) * pv_dta(ji,jj,jk)             
         END DO 
         !!
         !! Matrix inversion
         !! ----------------------------------------------------------
         DO ji = 1, jpi              
            zcff                 =  1._wp / z_elem_b( ji, 2 )
            zCF   (ji,   2     ) =   - zcff * z_elem_c( ji,     2       ) 
            v_abl (ji,jj,2,nt_a) =     zcff * v_abl   ( ji, jj, 2, nt_a ) 
         END DO             

         DO jk = 3, jpka            
            DO ji = 1, jpi                   
               zcff = 1._wp / ( z_elem_b( ji, jk ) + z_elem_a( ji, jk ) * zCF   (ji, jk-1 ) )
               zCF(ji,jk) = - zcff * z_elem_c( ji, jk )
               v_abl(ji,jj,jk,nt_a) = zcff * ( v_abl(ji,jj,jk  ,nt_a)   &
               &          - z_elem_a(ji, jk) * v_abl(ji,jj,jk-1,nt_a) )
            END DO
         END DO
            
         DO jk = jpkam1,2,-1            
            DO ji = 1, jpi    
               v_abl(ji,jj,jk,nt_a) = v_abl(ji,jj,jk,nt_a) + zCF(ji,jk) * v_abl(ji,jj,jk+1,nt_a)
            END DO
         END DO
         !                                          
      !-------------
      END DO             ! end outer loop 
      !-------------   
    
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  5 *** Apply nudging on the dynamics and the tracers 
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
      z_cft(:,:,:) = 0._wp      
      
      IF( nn_dyn_restore > 0  ) THEN
         !------------- 
         DO jk = 2, jpka    ! outer loop
         !-------------       
            DO_2D_01_01
               zcff1 = pblh( ji, jj )
               zsig  = ght_abl(jk) / MAX( jp_pblh_min,  MIN(  jp_pblh_max, zcff1  ) )                        
               zsig  =               MIN( jp_bmax    ,  MAX(         zsig, jp_bmin) ) 
               zmsk  = msk_abl(ji,jj)
               zcff2 = jp_alp3_dyn * zsig**3 + jp_alp2_dyn * zsig**2   &
                  &  + jp_alp1_dyn * zsig    + jp_alp0_dyn
               zcff  = (1._wp-zmsk) + zmsk * zcff2 * rn_Dt   ! zcff = 1 for masked points
                                                             ! rn_Dt = rDt_abl / nn_fsbc                          
               zcff  = zcff * rest_eq(ji,jj)
               z_cft( ji, jj, jk ) = zcff
               u_abl( ji, jj, jk, nt_a ) = (1._wp - zcff ) *  u_abl( ji, jj, jk, nt_a )   &
                  &                               + zcff   * pu_dta( ji, jj, jk       )                      
               v_abl( ji, jj, jk, nt_a ) = (1._wp - zcff ) *  v_abl( ji, jj, jk, nt_a )   &
                  &                               + zcff   * pv_dta( ji, jj, jk       )
            END_2D
         !-------------
         END DO             ! end outer loop
         !-------------               
      END IF

      !------------- 
      DO jk = 2, jpka    ! outer loop
      !-------------       
         DO_2D_11_11
            zcff1 = pblh( ji, jj )
            zsig  = ght_abl(jk) / MAX( jp_pblh_min,  MIN(  jp_pblh_max, zcff1  ) )
            zsig  =               MIN( jp_bmax    ,  MAX(         zsig, jp_bmin) ) 
            zmsk  = msk_abl(ji,jj)
            zcff2 = jp_alp3_tra * zsig**3 + jp_alp2_tra * zsig**2   &
               &  + jp_alp1_tra * zsig    + jp_alp0_tra
            zcff  = (1._wp-zmsk) + zmsk * zcff2 * rn_Dt   ! zcff = 1 for masked points
                                                          ! rn_Dt = rDt_abl / nn_fsbc                          
            !z_cft( ji, jj, jk ) = zcff
            tq_abl( ji, jj, jk, nt_a, jp_ta ) = (1._wp - zcff ) * tq_abl( ji, jj, jk, nt_a, jp_ta )   &
               &                                       + zcff   * pt_dta( ji, jj, jk )
            
            tq_abl( ji, jj, jk, nt_a, jp_qa ) = (1._wp - zcff ) * tq_abl( ji, jj, jk, nt_a, jp_qa )   &
               &                                       + zcff   * pq_dta( ji, jj, jk )
            
         END_2D
      !-------------
      END DO             ! end outer loop
      !------------- 	  
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  6 *** MPI exchanges  
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      CALL lbc_lnk_multi( 'ablmod',  u_abl(:,:,:,nt_a      ), 'T', -1.,  v_abl(:,:,:,nt_a      ), 'T', -1. )
      CALL lbc_lnk_multi( 'ablmod', tq_abl(:,:,:,nt_a,jp_ta), 'T',  1., tq_abl(:,:,:,nt_a,jp_qa), 'T',  1., kfillmode = jpfillnothing )   ! ++++ this should not be needed...
      !
      ! first ABL level
      IF ( iom_use("uz1_abl") ) CALL iom_put ( "uz1_abl",   u_abl(:,:,2,nt_a      ) )
      IF ( iom_use("vz1_abl") ) CALL iom_put ( "vz1_abl",   v_abl(:,:,2,nt_a      ) )
      IF ( iom_use("tz1_abl") ) CALL iom_put ( "tz1_abl",  tq_abl(:,:,2,nt_a,jp_ta) )
      IF ( iom_use("qz1_abl") ) CALL iom_put ( "qz1_abl",  tq_abl(:,:,2,nt_a,jp_qa) )
      IF ( iom_use("uz1_dta") ) CALL iom_put ( "uz1_dta",  pu_dta(:,:,2           ) )
      IF ( iom_use("vz1_dta") ) CALL iom_put ( "vz1_dta",  pv_dta(:,:,2           ) )
      IF ( iom_use("tz1_dta") ) CALL iom_put ( "tz1_dta",  pt_dta(:,:,2           ) )
      IF ( iom_use("qz1_dta") ) CALL iom_put ( "qz1_dta",  pq_dta(:,:,2           ) )
      ! all ABL levels
      IF ( iom_use("u_abl"  ) ) CALL iom_put ( "u_abl"  ,   u_abl(:,:,2:jpka,nt_a      ) )
      IF ( iom_use("v_abl"  ) ) CALL iom_put ( "v_abl"  ,   v_abl(:,:,2:jpka,nt_a      ) )
      IF ( iom_use("t_abl"  ) ) CALL iom_put ( "t_abl"  ,  tq_abl(:,:,2:jpka,nt_a,jp_ta) )
      IF ( iom_use("q_abl"  ) ) CALL iom_put ( "q_abl"  ,  tq_abl(:,:,2:jpka,nt_a,jp_qa) )
      IF ( iom_use("tke_abl") ) CALL iom_put ( "tke_abl", tke_abl(:,:,2:jpka,nt_a      ) )
      IF ( iom_use("avm_abl") ) CALL iom_put ( "avm_abl", avm_abl(:,:,2:jpka           ) )
      IF ( iom_use("avt_abl") ) CALL iom_put ( "avt_abl", avm_abl(:,:,2:jpka           ) )
      IF ( iom_use("mxl_abl") ) CALL iom_put ( "mxl_abl", mxl_abl(:,:,2:jpka           ) )
      IF ( iom_use("pblh"   ) ) CALL iom_put (  "pblh"  ,    pblh(:,:                  ) )
      ! debug (to be removed)
      IF ( iom_use("u_dta") ) CALL iom_put ( "u_dta",  pu_dta(:,:,2:jpka) )
      IF ( iom_use("v_dta") ) CALL iom_put ( "v_dta",  pv_dta(:,:,2:jpka) )
      IF ( iom_use("t_dta") ) CALL iom_put ( "t_dta",  pt_dta(:,:,2:jpka) )
      IF ( iom_use("q_dta") ) CALL iom_put ( "q_dta",  pq_dta(:,:,2:jpka) )
      IF ( iom_use("coeft") ) CALL iom_put ( "coeft",   z_cft(:,:,2:jpka) )
      IF( ln_geos_winds ) THEN
         IF ( iom_use("uz1_geo") ) CALL iom_put ( "uz1_geo", pgu_dta(:,:,2           ) )
         IF ( iom_use("vz1_geo") ) CALL iom_put ( "vz1_geo", pgv_dta(:,:,2           ) )
      END IF
      IF( ln_hpgls_frc ) THEN
         IF ( iom_use("uz1_geo") ) CALL iom_put ( "uz1_geo",  pgu_dta(:,:,2)/MAX(fft_abl(:,:),2.5e-5_wp)   )
         IF ( iom_use("vz1_geo") ) CALL iom_put ( "vz1_geo", -pgv_dta(:,:,2)/MAX(fft_abl(:,:),2.5e-5_wp)   )
      END IF
      !
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  7 *** Finalize flux computation
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

      DO_2D_11_11
         ztemp             = tq_abl  ( ji, jj, 2, nt_a, jp_ta ) 
         zhumi             = tq_abl  ( ji, jj, 2, nt_a, jp_qa ) 
         !zcff              = pslp_dta( ji, jj ) /   &              !<-- At this point ztemp and zhumi should not be zero ...
         !   &                        (  R_dry*ztemp * ( 1._wp + rctv0*zhumi )  )
         zcff              = rho_air( ztemp, zhumi, pslp_dta( ji, jj ) )
         psen ( ji, jj )   =      cp_air(zhumi) * zcff * psen(ji,jj) * ( psst(ji,jj) + rt0 - ztemp )
         pevp ( ji, jj )   = rn_efac*MAX( 0._wp,  zcff * pevp(ji,jj) * ( pssq(ji,jj)       - zhumi ) )
         rhoa( ji, jj )   = zcff              
      END_2D
      
      DO_2D_01_01
         zwnd_i(ji,jj) = u_abl(ji  ,jj,2,nt_a) - 0.5_wp * rn_vfac * ( pssu(ji  ,jj) + pssu(ji-1,jj) )  
         zwnd_j(ji,jj) = v_abl(ji,jj  ,2,nt_a) - 0.5_wp * rn_vfac * ( pssv(ji,jj  ) + pssv(ji,jj-1) ) 
      END_2D
      ! 
      CALL lbc_lnk_multi( 'ablmod', zwnd_i(:,:) , 'T', -1., zwnd_j(:,:) , 'T', -1. )
      !
      ! ... scalar wind ( = | U10m - U_oce | ) at T-point (masked)
      DO_2D_11_11
         zcff          = SQRT(  zwnd_i(ji,jj) * zwnd_i(ji,jj)   &
            &                 + zwnd_j(ji,jj) * zwnd_j(ji,jj)  )  ! * msk_abl(ji,jj)
         zztmp         = rhoa(ji,jj) * pcd_du(ji,jj)
         
         pwndm (ji,jj) =         zcff
         ptaum (ji,jj) = zztmp * zcff
         zwnd_i(ji,jj) = zztmp * zwnd_i(ji,jj)
         zwnd_j(ji,jj) = zztmp * zwnd_j(ji,jj)
      END_2D
      ! ... utau, vtau at U- and V_points, resp.
      !     Note the use of 0.5*(2-umask) in order to unmask the stress along coastlines
      !     Note the use of MAX(tmask(i,j),tmask(i+1,j) is to mask tau over ice shelves
      DO_2D_00_00
         zcff  = 0.5_wp * ( 2._wp - msk_abl(ji,jj)*msk_abl(ji+1,jj) )
         zztmp = MAX(msk_abl(ji,jj),msk_abl(ji+1,jj))
         ptaui(ji,jj) = zcff * zztmp * ( zwnd_i(ji,jj) + zwnd_i(ji+1,jj  ) )
         zcff  = 0.5_wp * ( 2._wp - msk_abl(ji,jj)*msk_abl(ji,jj+1) )
         zztmp = MAX(msk_abl(ji,jj),msk_abl(ji,jj+1))
         ptauj(ji,jj) = zcff * zztmp * ( zwnd_j(ji,jj) + zwnd_j(ji  ,jj+1) )
      END_2D
      !
      CALL lbc_lnk_multi( 'ablmod', ptaui(:,:), 'U', -1., ptauj(:,:), 'V', -1. )

      CALL iom_put( "taum_oce", ptaum )

      IF(sn_cfctl%l_prtctl) THEN
         CALL prt_ctl( tab2d_1=pwndm  , clinfo1=' abl_stp: wndm   : ' )
         CALL prt_ctl( tab2d_1=ptaui  , clinfo1=' abl_stp: utau   : ' )
         CALL prt_ctl( tab2d_2=ptauj  , clinfo2=          'vtau   : ' )
      ENDIF

#if defined key_si3
         ! ------------------------------------------------------------ !
         !    Wind stress relative to the moving ice ( U10m - U_ice )   !
         ! ------------------------------------------------------------ !
         DO_2D_00_00
            
            zztmp1 = 0.5_wp * ( u_abl(ji+1,jj,2,nt_a) + u_abl(ji,jj,2,nt_a) )
            zztmp2 = 0.5_wp * ( v_abl(ji,jj+1,2,nt_a) + v_abl(ji,jj,2,nt_a) )
   
            ptaui_ice(ji,jj) = 0.5_wp * (  rhoa(ji+1,jj) * pCd_du_ice(ji+1,jj)             &
               &                      +    rhoa(ji  ,jj) * pCd_du_ice(ji  ,jj)  )          &
               &         * ( zztmp1 - rn_vfac * pssu_ice(ji,jj) )
            ptauj_ice(ji,jj) = 0.5_wp * (  rhoa(ji,jj+1) * pCd_du_ice(ji,jj+1)             &
               &                      +    rhoa(ji,jj  ) * pCd_du_ice(ji,jj  )  )          &
               &         * ( zztmp2 - rn_vfac * pssv_ice(ji,jj) )
         END_2D
         CALL lbc_lnk_multi( 'ablmod', ptaui_ice, 'U', -1., ptauj_ice, 'V', -1. )
         !
         IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab2d_1=ptaui_ice  , clinfo1=' abl_stp: putaui : '   &
            &                                , tab2d_2=ptauj_ice  , clinfo2='          pvtaui : ' )
#endif
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  8 *** Swap time indices for the next timestep
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
      nt_n = 1 + MOD( kt  , 2)
      nt_a = 1 + MOD( kt+1, 2)
      !	  
!---------------------------------------------------------------------------------------------------
   END SUBROUTINE abl_stp
!===================================================================================================

















!===================================================================================================
   SUBROUTINE abl_zdf_tke( )
!---------------------------------------------------------------------------------------------------

      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE abl_zdf_tke  ***
      !!
      !! ** Purpose :   Time-step Turbulente Kinetic Energy (TKE) equation
      !!
      !! ** Method  : - source term due to shear
      !!              - source term due to stratification
      !!              - resolution of the TKE equation by inverting
      !!                a tridiagonal linear system
      !!
      !! ** Action  : - en : now turbulent kinetic energy)
      !!              - avmu, avmv : production of TKE by shear at u and v-points
      !!                (= Kz dz[Ub] * dz[Un] )
      !! ---------------------------------------------------------------------
      INTEGER                                 ::   ji, jj, jk, tind, jbak, jkup, jkdwn 
      INTEGER, DIMENSION(1:jpi          )     ::   ikbl
      REAL(wp)                                ::   zcff, zcff2, ztken, zesrf, zetop, ziRic, ztv
      REAL(wp)                                ::   zdU, zdV, zcff1,zshear,zbuoy,zsig, zustar2
      REAL(wp)                                ::   zdU2,zdV2      
      REAL(wp)                                ::   zwndi,zwndj
      REAL(wp), DIMENSION(1:jpi,      1:jpka) ::   zsh2
      REAL(wp), DIMENSION(1:jpi,1:jpj,1:jpka) ::   zbn2
      REAL(wp), DIMENSION(1:jpi,1:jpka  )     ::   zFC, zRH, zCF       
      REAL(wp), DIMENSION(1:jpi,1:jpka  )     ::   z_elem_a
      REAL(wp), DIMENSION(1:jpi,1:jpka  )     ::   z_elem_b
      REAL(wp), DIMENSION(1:jpi,1:jpka  )     ::   z_elem_c
      LOGICAL                                 ::   ln_Patankar    = .FALSE.
      LOGICAL                                 ::   ln_dumpvar     = .FALSE.
      LOGICAL , DIMENSION(1:jpi         )     ::   ln_foundl 
      !
      tind  = nt_n
      ziRic = 1._wp / rn_Ric
      ! if tind = nt_a it is required to apply lbc_lnk on u_abl(nt_a) and v_abl(nt_a)
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  Advance TKE equation to time n+1
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !-------------
      DO jj = 1, jpj    ! outer loop
      !-------------
         !
         ! Compute vertical shear         
         DO jk = 2, jpkam1
            DO ji = 1,jpi  
               zcff        = 1.0_wp / e3w_abl( jk )**2                
               zdU         = zcff* Avm_abl(ji,jj,jk) * (u_abl( ji, jj, jk+1, tind)-u_abl( ji, jj, jk  , tind) )**2                           
               zdV         = zcff* Avm_abl(ji,jj,jk) * (v_abl( ji, jj, jk+1, tind)-v_abl( ji, jj, jk  , tind) )**2
               zsh2(ji,jk) = zdU+zdV
            END DO
         END DO   
         !
         ! Compute brunt-vaisala frequency
         DO jk = 2, jpkam1
            DO ji = 1,jpi 
               zcff  = grav * itvref / e3w_abl( jk )  
               zcff1 =  tq_abl( ji, jj, jk+1, tind, jp_ta) - tq_abl( ji, jj, jk  , tind, jp_ta)
               zcff2 =  tq_abl( ji, jj, jk+1, tind, jp_ta) * tq_abl( ji, jj, jk+1, tind, jp_qa)        &
                  &   - tq_abl( ji, jj, jk  , tind, jp_ta) * tq_abl( ji, jj, jk  , tind, jp_qa)
               zbn2(ji,jj,jk) = zcff * ( zcff1 + rctv0 * zcff2 )  !<-- zbn2 defined on (2,jpi)
            END DO
         END DO 
         !
         ! Terms for the tridiagonal problem
         DO jk = 2, jpkam1
            DO ji = 1,jpi 
               zshear       =                 zsh2( ji,     jk )   ! zsh2 is already multiplied by Avm_abl at this point
               zsh2(ji,jk)  = zsh2( ji, jk ) / Avm_abl( ji, jj, jk )   ! reformulate zsh2 as a 'true' vertical shear for PBLH computation          
               zbuoy        = - Avt_abl( ji, jj, jk ) * zbn2( ji, jj, jk ) 
               
               z_elem_a( ji,     jk )   = - 0.5_wp * rDt_abl * rn_Sch * ( Avm_abl( ji, jj, jk   )+Avm_abl( ji, jj, jk-1 ) ) / e3t_abl( jk   ) ! lower-diagonal
               z_elem_c( ji,     jk )   = - 0.5_wp * rDt_abl * rn_Sch * ( Avm_abl( ji, jj, jk   )+Avm_abl( ji, jj, jk+1 ) ) / e3t_abl( jk+1 ) ! upper-diagonal          
               IF( (zbuoy + zshear) .gt. 0.) THEN    ! Patankar trick to avoid negative values of TKE
                  z_elem_b( ji,     jk )   = e3w_abl(jk) - z_elem_a( ji, jk ) - z_elem_c( ji, jk )   &
                     &                     + e3w_abl(jk) * rDt_abl * rn_Ceps * sqrt(tke_abl( ji, jj, jk, nt_n )) / mxl_abl(ji,jj,jk)     ! diagonal       
                  tke_abl( ji, jj, jk, nt_a )  = e3w_abl(jk) * ( tke_abl( ji, jj, jk, nt_n ) + rDt_abl * ( zbuoy + zshear ) )             ! right-hand-side
               ELSE
                  z_elem_b( ji,     jk )   = e3w_abl(jk) - z_elem_a( ji, jk ) - z_elem_c( ji, jk )   &
                     &                     + e3w_abl(jk) * rDt_abl * rn_Ceps * sqrt(tke_abl( ji, jj, jk, nt_n )) / mxl_abl(ji,jj,jk)   &  ! diagonal    
                     &                     - e3w_abl(jk) * rDt_abl * zbuoy   
                  tke_abl( ji, jj, jk, nt_a )  = e3w_abl(jk) * ( tke_abl( ji, jj, jk, nt_n ) + rDt_abl *  zshear )             ! right-hand-side                     
               END IF
            END DO
         END DO  
                            
         DO ji = 1,jpi    ! vector opt.            
            zesrf   =  MAX( 4.63_wp * ustar2(ji,jj), tke_min )           
            zetop   = tke_min             
            z_elem_a ( ji,     1    ) = 0._wp; z_elem_c ( ji,     1    ) = 0._wp; z_elem_b ( ji,     1    ) = 1._wp                                                   
            z_elem_a ( ji,     jpka ) = 0._wp; z_elem_c ( ji,     jpka ) = 0._wp; z_elem_b ( ji,     jpka ) = 1._wp            
            tke_abl( ji, jj,    1, nt_a ) = zesrf
            tke_abl( ji, jj, jpka, nt_a ) = zetop 
            zbn2(ji,jj,   1) = zbn2( ji,jj,   2) 
            zsh2(ji,      1) = zsh2( ji,      2)                                      
            zbn2(ji,jj,jpka) = zbn2( ji,jj,jpkam1) 
            zsh2(ji,   jpka) = zsh2( ji  , jpkam1)
         END DO        
         !!
         !! Matrix inversion
         !! ----------------------------------------------------------
         DO ji = 1,jpi
            zcff                  =  1._wp / z_elem_b( ji, 1 )
            zCF    (ji,   1     ) = - zcff * z_elem_c( ji,     1       ) 
            tke_abl(ji,jj,1,nt_a) =   zcff * tke_abl ( ji, jj, 1, nt_a ) 
         END DO             

         DO jk = 2, jpka            
            DO ji = 1,jpi
               zcff = 1._wp / ( z_elem_b( ji, jk ) + z_elem_a( ji, jk ) * zCF(ji, jk-1 ) )
               zCF(ji,jk) = - zcff * z_elem_c( ji, jk )
               tke_abl(ji,jj,jk,nt_a) =   zcff * ( tke_abl(ji,jj,jk  ,nt_a)   &
               &          - z_elem_a(ji, jk) * tke_abl(ji,jj,jk-1,nt_a) )
            END DO
         END DO
            
         DO jk = jpkam1,1,-1            
            DO ji = 1,jpi
               tke_abl(ji,jj,jk,nt_a) = tke_abl(ji,jj,jk,nt_a) + zCF(ji,jk) * tke_abl(ji,jj,jk+1,nt_a)
            END DO
         END DO
         
!!FL should not be needed because of Patankar procedure  
         tke_abl(2:jpi,jj,1:jpka,nt_a) = MAX( tke_abl(2:jpi,jj,1:jpka,nt_a), tke_min )

         !!
         !! Diagnose PBL height
         !! ----------------------------------------------------------
         
         
         !                                                       
         ! arrays zRH, zFC and zCF are available at this point
         ! and zFC(:, 1 ) = 0.
         ! diagnose PBL height based on zsh2 and zbn2
         zFC (  :  ,1) = 0._wp
         ikbl( 1:jpi ) = 0 
         
         DO jk = 2,jpka
            DO ji = 1, jpi 
               zcff  = ghw_abl( jk-1 )
               zcff1 = zcff / ( zcff + rn_epssfc * pblh ( ji, jj ) )
               zcff  = ghw_abl( jk   )
               zcff2 = zcff / ( zcff + rn_epssfc * pblh ( ji, jj ) )
               zFC( ji, jk ) = zFC( ji, jk-1) + 0.5_wp * e3t_abl( jk )*(                 &
                               zcff2 * ( zsh2( ji, jk  ) - ziRic * zbn2( ji, jj, jk   ) &
                           - rn_Cek  * ( fft_abl( ji, jj  ) * fft_abl( ji, jj ) ) ) &
                             + zcff1 * ( zsh2( ji, jk-1) - ziRic * zbn2( ji, jj, jk-1 ) &
                           - rn_Cek  * ( fft_abl( ji, jj  ) * fft_abl( ji, jj ) ) ) &
                           &                                                 )
               IF( ikbl(ji) == 0 .and. zFC( ji, jk ).lt.0._wp ) ikbl(ji)=jk
            END DO
         END DO
         !
         ! finalize the computation of the PBL height
         DO ji = 1, jpi
            jk = ikbl(ji)
            IF( jk > 2 ) THEN ! linear interpolation to get subgrid value of pblh
               pblh( ji, jj ) =  (  ghw_abl( jk-1 ) * zFC( ji, jk   )       &
                  &              -  ghw_abl( jk   ) * zFC( ji, jk-1 )       &
                  &              ) / ( zFC( ji, jk   ) - zFC( ji, jk-1 ) )
            ELSE IF( jk==2 ) THEN
               pblh( ji, jj ) = ghw_abl(2   )
            ELSE
               pblh( ji, jj ) = ghw_abl(jpka)
            END IF 
         END DO
      !-------------      
      END DO      
      !-------------
      !
      ! Optional : could add pblh smoothing if pblh is noisy horizontally ... 
      IF(ln_smth_pblh) THEN
         CALL lbc_lnk( 'ablmod', pblh, 'T', 1.)
         CALL smooth_pblh( pblh, msk_abl )
         CALL lbc_lnk( 'ablmod', pblh, 'T', 1.)	  
      ENDIF
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  Diagnostic mixing length computation
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !
      SELECT CASE ( nn_amxl )
      !
      CASE ( 0 )           ! Deardroff 80 length-scale bounded by the distance to surface and bottom         
#   define zlup zRH      
#   define zldw zFC 
         DO jj = 1, jpj     ! outer loop
            !
            DO ji = 1, jpi
               mxl_abl ( ji, jj,    1 )  = 0._wp
               mxl_abl ( ji, jj, jpka )  = mxl_min
               zldw( ji,        1 )  = 0._wp              
               zlup( ji,     jpka )  = 0._wp
            END DO
            !
            DO jk = 2, jpkam1   
               DO ji = 1, jpi  
                  zbuoy             = MAX( zbn2(ji, jj, jk), rsmall )  
                  mxl_abl( ji, jj, jk ) = MAX( mxl_min,  &
                     &               SQRT( 2._wp * tke_abl( ji, jj, jk, nt_a ) / zbuoy )   ) 
               END DO
            END DO   
            !
            ! Limit mxl
            DO jk = jpkam1,1,-1   
               DO ji = 1, jpi 
                  zlup(ji,jk) = MIN( zlup(ji,jk+1) + (ghw_abl(jk+1)-ghw_abl(jk)) , mxl_abl(ji, jj, jk) )
               END DO
            END DO   
            !
            DO jk = 2, jpka
               DO ji = 1, jpi  
                  zldw(ji,jk) = MIN( zldw(ji,jk-1) + (ghw_abl(jk)-ghw_abl(jk-1)) , mxl_abl(ji, jj, jk) )   
               END DO
            END DO 
            !
            DO jk = 1, jpka
               DO ji = 1, jpi
                  mxl_abl( ji, jj, jk ) = SQRT( zldw( ji, jk ) * zlup( ji, jk ) )   
               END DO
            END DO                         
            !
         END DO
#   undef zlup      
#   undef zldw  
         !
         !
      CASE ( 1 )             ! length-scale computed as the distance to the PBL height
         DO jj = 1,jpj      ! outer loop
            !
            DO ji = 1, jpi   ! vector opt.
               zcff = 1._wp / pblh( ji, jj )     ! inverse of hbl              
               DO jk = 1, jpka              
                  zsig  = MIN( zcff * ghw_abl( jk ), 1. )    
                  zcff1 = pblh( ji, jj )                 
                  mxl_abl( ji, jj, jk ) =  mxl_min                           &
                     &  +  zsig         * (  amx1*zcff1 + bmx1*mxl_min ) &
                     &  +  zsig *  zsig * (  amx2*zcff1 + bmx2*mxl_min ) &
                     &  +  zsig**3      * (  amx3*zcff1 + bmx3*mxl_min ) &
                     &  +  zsig**4      * (  amx4*zcff1 + bmx4*mxl_min ) &
                     &  +  zsig**5      * (  amx5*zcff1 + bmx5*mxl_min )
               END DO
            END DO 
            !            
         END DO            
         !
      CASE ( 2 )           ! Bougeault & Lacarrere 89 length-scale
         !
#   define zlup zRH      
#   define zldw zFC           
! zCF is used for matrix inversion
!             
       DO jj = 1, jpj      ! outer loop
         
         DO ji = 1, jpi           
            zlup( ji,    1 ) = mxl_min
            zldw( ji,    1 ) = mxl_min               
            zlup( ji, jpka ) = mxl_min
            zldw( ji, jpka ) = mxl_min                
         END DO            
         
         DO jk = 2,jpka-1
            DO ji = 1, jpi
               zlup(ji,jk) = ghw_abl(jpka) - ghw_abl(jk)
               zldw(ji,jk) = ghw_abl(jk  ) - ghw_abl( 1)        
            END DO
         END DO         
         !!
         !! BL89 search for lup
         !! ----------------------------------------------------------           
         DO jk=2,jpka-1 
            !
            DO ji = 1, jpi
               zCF(ji,1:jpka) = 0._wp
               zCF(ji,  jk  ) = - tke_abl( ji, jj, jk, nt_a )
               ln_foundl(ji ) = .false.
            END DO   
            !           
            DO jkup=jk+1,jpka-1
               DO ji = 1, jpi
                  zCF (ji,jkup) = zCF (ji,jkup-1) + 0.5_wp * e3t_abl(jkup) * &
                     &               ( zbn2(ji,jj,jkup  )*(ghw_abl(jkup  )-ghw_abl(jk)) &
                     &               + zbn2(ji,jj,jkup-1)*(ghw_abl(jkup-1)-ghw_abl(jk)) )
                  IF( zCF (ji,jkup) * zCF (ji,jkup-1) .le. 0._wp .and. .not. ln_foundl(ji) ) THEN
                     zcff2 = ghw_abl(jkup  ) - ghw_abl(jk)
                     zcff1 = ghw_abl(jkup-1) - ghw_abl(jk)                   
                     zcff  = ( zcff1 * zCF(ji,jkup) - zcff2 * zCF(ji,jkup-1) ) /   &
                        &    (         zCF(ji,jkup) -         zCF(ji,jkup-1) ) 
                     zlup(ji,jk) = zcff                
                     ln_foundl(ji) = .true.
                  END IF
               END DO
            END DO
            !
         END DO   
         !!
         !! BL89 search for ldwn
         !! ----------------------------------------------------------          
         DO jk=2,jpka-1         
            !
            DO ji = 1, jpi
               zCF(ji,1:jpka) = 0._wp
               zCF(ji,  jk  ) = - tke_abl( ji, jj, jk, nt_a )
               ln_foundl(ji ) = .false.
            END DO  
            !   
            DO jkdwn=jk-1,1,-1
               DO ji = 1, jpi             
                  zCF (ji,jkdwn) = zCF (ji,jkdwn+1) + 0.5_wp * e3t_abl(jkdwn+1)  &
                     &               * ( zbn2(ji,jj,jkdwn+1)*(ghw_abl(jk)-ghw_abl(jkdwn+1)) &
                                       + zbn2(ji,jj,jkdwn  )*(ghw_abl(jk)-ghw_abl(jkdwn  )) )   
                  IF(zCF (ji,jkdwn) * zCF (ji,jkdwn+1) .le. 0._wp  .and. .not. ln_foundl(ji) ) THEN 
                     zcff2 = ghw_abl(jk) - ghw_abl(jkdwn+1)
                     zcff1 = ghw_abl(jk) - ghw_abl(jkdwn  )              
                     zcff  = ( zcff1 * zCF(ji,jkdwn+1) - zcff2 * zCF(ji,jkdwn) ) /   &
                        &    (         zCF(ji,jkdwn+1) -         zCF(ji,jkdwn) ) 
                     zldw(ji,jk) = zcff 
                     ln_foundl(ji) = .true.               
                  END IF                                   
               END DO
            END DO
            !      
         END DO

         DO jk = 1, jpka
            DO ji = 1, jpi  
               mxl_abl( ji, jj, jk ) = MAX( SQRT( zldw( ji, jk ) * zlup( ji, jk ) ), mxl_min )  
            END DO
         END DO  

      END DO
#   undef zlup      
#   undef zldw           
         !
      END SELECT      
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                            !  Finalize the computation of turbulent visc./diff.
      !                            !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
      !-------------
      DO jj = 1, jpj     ! outer loop
      !-------------
         DO jk = 1, jpka   
            DO ji = 1, jpi  ! vector opt.
               zcff              = MAX( rn_phimax, rn_Ric * mxl_abl( ji, jj, jk ) * mxl_abl( ji, jj, jk )  &
                  &                                       * zbn2(ji, jj, jk) / tke_abl( ji, jj, jk, nt_a ) ) 
               zcff2             =  1. / ( 1. + zcff )   !<-- phi_z(z)
               zcff              = mxl_abl( ji, jj, jk ) * SQRT( tke_abl( ji, jj, jk, nt_a ) )
!!FL: MAX function probably useless because of the definition of mxl_min             
               Avm_abl( ji, jj, jk ) = MAX( rn_Cm * zcff         , avm_bak   )
               Avt_abl( ji, jj, jk ) = MAX( rn_Ct * zcff * zcff2 , avt_bak   )                              
            END DO
         END DO
      !-------------         
      END DO      
      !-------------

!---------------------------------------------------------------------------------------------------
   END SUBROUTINE abl_zdf_tke
!===================================================================================================


!===================================================================================================
   SUBROUTINE smooth_pblh( pvar2d, msk )
!---------------------------------------------------------------------------------------------------

      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE smooth_pblh  ***
      !!
      !! ** Purpose :   2D Hanning filter on atmospheric PBL height
      !!
      !! ---------------------------------------------------------------------
	  REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: msk   
	  REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: pvar2d
      INTEGER                                     :: ji,jj
	  REAL(wp)                                    :: smth_a, smth_b
	  REAL(wp), DIMENSION(jpi,jpj)                :: zdX,zdY,zFX,zFY
	  REAL(wp)                                    :: zumsk,zvmsk
      !!
      !!=========================================================
      !!
      !! Hanning filter
      smth_a = 1._wp / 8._wp
      smth_b = 1._wp / 4._wp
      !
      DO_2D_11_10
         zumsk = msk(ji,jj) * msk(ji+1,jj)
         zdX ( ji, jj ) = ( pvar2d( ji+1,jj ) - pvar2d( ji  ,jj ) ) * zumsk
      END_2D
      
	  DO_2D_10_11
         zvmsk = msk(ji,jj) * msk(ji,jj+1)
         zdY ( ji, jj ) = ( pvar2d( ji, jj+1 ) - pvar2d( ji  ,jj ) ) * zvmsk
	  END_2D
      
	  DO_2D_10_00
         zFY ( ji, jj  ) =   zdY ( ji, jj   )                        &
            & +  smth_a*  ( (zdX ( ji, jj+1 ) - zdX( ji-1, jj+1 ))   &
            &            -  (zdX ( ji, jj   ) - zdX( ji-1, jj   ))  )
	  END_2D

      DO_2D_00_10
         zFX( ji, jj  ) =    zdX( ji, jj   )                         &
           &    + smth_a*(  (zdY( ji+1, jj ) - zdY( ji+1, jj-1))     &
           &             -  (zdY( ji  , jj ) - zdY( ji  , jj-1)) )
      END_2D

      DO_2D_00_00
         pvar2d( ji  ,jj ) = pvar2d( ji  ,jj )              &
  &         + msk(ji,jj) * smth_b * (                       &
  &                  zFX( ji, jj ) - zFX( ji-1, jj )        &
  &                 +zFY( ji, jj ) - zFY( ji, jj-1 )  )
      END_2D
	  !!
!---------------------------------------------------------------------------------------------------
   END SUBROUTINE smooth_pblh
!===================================================================================================

!!======================================================================
END MODULE ablmod
