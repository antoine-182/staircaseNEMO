#define SPONGE && define SPONGE_TOP

MODULE agrif_oce_sponge
   !!======================================================================
   !!                   ***  MODULE  agrif_oce_interp  ***
   !! AGRIF: sponge package for the ocean dynamics (OPA)
   !!======================================================================
   !! History :  2.0  !  2002-06  (XXX)  Original cade
   !!             -   !  2005-11  (XXX) 
   !!            3.2  !  2009-04  (R. Benshila) 
   !!            3.6  !  2014-09  (R. Benshila) 
   !!----------------------------------------------------------------------
#if defined key_agrif
   !!----------------------------------------------------------------------
   !!   'key_agrif'                                              AGRIF zoom
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce
   USE dom_oce
   !
   USE in_out_manager
   USE agrif_oce
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE iom
   USE vremap

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_Sponge, Agrif_Sponge_Tra, Agrif_Sponge_Dyn
   PUBLIC interptsn_sponge, interpun_sponge, interpvn_sponge

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/NST 4.0 , NEMO Consortium (2018)
   !! $Id: agrif_oce_sponge.F90 12489 2020-02-28 15:55:11Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE Agrif_Sponge_Tra
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_Sponge_Tra ***
      !!----------------------------------------------------------------------
      REAL(wp) ::   zcoef   ! local scalar
      
      !!----------------------------------------------------------------------
      !
#if defined SPONGE
      !! Assume persistence:
      zcoef = REAL(Agrif_rhot()-1,wp)/REAL(Agrif_rhot())

      CALL Agrif_Sponge
      Agrif_SpecialValue    = 0._wp
      Agrif_UseSpecialValue = .TRUE.
      tabspongedone_tsn     = .FALSE.
      !
      CALL Agrif_Bc_Variable( tsn_sponge_id, calledweight=zcoef, procname=interptsn_sponge )
      !
      Agrif_UseSpecialValue = .FALSE.
#endif
      !
      CALL iom_put( 'agrif_spu', fspu(:,:))
      CALL iom_put( 'agrif_spv', fspv(:,:))
      !
   END SUBROUTINE Agrif_Sponge_Tra


   SUBROUTINE Agrif_Sponge_dyn
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_Sponge_dyn ***
      !!----------------------------------------------------------------------
      REAL(wp) ::   zcoef   ! local scalar
      !!----------------------------------------------------------------------
      !
#if defined SPONGE
      zcoef = REAL(Agrif_rhot()-1,wp)/REAL(Agrif_rhot())

      Agrif_SpecialValue=0.
      Agrif_UseSpecialValue = ln_spc_dyn
      !
      tabspongedone_u = .FALSE.
      tabspongedone_v = .FALSE.         
      CALL Agrif_Bc_Variable( un_sponge_id, calledweight=zcoef, procname=interpun_sponge )
      !
      tabspongedone_u = .FALSE.
      tabspongedone_v = .FALSE.
      CALL Agrif_Bc_Variable( vn_sponge_id, calledweight=zcoef, procname=interpvn_sponge )
      !
      Agrif_UseSpecialValue = .FALSE.
#endif
      !
      CALL iom_put( 'agrif_spt', fspt(:,:))
      CALL iom_put( 'agrif_spf', fspf(:,:))
      !
   END SUBROUTINE Agrif_Sponge_dyn


   SUBROUTINE Agrif_Sponge
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE  Agrif_Sponge ***
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, ind1, ind2
      INTEGER  ::   ispongearea, jspongearea
      REAL(wp) ::   z1_ispongearea, z1_jspongearea
      REAL(wp), DIMENSION(jpi,jpj) :: ztabramp
      REAL(wp), DIMENSION(jpjmax)  :: zmskwest,  zmskeast
      REAL(wp), DIMENSION(jpimax)  :: zmsknorth, zmsksouth
      !!----------------------------------------------------------------------
      !
      ! Sponge 1d example with:
      !      iraf = 3 ; nbghost = 3 ; nn_sponge_len = 2
      !                        
      !coarse :     U     T     U     T     U     T     U
      !|            |           |           |           |
      !fine :     t u t u t u t u t u t u t u t u t u t u t
      !sponge val:0   0   0   1  5/6 4/6 3/6 2/6 1/6  0   0
      !           |   ghost     | <-- sponge area  -- > |
      !           |   points    |                       |
      !                         |--> dynamical interface

#if defined SPONGE || defined SPONGE_TOP
      IF (( .NOT. spongedoneT ).OR.( .NOT. spongedoneU )) THEN
         !
         ! Retrieve masks at open boundaries:

         ! --- West --- !
         ztabramp(:,:) = 0._wp
         ind1 = 1+nbghostcells
         DO ji = mi0(ind1), mi1(ind1)                
            ztabramp(ji,:) = ssumask(ji,:)
         END DO
         !
         zmskwest(:) = 0._wp
         zmskwest(1:jpj) = MAXVAL(ztabramp(:,:), dim=1)

         ! --- East --- !
         ztabramp(:,:) = 0._wp
         ind1 = jpiglo - nbghostcells - 1
         DO ji = mi0(ind1), mi1(ind1)                 
            ztabramp(ji,:) = ssumask(ji,:)
         END DO
         !
         zmskeast(:) = 0._wp
         zmskeast(1:jpj) = MAXVAL(ztabramp(:,:), dim=1)

         ! --- South --- !
         ztabramp(:,:) = 0._wp
         ind1 = 1+nbghostcells
         DO jj = mj0(ind1), mj1(ind1)                 
            ztabramp(:,jj) = ssvmask(:,jj)
         END DO
         !
         zmsksouth(:) = 0._wp
         zmsksouth(1:jpi) = MAXVAL(ztabramp(:,:), dim=2)

         ! --- North --- !
         ztabramp(:,:) = 0._wp
         ind1 = jpjglo - nbghostcells - 1
         DO jj = mj0(ind1), mj1(ind1)                 
            ztabramp(:,jj) = ssvmask(:,jj)
         END DO
         !
         zmsknorth(:) = 0._wp
         zmsknorth(1:jpi) = MAXVAL(ztabramp(:,:), dim=2)
         ! JC: SPONGE MASKING TO BE SORTED OUT:
         zmskwest(:)  = 1._wp
         zmskeast(:)  = 1._wp
         zmsknorth(:) = 1._wp
         zmsksouth(:) = 1._wp
#if defined key_mpp_mpi
!         CALL mpp_max( 'AGRIF_sponge', zmskwest(:) , jpjmax )
!         CALL mpp_max( 'AGRIF_Sponge', zmskeast(:) , jpjmax )
!         CALL mpp_max( 'AGRIF_Sponge', zmsksouth(:), jpimax )
!         CALL mpp_max( 'AGRIF_Sponge', zmsknorth(:), jpimax )
#endif

         ! Define ramp from boundaries towards domain interior at T-points
         ! Store it in ztabramp

         ispongearea  = nn_sponge_len * Agrif_irhox()
         z1_ispongearea = 1._wp / REAL( ispongearea )
         jspongearea  = nn_sponge_len * Agrif_irhoy()
         z1_jspongearea = 1._wp / REAL( jspongearea )
         
         ztabramp(:,:) = 0._wp

         ! Trick to remove sponge in 2DV domains:
         IF ( nbcellsx <= 3 ) ispongearea = -1
         IF ( nbcellsy <= 3 ) jspongearea = -1

         ! --- West --- !
         ind1 = 1+nbghostcells
         ind2 = 1+nbghostcells + ispongearea 
         DO ji = mi0(ind1), mi1(ind2)   
            DO jj = 1, jpj               
               ztabramp(ji,jj) = REAL( ind2 - mig(ji) ) * z1_ispongearea * zmskwest(jj)
            END DO
         END DO

         ! ghost cells:
         ind1 = 1
         ind2 = nbghostcells + 1
         DO ji = mi0(ind1), mi1(ind2)   
            DO jj = 1, jpj               
               ztabramp(ji,jj) = zmskwest(jj)
            END DO
         END DO

         ! --- East --- !
         ind1 = jpiglo - nbghostcells - ispongearea
         ind2 = jpiglo - nbghostcells
         DO ji = mi0(ind1), mi1(ind2)
            DO jj = 1, jpj
               ztabramp(ji,jj) = MAX( ztabramp(ji,jj), REAL( mig(ji) - ind1 ) * z1_ispongearea) * zmskeast(jj)
            ENDDO
         END DO

         ! ghost cells:
         ind1 = jpiglo - nbghostcells
         ind2 = jpiglo
         DO ji = mi0(ind1), mi1(ind2)
            DO jj = 1, jpj
               ztabramp(ji,jj) = zmskeast(jj)
            ENDDO
         END DO

         ! --- South --- !
         ind1 = 1+nbghostcells
         ind2 = 1+nbghostcells + jspongearea
         DO jj = mj0(ind1), mj1(ind2) 
            DO ji = 1, jpi
               ztabramp(ji,jj) = MAX( ztabramp(ji,jj), REAL( ind2 - mjg(jj) ) * z1_jspongearea) * zmsksouth(ji)
            END DO
         END DO

         ! ghost cells:
         ind1 = 1
         ind2 = nbghostcells + 1
         DO jj = mj0(ind1), mj1(ind2) 
            DO ji = 1, jpi
               ztabramp(ji,jj) = zmsksouth(ji)
            END DO
         END DO

         ! --- North --- !
         ind1 = jpjglo - nbghostcells - jspongearea
         ind2 = jpjglo - nbghostcells
         DO jj = mj0(ind1), mj1(ind2)
            DO ji = 1, jpi
               ztabramp(ji,jj) = MAX( ztabramp(ji,jj), REAL( mjg(jj) - ind1 ) * z1_jspongearea) * zmsknorth(ji)
            END DO
         END DO

         ! ghost cells:
         ind1 = jpjglo - nbghostcells
         ind2 = jpjglo
         DO jj = mj0(ind1), mj1(ind2)
            DO ji = 1, jpi
               ztabramp(ji,jj) = zmsknorth(ji)
            END DO
         END DO

      ENDIF

      ! Tracers
      IF( .NOT. spongedoneT ) THEN
         fspu(:,:) = 0._wp
         fspv(:,:) = 0._wp
         DO_2D_00_00
            fspu(ji,jj) = 0.5_wp * ( ztabramp(ji,jj) + ztabramp(ji+1,jj  ) ) * ssumask(ji,jj)
            fspv(ji,jj) = 0.5_wp * ( ztabramp(ji,jj) + ztabramp(ji  ,jj+1) ) * ssvmask(ji,jj)
         END_2D
         CALL lbc_lnk( 'agrif_Sponge', fspu, 'U', 1. )   ! Lateral boundary conditions
         CALL lbc_lnk( 'agrif_Sponge', fspv, 'V', 1. )

         spongedoneT = .TRUE.
      ENDIF

      ! Dynamics
      IF( .NOT. spongedoneU ) THEN
         fspt(:,:) = 0._wp
         fspf(:,:) = 0._wp
         DO_2D_00_00
            fspt(ji,jj) = ztabramp(ji,jj) * ssmask(ji,jj)
            fspf(ji,jj) = 0.25_wp * ( ztabramp(ji  ,jj  ) + ztabramp(ji  ,jj+1)   &
                                  &  +ztabramp(ji+1,jj+1) + ztabramp(ji+1,jj  ) ) &
                                  &  * ssvmask(ji,jj) * ssvmask(ji,jj+1)
         END_2D
         CALL lbc_lnk( 'agrif_Sponge', fspt, 'T', 1. )   ! Lateral boundary conditions
         CALL lbc_lnk( 'agrif_Sponge', fspf, 'F', 1. )
         
         spongedoneU = .TRUE.
      ENDIF

#if defined key_vertical
      ! Remove vertical interpolation where not needed:
      DO_2D_00_00
         IF ((fspu(ji-1,jj)==0._wp).AND.(fspu(ji,jj)==0._wp).AND. &
         &   (fspv(ji,jj-1)==0._wp).AND.(fspv(ji,jj)==0._wp)) mbkt_parent(ji,jj) = 0
!
         IF ((fspt(ji+1,jj)==0._wp).AND.(fspt(ji,jj)==0._wp).AND. &
         &   (fspf(ji,jj-1)==0._wp).AND.(fspf(ji,jj)==0._wp)) mbku_parent(ji,jj) = 0
!
         IF ((fspt(ji,jj+1)==0._wp).AND.(fspt(ji,jj)==0._wp).AND. &
         &   (fspf(ji-1,jj)==0._wp).AND.(fspf(ji,jj)==0._wp)) mbkv_parent(ji,jj) = 0
!
         IF ( ssmask(ji,jj) == 0._wp) mbkt_parent(ji,jj) = 0
         IF (ssumask(ji,jj) == 0._wp) mbku_parent(ji,jj) = 0
         IF (ssvmask(ji,jj) == 0._wp) mbkv_parent(ji,jj) = 0
      END_2D
      !
      ztabramp(:,:) = REAL( mbkt_parent(:,:), wp )   ;   CALL lbc_lnk( 'Agrif_Sponge', ztabramp, 'T', 1. )
      mbkt_parent(:,:) = NINT( ztabramp(:,:) )
      ztabramp(:,:) = REAL( mbku_parent(:,:), wp )   ;   CALL lbc_lnk( 'Agrif_Sponge', ztabramp, 'U', 1. )
      mbku_parent(:,:) = NINT( ztabramp(:,:) )
      ztabramp(:,:) = REAL( mbkv_parent(:,:), wp )   ;   CALL lbc_lnk( 'Agrif_Sponge', ztabramp, 'V', 1. )
      mbkv_parent(:,:) = NINT( ztabramp(:,:) )
#endif
      !
#endif
      !
   END SUBROUTINE Agrif_Sponge

   SUBROUTINE interptsn_sponge( tabres, i1, i2, j1, j2, k1, k2, n1, n2, before )
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE interptsn_sponge ***
      !!----------------------------------------------------------------------
      INTEGER                                     , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2, n1, n2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) ::   tabres
      LOGICAL                                     , INTENT(in   ) ::   before
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::   iku, ikv
      REAL(wp) :: ztsa, zabe1, zabe2, zbtr, zhtot, ztrelax
      REAL(wp), DIMENSION(i1:i2,j1:j2,jpk) :: ztu, ztv
      REAL(wp), DIMENSION(i1:i2,j1:j2,jpk,n1:n2) ::tsbdiff
      ! vertical interpolation:
      REAL(wp), DIMENSION(i1:i2,j1:j2,jpk,n1:n2) ::tabres_child
      REAL(wp), DIMENSION(k1:k2,n1:n2-1) :: tabin
      REAL(wp), DIMENSION(k1:k2) :: h_in
      REAL(wp), DIMENSION(1:jpk) :: h_out
      INTEGER :: N_in, N_out
      !!----------------------------------------------------------------------
      !
      IF( before ) THEN
         DO jn = 1, jpts
            DO jk=k1,k2
               DO jj=j1,j2
                  DO ji=i1,i2
                     tabres(ji,jj,jk,jn) = ts(ji,jj,jk,jn,Kbb_a)
                  END DO
               END DO
            END DO
         END DO

# if defined key_vertical
        ! Interpolate thicknesses
        ! Warning: these are masked, hence extrapolated prior interpolation.
        DO jk=k1,k2
           DO jj=j1,j2
              DO ji=i1,i2
                  tabres(ji,jj,jk,jpts+1) = tmask(ji,jj,jk) * e3t(ji,jj,jk,Kbb_a)
              END DO
           END DO
        END DO

        ! Extrapolate thicknesses in partial bottom cells:
        ! Set them to Agrif_SpecialValue (0.). Correct bottom thicknesses are retrieved later on
        IF (ln_zps) THEN
           DO jj=j1,j2
              DO ji=i1,i2
                  jk = mbkt(ji,jj)
                  tabres(ji,jj,jk,jpts+1) = 0._wp
              END DO
           END DO           
        END IF
     
        ! Save ssh at last level:
        IF (.NOT.ln_linssh) THEN
           tabres(i1:i2,j1:j2,k2,jpts+1) = ssh(i1:i2,j1:j2,Kbb_a)*tmask(i1:i2,j1:j2,1) 
        ELSE
           tabres(i1:i2,j1:j2,k2,jpts+1) = 0._wp
        END IF      
# endif

      ELSE   
         !
# if defined key_vertical

         IF (ln_linssh) tabres(i1:i2,j1:j2,k2,n2) = 0._wp

         DO jj=j1,j2
            DO ji=i1,i2
               tabres_child(ji,jj,:,:) = 0._wp 
               N_in = mbkt_parent(ji,jj)
               zhtot = 0._wp
               DO jk=1,N_in !k2 = jpk of parent grid
                  IF (jk==N_in) THEN
                     h_in(jk) = ht0_parent(ji,jj) + tabres(ji,jj,k2,n2) - zhtot
                  ELSE
                     h_in(jk) = tabres(ji,jj,jk,n2)
                  ENDIF
                  zhtot = zhtot + h_in(jk)
                  tabin(jk,:) = tabres(ji,jj,jk,n1:n2-1)
               END DO
               N_out = 0
               DO jk=1,jpk ! jpk of child grid
                  IF (tmask(ji,jj,jk) == 0) EXIT 
                  N_out = N_out + 1
                  h_out(jk) = e3t(ji,jj,jk,Kbb_a) !Child grid scale factors. Could multiply by e1e2t here instead of division above
               ENDDO

               ! Account for small differences in free-surface
               IF ( sum(h_out(1:N_out)) > sum(h_in(1:N_in) )) THEN
                  h_out(1) = h_out(1) - ( sum(h_out(1:N_out))-sum(h_in(1:N_in)) )
               ELSE
                  h_in(1)   = h_in(1)   - (sum(h_in(1:N_in))-sum(h_out(1:N_out)) )
               ENDIF
               IF (N_in*N_out > 0) THEN
                  CALL reconstructandremap(tabin(1:N_in,1:jpts),h_in(1:N_in),tabres_child(ji,jj,1:N_out,1:jpts),h_out(1:N_out),N_in,N_out,jpts)
               ENDIF
            ENDDO
         ENDDO
# endif

         DO jj=j1,j2
            DO ji=i1,i2
               DO jk=1,jpkm1
# if defined key_vertical
                  tsbdiff(ji,jj,jk,1:jpts) = (ts(ji,jj,jk,1:jpts,Kbb_a) - tabres_child(ji,jj,jk,1:jpts)) * tmask(ji,jj,jk)
# else
                  tsbdiff(ji,jj,jk,1:jpts) = (ts(ji,jj,jk,1:jpts,Kbb_a) - tabres(ji,jj,jk,1:jpts))*tmask(ji,jj,jk)
# endif
               ENDDO
            ENDDO
         ENDDO

         !* set relaxation time scale
         IF( l_1st_euler .AND. lk_agrif_fstep ) THEN   ;   ztrelax =   rn_trelax_tra  / (        rn_Dt )
         ELSE                                          ;   ztrelax =   rn_trelax_tra  / (2._wp * rn_Dt )
         ENDIF

         DO jn = 1, jpts            
            DO jk = 1, jpkm1
               ztu(i1:i2,j1:j2,jk) = 0._wp
               DO jj = j1,j2
                  DO ji = i1,i2-1
                     zabe1 = rn_sponge_tra * fspu(ji,jj) * umask(ji,jj,jk) * e2_e1u(ji,jj) * e3u(ji,jj,jk,Kmm_a)
                     ztu(ji,jj,jk) = zabe1 * ( tsbdiff(ji+1,jj  ,jk,jn) - tsbdiff(ji,jj,jk,jn) ) 
                  END DO
               END DO
               ztv(i1:i2,j1:j2,jk) = 0._wp
               DO ji = i1,i2
                  DO jj = j1,j2-1
                     zabe2 = rn_sponge_tra * fspv(ji,jj) * vmask(ji,jj,jk) * e1_e2v(ji,jj) * e3v(ji,jj,jk,Kmm_a)
                     ztv(ji,jj,jk) = zabe2 * ( tsbdiff(ji  ,jj+1,jk,jn) - tsbdiff(ji,jj,jk,jn) )
                  END DO
               END DO
               !
               IF( ln_zps ) THEN      ! set gradient at partial step level
                  DO jj = j1,j2
                     DO ji = i1,i2
                        ! last level
                        iku = mbku(ji,jj)
                        ikv = mbkv(ji,jj)
                        IF( iku == jk )   ztu(ji,jj,jk) = 0._wp
                        IF( ikv == jk )   ztv(ji,jj,jk) = 0._wp
                     END DO
                  END DO
               ENDIF
            END DO
            !
            DO jk = 1, jpkm1
               DO jj = j1+1,j2-1
                  DO ji = i1+1,i2-1
                     IF (.NOT. tabspongedone_tsn(ji,jj)) THEN 
                        zbtr = r1_e1e2t(ji,jj) / e3t(ji,jj,jk,Kmm_a)
                        ! horizontal diffusive trends
                        ztsa = zbtr * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  ) &
                             &  - ztrelax * fspt(ji,jj) * tsbdiff(ji,jj,jk,jn) 
                        ! add it to the general tracer trends
                        ts(ji,jj,jk,jn,Krhs_a) = ts(ji,jj,jk,jn,Krhs_a) + ztsa
                     ENDIF
                  END DO
               END DO
            END DO
            !
         END DO
         !
         tabspongedone_tsn(i1+1:i2-1,j1+1:j2-1) = .TRUE.
         !
      ENDIF
      !
   END SUBROUTINE interptsn_sponge

   SUBROUTINE interpun_sponge(tabres,i1,i2,j1,j2,k1,k2,m1,m2, before)
      !!---------------------------------------------
      !!   *** ROUTINE interpun_sponge ***
      !!---------------------------------------------    
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,m1,m2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,m1:m2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before

      INTEGER :: ji,jj,jk,jmax

      ! sponge parameters 
      REAL(wp) :: ze2u, ze1v, zua, zva, zbtr, zhtot, ztrelax
      REAL(wp), DIMENSION(i1:i2,j1:j2,1:jpk) :: ubdiff
      REAL(wp), DIMENSION(i1:i2,j1:j2,1:jpk) :: rotdiff, hdivdiff
      ! vertical interpolation:
      REAL(wp), DIMENSION(i1:i2,j1:j2,1:jpk) :: tabres_child
      REAL(wp), DIMENSION(k1:k2) :: tabin, h_in
      REAL(wp), DIMENSION(1:jpk) :: h_out
      INTEGER ::N_in, N_out
      !!---------------------------------------------    
      !
      IF( before ) THEN
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj,jk,m1) = uu(ji,jj,jk,Kbb_a)
# if defined key_vertical
                  tabres(ji,jj,jk,m2) = e3u(ji,jj,jk,Kbb_a)*umask(ji,jj,jk)
# endif
               END DO
            END DO
         END DO

# if defined key_vertical
         ! Extrapolate thicknesses in partial bottom cells:
         ! Set them to Agrif_SpecialValue (0.). Correct bottom thicknesses are retrieved later on
         IF (ln_zps) THEN
            DO jj=j1,j2
               DO ji=i1,i2
                  jk = mbku(ji,jj)
                  tabres(ji,jj,jk,m2) = 0._wp
               END DO
            END DO           
         END IF
        ! Save ssh at last level:
        tabres(i1:i2,j1:j2,k2,m2) = 0._wp
        IF (.NOT.ln_linssh) THEN
           ! This vertical sum below should be replaced by the sea-level at U-points (optimization):
           DO jk=1,jpk
              tabres(i1:i2,j1:j2,k2,m2) = tabres(i1:i2,j1:j2,k2,m2) + e3u(i1:i2,j1:j2,jk,Kbb_a) * umask(i1:i2,j1:j2,jk)
           END DO
           tabres(i1:i2,j1:j2,k2,m2) = tabres(i1:i2,j1:j2,k2,m2) - hu_0(i1:i2,j1:j2)
        END IF 
# endif

      ELSE

# if defined key_vertical
         IF (ln_linssh) tabres(i1:i2,j1:j2,k2,m2) = 0._wp

         DO jj=j1,j2
            DO ji=i1,i2
               tabres_child(ji,jj,:) = 0._wp
               N_in = mbku_parent(ji,jj)
               zhtot = 0._wp
               DO jk=1,N_in
                  IF (jk==N_in) THEN
                     h_in(jk) = hu0_parent(ji,jj) + tabres(ji,jj,k2,m2) - zhtot
                  ELSE
                     h_in(jk) = tabres(ji,jj,jk,m2)
                  ENDIF
                  zhtot = zhtot + h_in(jk)
                  tabin(jk) = tabres(ji,jj,jk,m1)
               ENDDO
               !         
               N_out = 0
               DO jk=1,jpk
                  IF (umask(ji,jj,jk) == 0) EXIT
                  N_out = N_out + 1
                  h_out(N_out) = e3u(ji,jj,jk,Kbb_a)
               ENDDO

               ! Account for small differences in free-surface
               IF ( sum(h_out(1:N_out)) > sum(h_in(1:N_in) )) THEN
                  h_out(1) = h_out(1) - ( sum(h_out(1:N_out))-sum(h_in(1:N_in)) )
               ELSE
                  h_in(1)   = h_in(1)   - (sum(h_in(1:N_in))-sum(h_out(1:N_out)) )
               ENDIF
                  
               IF (N_in * N_out > 0) THEN
                  CALL reconstructandremap(tabin(1:N_in),h_in(1:N_in),tabres_child(ji,jj,1:N_out),h_out(1:N_out),N_in,N_out,1)
               ENDIF 
            ENDDO
         ENDDO

         ubdiff(i1:i2,j1:j2,:) = (uu(i1:i2,j1:j2,:,Kbb_a) - tabres_child(i1:i2,j1:j2,:))*umask(i1:i2,j1:j2,:)
#else
         ubdiff(i1:i2,j1:j2,:) = (uu(i1:i2,j1:j2,:,Kbb_a) - tabres(i1:i2,j1:j2,:,1))*umask(i1:i2,j1:j2,:)
#endif
         !* set relaxation time scale
         IF( l_1st_euler .AND. lk_agrif_fstep ) THEN   ;   ztrelax =   rn_trelax_dyn  / (        rn_Dt )
         ELSE                                          ;   ztrelax =   rn_trelax_dyn  / (2._wp * rn_Dt )
         ENDIF
         !
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============

            !                                             ! --------
            ! Horizontal divergence                       !   div
            !                                             ! --------
            DO jj = j1,j2
               DO ji = i1+1,i2   ! vector opt.
                  zbtr = r1_e1e2t(ji,jj) / e3t(ji,jj,jk,Kbb_a) * rn_sponge_dyn * fspt(ji,jj)
                  hdivdiff(ji,jj,jk) = (  e2u(ji  ,jj)*e3u(ji  ,jj,jk,Kbb_a) * ubdiff(ji  ,jj,jk) &
                                     &   -e2u(ji-1,jj)*e3u(ji-1,jj,jk,Kbb_a) * ubdiff(ji-1,jj,jk) ) * zbtr
               END DO
            END DO

            DO jj = j1,j2-1
               DO ji = i1,i2   ! vector opt.
                  zbtr = r1_e1e2f(ji,jj) * e3f(ji,jj,jk) * rn_sponge_dyn * fspf(ji,jj)
                  rotdiff(ji,jj,jk) = ( -e1u(ji,jj+1) * ubdiff(ji,jj+1,jk)   &
                                    &   +e1u(ji,jj  ) * ubdiff(ji,jj  ,jk) ) * fmask(ji,jj,jk) * zbtr 
               END DO
            END DO
         END DO
         !
         DO jj = j1+1, j2-1
            DO ji = i1+1, i2-1   ! vector opt.

               IF (.NOT. tabspongedone_u(ji,jj)) THEN
                  DO jk = 1, jpkm1                                 ! Horizontal slab
                     ze2u = rotdiff (ji,jj,jk)
                     ze1v = hdivdiff(ji,jj,jk)
                     ! horizontal diffusive trends
                     zua = - ( ze2u - rotdiff (ji,jj-1,jk) ) / ( e2u(ji,jj) * e3u(ji,jj,jk,Kmm_a) )   &
                         & + ( hdivdiff(ji+1,jj,jk) - ze1v ) * r1_e1u(ji,jj) & 
                         & - ztrelax  * fspu(ji,jj) * ubdiff(ji,jj,jk)

                     ! add it to the general momentum trends
                     uu(ji,jj,jk,Krhs_a) = uu(ji,jj,jk,Krhs_a) + zua                                 
                  END DO
               ENDIF

            END DO
         END DO

         tabspongedone_u(i1+1:i2-1,j1+1:j2-1) = .TRUE.

         jmax = j2-1
         IF ((nbondj == 1).OR.(nbondj == 2)) jmax = MIN(jmax,nlcj-nbghostcells-2)   ! North

         DO jj = j1+1, jmax
            DO ji = i1+1, i2   ! vector opt.

               IF (.NOT. tabspongedone_v(ji,jj)) THEN
                  DO jk = 1, jpkm1                                 ! Horizontal slab
                     ze2u = rotdiff (ji,jj,jk)
                     ze1v = hdivdiff(ji,jj,jk)

                     ! horizontal diffusive trends
                     zva = + ( ze2u - rotdiff (ji-1,jj,jk) ) / ( e1v(ji,jj) * e3v(ji,jj,jk,Kmm_a) )   &
                           + ( hdivdiff(ji,jj+1,jk) - ze1v ) * r1_e2v(ji,jj)

                     ! add it to the general momentum trends
                     vv(ji,jj,jk,Krhs_a) = vv(ji,jj,jk,Krhs_a) + zva
                  END DO
               ENDIF
               !
            END DO
         END DO
         !
         tabspongedone_v(i1+1:i2,j1+1:jmax) = .TRUE.
         !
      ENDIF
      !
   END SUBROUTINE interpun_sponge

   SUBROUTINE interpvn_sponge(tabres,i1,i2,j1,j2,k1,k2,m1,m2, before,nb,ndir)
      !!---------------------------------------------
      !!   *** ROUTINE interpvn_sponge ***
      !!--------------------------------------------- 
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,m1,m2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,m1:m2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      INTEGER, INTENT(in) :: nb , ndir
      !
      INTEGER  ::   ji, jj, jk, imax
      REAL(wp) ::   ze2u, ze1v, zua, zva, zbtr, zhtot, ztrelax
      REAL(wp), DIMENSION(i1:i2,j1:j2,1:jpk) :: vbdiff
      REAL(wp), DIMENSION(i1:i2,j1:j2,1:jpk) :: rotdiff, hdivdiff
      ! vertical interpolation:
      REAL(wp), DIMENSION(i1:i2,j1:j2,1:jpk) :: tabres_child
      REAL(wp), DIMENSION(k1:k2) :: tabin, h_in
      REAL(wp), DIMENSION(1:jpk) :: h_out
      INTEGER :: N_in, N_out
      !!--------------------------------------------- 
      
      IF( before ) THEN 
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj,jk,m1) = vv(ji,jj,jk,Kbb_a)
# if defined key_vertical
                  tabres(ji,jj,jk,m2) = vmask(ji,jj,jk) * e3v(ji,jj,jk,Kbb_a)
# endif
               END DO
            END DO
         END DO

# if defined key_vertical
         ! Extrapolate thicknesses in partial bottom cells:
         ! Set them to Agrif_SpecialValue (0.). Correct bottom thicknesses are retrieved later on
         IF (ln_zps) THEN
            DO jj=j1,j2
               DO ji=i1,i2
                  jk = mbkv(ji,jj)
                  tabres(ji,jj,jk,m2) = 0._wp
               END DO
            END DO           
         END IF
        ! Save ssh at last level:
        tabres(i1:i2,j1:j2,k2,m2) = 0._wp
        IF (.NOT.ln_linssh) THEN
           ! This vertical sum below should be replaced by the sea-level at V-points (optimization):
           DO jk=1,jpk
              tabres(i1:i2,j1:j2,k2,m2) = tabres(i1:i2,j1:j2,k2,m2) + e3v(i1:i2,j1:j2,jk,Kbb_a) * vmask(i1:i2,j1:j2,jk)
           END DO
           tabres(i1:i2,j1:j2,k2,m2) = tabres(i1:i2,j1:j2,k2,m2) - hv_0(i1:i2,j1:j2)
        END IF 
# endif

      ELSE

# if defined key_vertical
         IF (ln_linssh) tabres(i1:i2,j1:j2,k2,m2) = 0._wp
         DO jj=j1,j2
            DO ji=i1,i2
               tabres_child(ji,jj,:) = 0._wp
               N_in = mbkv_parent(ji,jj)
               zhtot = 0._wp
               DO jk=1,N_in
                  IF (jk==N_in) THEN
                     h_in(jk) = hv0_parent(ji,jj) + tabres(ji,jj,k2,m2) - zhtot
                  ELSE
                     h_in(jk) = tabres(ji,jj,jk,m2)
                  ENDIF
                  zhtot = zhtot + h_in(jk)
                  tabin(jk) = tabres(ji,jj,jk,m1)
               ENDDO
               !          
               N_out = 0
               DO jk=1,jpk
                  IF (vmask(ji,jj,jk) == 0) EXIT
                  N_out = N_out + 1
                  h_out(N_out) = e3v(ji,jj,jk,Kbb_a)
               ENDDO

               ! Account for small differences in free-surface
               IF ( sum(h_out(1:N_out)) > sum(h_in(1:N_in) )) THEN
                  h_out(1) = h_out(1) - ( sum(h_out(1:N_out))-sum(h_in(1:N_in)) )
               ELSE
                  h_in(1)   = h_in(1) - (  sum(h_in(1:N_in))-sum(h_out(1:N_out)) )
               ENDIF
         
               IF (N_in * N_out > 0) THEN
                  CALL reconstructandremap(tabin(1:N_in),h_in(1:N_in),tabres_child(ji,jj,1:N_out),h_out(1:N_out),N_in,N_out,1)
               ENDIF
            ENDDO
         ENDDO

         vbdiff(i1:i2,j1:j2,:) = (vv(i1:i2,j1:j2,:,Kbb_a) - tabres_child(i1:i2,j1:j2,:))*vmask(i1:i2,j1:j2,:)  
# else
         vbdiff(i1:i2,j1:j2,:) = (vv(i1:i2,j1:j2,:,Kbb_a) - tabres(i1:i2,j1:j2,:,1))*vmask(i1:i2,j1:j2,:)
# endif
         !* set relaxation time scale
         IF( l_1st_euler .AND. lk_agrif_fstep ) THEN   ;   ztrelax =   rn_trelax_dyn  / (        rn_Dt )
         ELSE                                          ;   ztrelax =   rn_trelax_dyn  / (2._wp * rn_Dt )
         ENDIF
         !
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============

            !                                             ! --------
            ! Horizontal divergence                       !   div
            !                                             ! --------
            DO jj = j1+1,j2
               DO ji = i1,i2   ! vector opt.
                  zbtr = r1_e1e2t(ji,jj) / e3t(ji,jj,jk,Kbb_a) * rn_sponge_dyn * fspt(ji,jj)
                  hdivdiff(ji,jj,jk) = ( e1v(ji,jj  ) * e3v(ji,jj  ,jk,Kbb_a) * vbdiff(ji,jj  ,jk)  &
                                     &  -e1v(ji,jj-1) * e3v(ji,jj-1,jk,Kbb_a) * vbdiff(ji,jj-1,jk)  ) * zbtr
               END DO
            END DO
            DO jj = j1,j2
               DO ji = i1,i2-1   ! vector opt.
                  zbtr = r1_e1e2f(ji,jj) * e3f(ji,jj,jk) * rn_sponge_dyn * fspf(ji,jj)
                  rotdiff(ji,jj,jk) = ( e2v(ji+1,jj) * vbdiff(ji+1,jj,jk) & 
                                    &  -e2v(ji  ,jj) * vbdiff(ji  ,jj,jk)  ) * fmask(ji,jj,jk) * zbtr
               END DO
            END DO
         END DO

         !                                                ! ===============
         !                                                

         imax = i2 - 1
         IF ((nbondi == 1).OR.(nbondi == 2))   imax = MIN(imax,nlci-nbghostcells-2)   ! East

         DO jj = j1+1, j2
            DO ji = i1+1, imax   ! vector opt.
               IF( .NOT. tabspongedone_u(ji,jj) ) THEN
                  DO jk = 1, jpkm1
                     uu(ji,jj,jk,Krhs_a) = uu(ji,jj,jk,Krhs_a)                                                               &
                        & - ( rotdiff (ji  ,jj,jk) - rotdiff (ji,jj-1,jk)) / ( e2u(ji,jj) * e3u(ji,jj,jk,Kmm_a) )  &
                        & + ( hdivdiff(ji+1,jj,jk) - hdivdiff(ji,jj  ,jk)) * r1_e1u(ji,jj)
                  END DO
               ENDIF
            END DO
         END DO
         !
         tabspongedone_u(i1+1:imax,j1+1:j2) = .TRUE.
         !
         DO jj = j1+1, j2-1
            DO ji = i1+1, i2-1   ! vector opt.
               IF( .NOT. tabspongedone_v(ji,jj) ) THEN
                  DO jk = 1, jpkm1
                     vv(ji,jj,jk,Krhs_a) = vv(ji,jj,jk,Krhs_a)                                                                  &
                        &  + ( rotdiff (ji,jj  ,jk) - rotdiff (ji-1,jj,jk) ) / ( e1v(ji,jj) * e3v(ji,jj,jk,Kmm_a) )   &
                        &  + ( hdivdiff(ji,jj+1,jk) - hdivdiff(ji  ,jj,jk) ) * r1_e2v(ji,jj)                      &
                        &  - ztrelax * fspv(ji,jj) * vbdiff(ji,jj,jk)
                  END DO
               ENDIF
            END DO
         END DO
         tabspongedone_v(i1+1:i2-1,j1+1:j2-1) = .TRUE.
      ENDIF
      !
   END SUBROUTINE interpvn_sponge

#else
   !!----------------------------------------------------------------------
   !!   Empty module                                          no AGRIF zoom
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE agrif_oce_sponge_empty
      WRITE(*,*)  'agrif_oce_sponge : You should not have seen this print! error?'
   END SUBROUTINE agrif_oce_sponge_empty
#endif

   !!======================================================================
END MODULE agrif_oce_sponge
