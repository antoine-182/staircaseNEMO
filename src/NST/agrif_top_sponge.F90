#define SPONGE_TOP

MODULE agrif_top_sponge
   !!======================================================================
   !!                ***  MODULE agrif_top_sponge  ***
   !! AGRIF :   sponge layer pakage for passive tracers (TOP)
   !!======================================================================
   !! History :  2.0  ! 2006-08  (R. Benshila, L. Debreu)  Original code
   !!----------------------------------------------------------------------
#if defined key_agrif && defined key_top
   !!----------------------------------------------------------------------
   !!   Agrif_Sponge_trc : 
   !!   interptrn_sponge :  
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_trc
   USE oce
   USE trc
   USE dom_oce
   USE agrif_oce
   USE agrif_oce_sponge
   USE vremap
   !
   USE in_out_manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_Sponge_trc, interptrn_sponge

   !!----------------------------------------------------------------------
   !! NEMO/NST 4.0 , NEMO Consortium (2018)
   !! $Id: agrif_top_sponge.F90 12489 2020-02-28 15:55:11Z davestorkey $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE Agrif_Sponge_trc
      !!----------------------------------------------------------------------
      !!                   *** ROUTINE Agrif_Sponge_Trc ***
      !!----------------------------------------------------------------------
      REAL(wp) ::   zcoef   ! local scalar
      !!----------------------------------------------------------------------
      !
#if defined SPONGE_TOP
!! Assume persistence 
      zcoef = REAL(Agrif_rhot()-1,wp)/REAL(Agrif_rhot())
      CALL Agrif_sponge
      Agrif_SpecialValue    = 0._wp
      Agrif_UseSpecialValue = .TRUE.
      tabspongedone_trn     = .FALSE.
      CALL Agrif_Bc_Variable( trn_sponge_id, calledweight=zcoef, procname=interptrn_sponge )
      Agrif_UseSpecialValue = .FALSE.
#endif
      !
   END SUBROUTINE Agrif_Sponge_Trc


   SUBROUTINE interptrn_sponge( tabres, i1, i2, j1, j2, k1, k2, n1, n2, before )
      !!----------------------------------------------------------------------
      !!                   *** ROUTINE interptrn_sponge ***
      !!----------------------------------------------------------------------
      INTEGER                                     , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2, n1, n2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) ::   tabres
      LOGICAL                                     , INTENT(in   ) ::   before
      !
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      REAL(wp) ::   zabe1, zabe2, ztrelax
      REAL(wp), DIMENSION(i1:i2,j1:j2)               ::   ztu, ztv
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,1:jptra) ::   trbdiff
      ! vertical interpolation:
      REAL(wp), DIMENSION(i1:i2,j1:j2,jpk,1:jptra) ::tabres_child
      REAL(wp), DIMENSION(k1:k2,1:jptra) :: tabin
      REAL(wp), DIMENSION(k1:k2) :: h_in
      REAL(wp), DIMENSION(1:jpk) :: h_out
      INTEGER :: N_in, N_out
      REAL(wp) :: h_diff
      !!----------------------------------------------------------------------
      !
      IF( before ) THEN
         DO jn = 1, jptra
            DO jk=k1,k2
               DO jj=j1,j2
                  DO ji=i1,i2
                     tabres(ji,jj,jk,jn) = tr(ji,jj,jk,jn,Kbb_a)
                  END DO
               END DO
            END DO
         END DO

# if defined key_vertical
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  tabres(ji,jj,jk,jpts+1) = tmask(ji,jj,jk) * e3t(ji,jj,jk,Kbb_a) 
               END DO
            END DO
         END DO
# endif
      ELSE      
# if defined key_vertical
         tabres_child(:,:,:,:) = 0.
         DO jj=j1,j2
            DO ji=i1,i2
               N_in = 0
               DO jk=k1,k2 !k2 = jpk of parent grid
                  IF (tabres(ji,jj,jk,n2) == 0) EXIT
                  N_in = N_in + 1
                  tabin(jk,:) = tabres(ji,jj,jk,n1:n2-1)
                  h_in(N_in) = tabres(ji,jj,jk,n2)
               END DO
               N_out = 0
               DO jk=1,jpk ! jpk of child grid
                  IF (tmask(ji,jj,jk) == 0) EXIT 
                  N_out = N_out + 1
                  h_out(jk) = e3t(ji,jj,jk,Kbb_a) !Child grid scale factors. Could multiply by e1e2t here instead of division above
               ENDDO
               IF (N_in > 0) THEN
                  CALL reconstructandremap(tabin(1:N_in,1:jptra),h_in,tabres_child(ji,jj,1:N_out,1:jptra),h_out,N_in,N_out,jptra)
               ENDIF
            ENDDO
         ENDDO
# endif

         DO jj=j1,j2
            DO ji=i1,i2
               DO jk=1,jpkm1
# if defined key_vertical
                  trbdiff(ji,jj,jk,1:jptra) = tr(ji,jj,jk,1:jptra,Kbb_a) - tabres_child(ji,jj,jk,1:jptra)
# else
                  trbdiff(ji,jj,jk,1:jptra) = tr(ji,jj,jk,1:jptra,Kbb_a) - tabres(ji,jj,jk,1:jptra)
# endif
               ENDDO
            ENDDO
         ENDDO

         !* set relaxation time scale
         IF( l_1st_euler .AND. lk_agrif_fstep ) THEN   ;   ztrelax =   rn_trelax_tra  / (        rn_Dt )
         ELSE                                          ;   ztrelax =   rn_trelax_tra  / (2._wp * rn_Dt )
         ENDIF

         DO jn = 1, jptra
            DO jk = 1, jpkm1
               DO jj = j1,j2-1
                  DO ji = i1,i2-1
                     zabe1 = rn_sponge_tra * fspu(ji,jj) * e2_e1u(ji,jj) * e3u(ji,jj,jk,Kmm_a) * umask(ji,jj,jk)
                     zabe2 = rn_sponge_tra * fspv(ji,jj) * e1_e2v(ji,jj) * e3v(ji,jj,jk,Kmm_a) * vmask(ji,jj,jk)
                     ztu(ji,jj) = zabe1 * ( trbdiff(ji+1,jj  ,jk,jn) - trbdiff(ji,jj,jk,jn) )
                     ztv(ji,jj) = zabe2 * ( trbdiff(ji  ,jj+1,jk,jn) - trbdiff(ji,jj,jk,jn) )
                  END DO
               END DO
               !
               DO jj = j1+1,j2-1
                  DO ji = i1+1,i2-1
                     IF( .NOT. tabspongedone_trn(ji,jj) ) THEN 
                        tr(ji,jj,jk,jn,Krhs_a) = tr(ji,jj,jk,jn,Krhs_a) + (  ztu(ji,jj) - ztu(ji-1,jj  )     &
                           &                                   + ztv(ji,jj) - ztv(ji  ,jj-1)  )  &
                           &                                * r1_e1e2t(ji,jj) / e3t(ji,jj,jk,Kmm_a)  &
                           &                                - ztrelax * fspt(ji,jj) * trbdiff(ji,jj,jk,jn)
                     ENDIF
                  END DO
               END DO
            END DO
            !
         END DO
         !
         tabspongedone_trn(i1+1:i2-1,j1+1:j2-1) = .TRUE.
      ENDIF
      !                 
   END SUBROUTINE interptrn_sponge

#else
   !!----------------------------------------------------------------------
   !!   Empty module                                           no TOP AGRIF
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE agrif_top_sponge_empty
      WRITE(*,*)  'agrif_top_sponge : You should not have seen this print! error?'
   END SUBROUTINE agrif_top_sponge_empty
#endif

   !!======================================================================
END MODULE agrif_top_sponge
