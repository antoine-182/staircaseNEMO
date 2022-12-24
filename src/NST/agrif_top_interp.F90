MODULE agrif_top_interp
   !!======================================================================
   !!                   ***  MODULE  agrif_top_interp  ***
   !! AGRIF: interpolation package for TOP
   !!======================================================================
   !! History :  2.0  !  ??? 
   !!----------------------------------------------------------------------
#if defined key_agrif && defined key_top
   !!----------------------------------------------------------------------
   !!   'key_agrif'                                              AGRIF zoom
   !!   'key_top'                                           on-line tracers
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce
   USE dom_oce      
   USE agrif_oce
   USE agrif_top_sponge
   USE par_trc
   USE trc
   USE vremap
   !
   USE lib_mpp     ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC Agrif_trc, interptrn

  !!----------------------------------------------------------------------
   !! NEMO/NST 4.0 , NEMO Consortium (2018)
   !! $Id: agrif_top_interp.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE Agrif_trc
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE Agrif_trc  ***
      !!----------------------------------------------------------------------
      !
      IF( Agrif_Root() )   RETURN
      !
      Agrif_SpecialValue    = 0._wp
      Agrif_UseSpecialValue = .TRUE.
      !
      CALL Agrif_Bc_variable( trn_id, procname=interptrn )
      Agrif_UseSpecialValue = .FALSE.
      !
   END SUBROUTINE Agrif_trc

   SUBROUTINE interptrn( ptab, i1, i2, j1, j2, k1, k2, n1, n2, before )
      !!----------------------------------------------------------------------
      !!                  *** ROUTINE interptrn ***
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) ::   ptab
      INTEGER                                     , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2, n1, n2
      LOGICAL                                     , INTENT(in   ) ::   before
      !
      INTEGER  ::   ji, jj, jk, jn, ibdy, jbdy   ! dummy loop indices
      INTEGER  ::   imin, imax, jmin, jmax, N_in, N_out
      REAL(wp) ::   zrho, z1, z2, z3, z4, z5, z6, z7

      ! vertical interpolation:
      REAL(wp), DIMENSION(i1:i2,j1:j2,1:jpk,1:jptra) :: ptab_child
      REAL(wp), DIMENSION(k1:k2,1:jptra) :: tabin
      REAL(wp), DIMENSION(k1:k2) :: h_in
      REAL(wp), DIMENSION(1:jpk) :: h_out
      !!----------------------------------------------------------------------

      IF( before ) THEN         
         DO jn = 1,jptra
            DO jk=k1,k2
               DO jj=j1,j2
                 DO ji=i1,i2
                       ptab(ji,jj,jk,jn) = tr(ji,jj,jk,jn,Kmm_a)
                 END DO
              END DO
           END DO
        END DO

# if defined key_vertical
        DO jk=k1,k2
           DO jj=j1,j2
              DO ji=i1,i2
                 ptab(ji,jj,jk,jptra+1) = tmask(ji,jj,jk) * e3t(ji,jj,jk,Kmm_a) 
              END DO
           END DO
        END DO
# endif
      ELSE 

# if defined key_vertical
         DO jj=j1,j2
            DO ji=i1,i2
               ptab_child(ji,jj,:) = 0._wp
               N_in = 0
               DO jk=k1,k2 !k2 = jpk of parent grid
                  IF (ptab(ji,jj,jk,n2) == 0) EXIT
                  N_in = N_in + 1
                  tabin(jk,:) = ptab(ji,jj,jk,n1:n2-1)
                  h_in(N_in) = ptab(ji,jj,jk,n2)
               END DO
               N_out = 0
               DO jk=1,jpk ! jpk of child grid
                  IF (tmask(ji,jj,jk) == 0) EXIT 
                  N_out = N_out + 1
                  h_out(jk) = e3t(ji,jj,jk,Krhs_a)
               ENDDO
               IF (N_in > 0) THEN
                  CALL reconstructandremap(tabin(1:N_in,1:jptra),h_in,ptab_child(ji,jj,1:N_out,1:jptra),h_out,N_in,N_out,jptra)
               ENDIF
            ENDDO
         ENDDO
# else
         ptab_child(i1:i2,j1:j2,1:jpk,1:jptra) = ptab(i1:i2,j1:j2,1:jpk,1:jptra)
# endif
         !
         DO jn=1, jptra
            tr(i1:i2,j1:j2,1:jpk,jn,Krhs_a)=ptab_child(i1:i2,j1:j2,1:jpk,jn)*tmask(i1:i2,j1:j2,1:jpk) 
         END DO

      ENDIF
      !
   END SUBROUTINE interptrn

#else
   !!----------------------------------------------------------------------
   !!   Empty module                                           no TOP AGRIF
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE Agrif_TOP_Interp_empty
      !!---------------------------------------------
      !!   *** ROUTINE agrif_Top_Interp_empty ***
      !!---------------------------------------------
      WRITE(*,*)  'agrif_top_interp : You should not have seen this print! error?'
   END SUBROUTINE Agrif_TOP_Interp_empty
#endif

   !!======================================================================
END MODULE agrif_top_interp
