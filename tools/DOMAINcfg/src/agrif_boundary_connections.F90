#if defined key_agrif
subroutine agrif_boundary_connections
use agrif_profiles
use agrif_parameters
implicit none
external connect_e3t_copy, connect_e3t_connect, connect_bottom_level,connect_e3u, connect_e3v
if (agrif_root()) return

call Agrif_connection()

call Agrif_Bc_variable(bottom_level_id, procname = connect_bottom_level)

call Agrif_Bc_variable(e3t_copy_id, procname = connect_e3t_copy)

! Agrif_UseSpecialValue = .TRUE.
! Agrif_SpecialValue = 0.
! call Agrif_Bc_variable(e3t_connect_id, procname = connect_e3t_connect)
! Agrif_UseSpecialValue = .FALSE.


      Agrif_SpecialValue    = 0.
      Agrif_UseSpecialValue = ln_spc_dyn
      !
      CALL Agrif_Bc_variable( e3u_id, procname=connect_e3u )
      CALL Agrif_Bc_variable( e3v_id, procname=connect_e3v )
      !
      Agrif_UseSpecialValue = .FALSE.
      
end subroutine agrif_boundary_connections


    SUBROUTINE connect_e3t_copy( ptab, i1, i2, j1, j2, k1, k2, before, nb,ndir)
    USE dom_oce
    USE domzgr
    USE agrif_parameters
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpsshn  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2
      REAL, DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      INTEGER                         , INTENT(in   ) ::   nb , ndir
      LOGICAL  ::   western_side, eastern_side,northern_side,southern_side
      !
      !!---------------------------------------------------------------------- 
      INTEGER :: ji,jj,jk      
      !
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

      IF( before) THEN
         ptab(i1:i2,j1:j2,k1:k2) = e3t_0(i1:i2,j1:j2,k1:k2)
      ELSE
         e3t_0(i1:i2,j1:j2,k1:k2) = ptab(i1:i2,j1:j2,k1:k2)
      ENDIF
      !
   END SUBROUTINE connect_e3t_copy
   
    SUBROUTINE connect_bottom_level( ptab, i1, i2, j1, j2, before, nb,ndir)
    USE dom_oce
    USE domzgr
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpsshn  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL, DIMENSION(i1:i2,j1:j2)    , INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      INTEGER                         , INTENT(in   ) ::   nb , ndir
      LOGICAL  ::   western_side, eastern_side,northern_side,southern_side
      !
      !!---------------------------------------------------------------------- 
      INTEGER :: ji,jj     
      !
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

      IF( before) THEN
         ptab(i1:i2,j1:j2) = mbkt(i1:i2,j1:j2)*ssmask(i1:i2,j1:j2)
      ELSE
         mbkt(i1:i2,j1:j2) = nint(ptab(i1:i2,j1:j2))
         
         WHERE (mbkt(i1:i2,j1:j2)==0)
           ssmask(i1:i2,j1:j2) = 0.
         END WHERE
           
      ENDIF
      !
   END SUBROUTINE connect_bottom_level
   
    SUBROUTINE connect_e3t_connect( ptab, i1, i2, j1, j2, k1, k2, before, nb,ndir)
    USE dom_oce
    USE domzgr
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpsshn  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2
      REAL, DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      INTEGER                         , INTENT(in   ) ::   nb , ndir
      LOGICAL  ::   western_side, eastern_side,northern_side,southern_side
      !
      !!---------------------------------------------------------------------- 
      INTEGER :: ji,jj,jk    
      REAL,DIMENSION(i1:i2,j1:j2) :: bathy_local      
      !
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

      IF( before) THEN
         do jk=1,jpk
         do jj=j1,j2
         do ji=i1,i2
          if (mbkt(ji,jj)>=jk) then
            ptab(ji,jj,jk) = e3t_0(ji,jj,jk)
          else
            ptab(ji,jj,jk) = 0.
          endif
         enddo
         enddo
         enddo
         
         do jj=j1,j2
         do ji=i1,i2
           ptab(ji,jj,jpk+1) = SUM ( e3t_0(ji,jj, 1:mbkt(ji,jj) ) ) * ssmask(ji,jj)
         enddo
         enddo

      ELSE
         do jj=j1,j2
         do ji=i1,i2
           bathy_local (ji,jj) = SUM ( e3t_0(ji,jj, 1:mbkt(ji,jj) ) ) * ssmask(ji,jj)
           print *,'ji = ',ji,jj,bathy_local(ji,jj),ptab(ji,jj,jpk+1)
         enddo
         enddo
         
         ! DO jk=1,jpk
           ! DO jj=j1,j2
             ! DO ji=i1,i2
                 ! e3t_0(ji,jj,jk) = MAX(ptab(ji,jj,jk),MIN(e3zps_min, e3t_1d(jk)*e3zps_rat ))
                 ! e3t_0(ji,jj,jk) = MIN(e3t_0(ji,jj,jk),e3t_1d(jk))
             ! ENDDO
           ! ENDDO
         ! ENDDO
      ENDIF
      !
   END SUBROUTINE connect_e3t_connect
   
   SUBROUTINE connect_e3u( ptab, i1, i2, j1, j2, k1, k2,before, nb, ndir )
   USE dom_oce
      !!----------------------------------------------------------------------
      !!                  *** ROUTINE interpun ***
      !!---------------------------------------------    
      !!
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      INTEGER, INTENT(in) :: nb , ndir
      !!
      INTEGER :: ji,jj,jk
      REAL(wp) :: zrhoy
      ! vertical interpolation:
      REAL(wp), DIMENSION(k1:k2) :: tabin, h_in
      REAL(wp), DIMENSION(1:jpk) :: h_out
      INTEGER  :: N_in, N_out, iref
      REAL(wp) :: h_diff
      LOGICAL  :: western_side, eastern_side
      !!---------------------------------------------    
      !
      IF (before) THEN 
         DO jk=1,jpk
            DO jj=j1,j2
               DO ji=i1,i2
                 if (min(mbkt(ji,jj),mbkt(ji+1,jj))<jk) then
                  ptab(ji,jj,jk) = 0.
                 else
                  ptab(ji,jj,jk) = e2u(ji,jj) * e3u_0(ji,jj,jk)
                 endif
# if defined key_vertical
                  ptab(ji,jj,jk,2) = (umask(ji,jj,jk) * e2u(ji,jj) * e3u_n(ji,jj,jk))
# endif
               END DO
            END DO
         END DO
      ELSE
         zrhoy = Agrif_rhoy()
# if defined key_vertical
! VERTICAL REFINEMENT BEGIN
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)

         DO ji=i1,i2
            iref = ji
            IF (western_side) iref = MAX(2,ji)
            IF (eastern_side) iref = MIN(nlci-2,ji)
            DO jj=j1,j2
               N_in = 0
               DO jk=k1,k2
                  IF (ptab(ji,jj,jk,2) == 0) EXIT
                  N_in = N_in + 1
                  tabin(jk) = ptab(ji,jj,jk,1)/ptab(ji,jj,jk,2)
                  h_in(N_in) = ptab(ji,jj,jk,2)/(e2u(ji,jj)*zrhoy) 
              ENDDO
         
              IF (N_in == 0) THEN
                 ua(ji,jj,:) = 0._wp
                 CYCLE
              ENDIF
         
              N_out = 0
              DO jk=1,jpk
                 if (umask(iref,jj,jk) == 0) EXIT
                 N_out = N_out + 1
                 h_out(N_out) = e3u_a(iref,jj,jk)
              ENDDO
         
              IF (N_out == 0) THEN
                 ua(ji,jj,:) = 0._wp
                 CYCLE
              ENDIF
         
              IF (N_in * N_out > 0) THEN
                 h_diff = sum(h_out(1:N_out))-sum(h_in(1:N_in))
! Should be able to remove the next IF/ELSEIF statement once scale factors are dealt with properly
                 if (h_diff < -1.e4) then
                    print *,'CHECK YOUR BATHY ...', h_diff, sum(h_out(1:N_out)), sum(h_in(1:N_in))
!                    stop
                 endif
              ENDIF
              call reconstructandremap(tabin(1:N_in),h_in(1:N_in),ua(ji,jj,1:N_out),h_out(1:N_out),N_in,N_out)
            ENDDO
         ENDDO

# else
         DO jk = 1, jpkm1
            DO jj=j1,j2
            do ji=i1,i2
              if (min(mbkt(ji+1,jj),mbkt(ji,jj))<jk) then
                e3u_0(ji,jj,jk)=e3t_1d(jk)
              else
                e3u_0(ji,jj,jk) = MAX(ptab(ji,jj,jk) / ( zrhoy * e2u(ji,jj) ),MIN(e3zps_min,e3t_1d(jk)*e3zps_rat))
              endif
            enddo
               
            END DO
         END DO
# endif

      ENDIF
      ! 
   END SUBROUTINE connect_e3u
   
   SUBROUTINE connect_e3v( ptab, i1, i2, j1, j2, k1, k2, before, nb, ndir )
   USE dom_oce
      !!----------------------------------------------------------------------
      !!                  *** ROUTINE interpvn ***
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      INTEGER, INTENT(in) :: nb , ndir
      !
      INTEGER :: ji,jj,jk
      REAL(wp) :: zrhox
      ! vertical interpolation:
      REAL(wp), DIMENSION(k1:k2) :: tabin, h_in
      REAL(wp), DIMENSION(1:jpk) :: h_out
      INTEGER  :: N_in, N_out, jref
      REAL(wp) :: h_diff
      LOGICAL  :: northern_side,southern_side
      !!---------------------------------------------    
      !      
      IF (before) THEN          
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                 if (min(mbkt(ji,jj),mbkt(ji,jj+1))<jk) then
                  ptab(ji,jj,jk) = 0.
                 else
                  ptab(ji,jj,jk) = (e1v(ji,jj) * e3v_0(ji,jj,jk))
                 endif
# if defined key_vertical
                  ptab(ji,jj,jk,2) = vmask(ji,jj,jk) * e1v(ji,jj) * e3v_n(ji,jj,jk)
# endif
               END DO
            END DO
         END DO
      ELSE       
         zrhox = Agrif_rhox()
# if defined key_vertical

         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)

         DO jj=j1,j2
            jref = jj
            IF (southern_side) jref = MAX(2,jj)
            IF (northern_side) jref = MIN(nlcj-2,jj)
            DO ji=i1,i2
               N_in = 0
               DO jk=k1,k2
                  if (ptab(ji,jj,jk,2) == 0) EXIT
                  N_in = N_in + 1
                  tabin(jk) = ptab(ji,jj,jk,1)/ptab(ji,jj,jk,2)
                  h_in(N_in) = ptab(ji,jj,jk,2)/(e1v(ji,jj)*zrhox)
               END DO
               IF (N_in == 0) THEN
                  va(ji,jj,:) = 0._wp
                  CYCLE
               ENDIF
         
               N_out = 0
               DO jk=1,jpk
                  if (vmask(ji,jref,jk) == 0) EXIT
                  N_out = N_out + 1
                  h_out(N_out) = e3v_a(ji,jref,jk)
               END DO
               IF (N_out == 0) THEN
                 va(ji,jj,:) = 0._wp
                 CYCLE
               ENDIF
               call reconstructandremap(tabin(1:N_in),h_in(1:N_in),va(ji,jj,1:N_out),h_out(1:N_out),N_in,N_out)
            END DO
         END DO
# else
         DO jk = 1, jpkm1
          DO jj=j1,j2
          DO ji=i1,i2
              if (min(mbkt(ji,jj),mbkt(ji,jj+1))<jk) then
                e3v_0(ji,jj,jk)=e3t_1d(jk)
              else
                e3v_0(ji,jj,jk) = MAX(ptab(ji,jj,jk) / ( zrhox * e1v(ji,jj) ),MIN(e3zps_min,e3t_1d(jk)*e3zps_rat))
              endif
          ENDDO
          ENDDO
         END DO
# endif
      ENDIF
      !        
   END SUBROUTINE connect_e3v

#else
subroutine agrif_boundary_connections_empty
end subroutine agrif_boundary_connections_empty
#endif