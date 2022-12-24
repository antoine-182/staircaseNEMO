#if defined key_agrif
subroutine agrif_update_all
USE agrif_parameters
USE agrif_profiles
external update_bottom_level, update_e3t, update_e3u, update_e3v

if (agrif_root()) return
call agrif_update_variable(bottom_level_id,locupdate=(/npt_copy,0/),procname = update_bottom_level)

      Agrif_UseSpecialValueInUpdate = .TRUE.
      Agrif_SpecialValueFineGrid    = 0._wp
      
call agrif_update_variable(e3t_id,procname = update_e3t)
      Agrif_UseSpecialValueInUpdate = .FALSE.

call agrif_update_variable(e3u_id,procname = update_e3u)
call agrif_update_variable(e3v_id,procname = update_e3v)
      
end subroutine agrif_update_all

    SUBROUTINE update_bottom_level( ptab, i1, i2, j1, j2, before, nb,ndir)
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
   END SUBROUTINE update_bottom_level
   
   SUBROUTINE update_e3t( tabres, i1, i2, j1, j2, k1, k2,  before )
   USE dom_oce
   implicit none
      !!---------------------------------------------
      !!           *** update_e3t updateT ***
      !!---------------------------------------------
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji,jj,jk
      !!---------------------------------------------
      !
      IF (before) THEN
            DO jk=k1,k2
               DO jj=j1,j2
                  DO ji=i1,i2
                   if (mbkt(ji,jj) < jk) then
                     tabres(ji,jj,jk) = e3t_0(ji,jj,jk)
                   else
                     tabres(ji,jj,jk) = 0.
                   endif
                  END DO
               END DO
            END DO
      ELSE
            DO jk=k1,k2
               DO jj=j1,j2
                  DO ji=i1,i2
                   if (mbkt(ji,jj) < jk) then
                     e3t_0(ji,jj,jk) = e3t_1d(jk)
                   else
                     e3t_0(ji,jj,jk) = MAX(tabres(ji,jj,jk),MIN(e3zps_min,e3t_1d(jk)*e3zps_rat))
                   endif
                  END DO
               END DO
            END DO
         !
      ENDIF
      ! 
   END SUBROUTINE update_e3t
   
   SUBROUTINE update_e3u( tabres, i1, i2, j1, j2, k1, k2, before )
   USE dom_oce
   implicit none
      !!---------------------------------------------
      !!           *** ROUTINE update_e3u ***
      !!---------------------------------------------
      INTEGER                               , INTENT(in   ) :: i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      LOGICAL                                     , INTENT(in   ) :: before
      !
      INTEGER  :: ji, jj, jk
      REAL :: zrhoy
      !!---------------------------------------------
      ! 
      IF( before ) THEN
         zrhoy = Agrif_Rhoy()
         DO jk = k1, k2
          do jj=j1,j2
          do ji=i1,i2
           if (min(mbkt(ji,jj),mbkt(ji+1,jj))<jk) then
            tabres(ji,jj,jk) = zrhoy * e2u(ji,jj) * MIN(e3zps_min,e3t_1d(jk)*e3zps_rat)
           else
            tabres(ji,jj,jk) = zrhoy * e2u(ji,jj) * e3u_0(ji,jj,jk)
           endif
          enddo
          enddo
         END DO
      ELSE
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                 if (min(mbkt(ji,jj),mbkt(ji+1,jj))<jk) then
                   e3u_0(ji,jj,jk)=e3t_1d(jk)
                 else
                   e3u_0(ji,jj,jk) = MAX(tabres(ji,jj,jk) / e2u(ji,jj),MIN(e3zps_min,e3t_1d(jk)*e3zps_rat))
                 endif
               END DO
            END DO
         END DO
         !
      ENDIF
      ! 
   END SUBROUTINE update_e3u
   
   SUBROUTINE update_e3v( tabres, i1, i2, j1, j2, k1, k2, before )
   USE dom_oce
   implicit none
      !!---------------------------------------------
      !!           *** ROUTINE update_e3v ***
      !!---------------------------------------------
      INTEGER                               , INTENT(in   ) :: i1, i2, j1, j2, k1, k2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: tabres
      LOGICAL                                     , INTENT(in   ) :: before
      !
      INTEGER  :: ji, jj, jk
      REAL :: zrhox
      !!---------------------------------------------
      ! 
      IF( before ) THEN
         zrhox = Agrif_Rhox()
         DO jk = k1, k2
          do jj=j1,j2
          do ji=i1,i2
           if (min(mbkt(ji,jj),mbkt(ji,jj+1))<jk) then
            tabres(ji,jj,jk) = zrhox * e1v(ji,jj) * MIN(e3zps_min,e3t_1d(jk)*e3zps_rat)
           else
            tabres(ji,jj,jk) = zrhox * e1v(ji,jj) * e3v_0(ji,jj,jk)
           endif
          enddo
          enddo
         END DO
      ELSE
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                 if (min(mbkt(ji,jj),mbkt(ji,jj+1))<jk) then
                   e3v_0(ji,jj,jk)=e3t_1d(jk)
                 else
                   e3v_0(ji,jj,jk) = MAX(tabres(ji,jj,jk) / e1v(ji,jj),MIN(e3zps_min,e3t_1d(jk)*e3zps_rat))
                 endif
               END DO
            END DO
         END DO
         !
      ENDIF
      ! 
   END SUBROUTINE update_e3v
   
#else
subroutine agrif_update_all_empty
end subroutine agrif_update_all_empty
#endif