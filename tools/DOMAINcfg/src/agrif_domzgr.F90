#if defined key_agrif
subroutine agrif_domzgr
end subroutine agrif_domzgr

subroutine agrif_create_bathy_meter
use agrif_profiles
external :: init_bathy

call Agrif_Init_variable(bathy_id, procname = init_bathy)

end subroutine agrif_create_bathy_meter

    SUBROUTINE init_bathy( ptab, i1, i2, j1, j2, before, nb,ndir)
    use dom_oce
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpsshn  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL, DIMENSION(i1:i2,j1:j2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      INTEGER                         , INTENT(in   ) ::   nb , ndir
      LOGICAL  ::   western_side, eastern_side,northern_side,southern_side
      !
      !!----------------------------------------------------------------------  
      !
         western_side  = (nb == 1).AND.(ndir == 1)
         eastern_side  = (nb == 1).AND.(ndir == 2)
         southern_side = (nb == 2).AND.(ndir == 1)
         northern_side = (nb == 2).AND.(ndir == 2)
      IF( before) THEN
         ptab(i1:i2,j1:j2) = bathy(i1:i2,j1:j2)
      ELSE
         bathy(i1:i2,j1:j2)=ptab
      ENDIF
      !
   END SUBROUTINE init_bathy
#else
subroutine agrif_domzgr_empty
end subroutine agrif_domzgr_empty
#endif