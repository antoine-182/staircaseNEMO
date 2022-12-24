#ifdef key_agrif
module agrif_parameters
USE par_kind

INTEGER :: nn_cln_update
LOGICAL :: ln_spc_dyn
REAL(wp) :: rn_sponge_tra
REAL(wp) :: rn_sponge_dyn
LOGICAL :: ln_chk_bathy
INTEGER :: npt_copy
INTEGER :: npt_connect

REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:) ::   ztabramp

end module agrif_parameters
#else
subroutine agrif_parameters_empty
end subroutine agrif_parameters_empty
#endif