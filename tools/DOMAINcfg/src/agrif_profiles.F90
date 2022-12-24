module agrif_profiles

integer :: glamt_id, glamu_id, glamv_id,glamf_id
integer :: gphit_id, gphiu_id, gphiv_id,gphif_id
integer :: e1t_id, e1u_id, e1v_id, e1f_id
integer :: e2t_id, e2u_id, e2v_id, e2f_id


integer :: bathy_id

! Vertical scale factors

integer :: e3t_id
integer :: e3t_copy_id
integer :: e3t_connect_id
integer :: e3u_id, e3v_id

! Bottom level
integer :: bottom_level_id
end module agrif_profiles