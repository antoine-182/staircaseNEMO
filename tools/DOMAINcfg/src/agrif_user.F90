#if defined key_agrif
subroutine agrif_initworkspace()
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_InitWorkspace ***
      !!----------------------------------------------------------------------
end subroutine agrif_initworkspace
SUBROUTINE Agrif_InitValues
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_InitValues ***
      !!
      !! ** Purpose :: Declaration of variables to be interpolated
      !!----------------------------------------------------------------------
   USE Agrif_Util
   USE dom_oce
   USE nemogcm
   USE domain
   !!
   IMPLICIT NONE
   
 
   CALL nemo_init       !* Initializations of each fine grid

   CALL dom_nam
   CALL cfg_write         ! create the configuration file
   
END SUBROUTINE Agrif_InitValues

SUBROUTINE Agrif_InitValues_cont

use dom_oce
    integer :: irafx, irafy
    logical :: ln_perio
    integer nx,ny

irafx = agrif_irhox()
irafy = agrif_irhoy()

nx=nlci ; ny=nlcj

   !       IF(jperio /=1 .AND. jperio/=4 .AND. jperio/=6 ) THEN
   !          nx = (nbcellsx)+2*nbghostcellsfine+2
   !          ny = (nbcellsy)+2*nbghostcellsfine+2
   !          nbghostcellsfine_tot_x= nbghostcellsfine_x +1
   !          nbghostcellsfine_tot_y= nbghostcellsfine_y +1
   !       ELSE
   !         nx = (nbcellsx)+2*nbghostcellsfine_x 
   !         ny = (nbcellsy)+2*nbghostcellsfine+2
   !         nbghostcellsfine_tot_x= 1
   !         nbghostcellsfine_tot_y= nbghostcellsfine_y +1
   !      ENDIF
   !    ELSE
   !       nbghostcellsfine = 0
   !       nx = nbcellsx+irafx
   !       ny = nbcellsy+irafy
       
  WRITE(*,*) ' '
  WRITE(*,*)'Size of the High resolution grid: ',nx,' x ',ny
  WRITE(*,*) ' '
       
       call agrif_init_lonlat()
       ln_perio=.FALSE. 
       if( jperio ==1 .OR. jperio==2 .OR. jperio==4) ln_perio=.TRUE.

       where (glamt < -180) glamt = glamt +360.
       if (ln_perio) then
         glamt(1,:)=glamt(nx-1,:)
         glamt(nx,:)=glamt(2,:)
       endif
 
       where (glamu < -180) glamu = glamu +360.
       if (ln_perio) then
         glamu(1,:)=glamu(nx-1,:)
         glamu(nx,:)=glamu(2,:)
       endif

      where (glamv < -180) glamv = glamv +360.
       if (ln_perio) then
         glamv(1,:)=glamv(nx-1,:)
         glamv(nx,:)=glamv(2,:)
       endif

      where (glamf < -180) glamf = glamf +360.
       if (ln_perio) then
         glamf(1,:)=glamf(nx-1,:)
         glamf(nx,:)=glamf(2,:)
       endif

       call agrif_init_scales()

        
END SUBROUTINE Agrif_InitValues_cont  


subroutine agrif_declare_var()
use par_oce
use dom_oce
use agrif_profiles
use agrif_parameters

   IMPLICIT NONE
   
integer :: ind1, ind2, ind3
integer nx,ny
integer nbghostcellsfine_tot_x, nbghostcellsfine_tot_y
INTEGER :: irafx
!!----------------------------------------------------------------------

   ! 1. Declaration of the type of variable which have to be interpolated
   !---------------------------------------------------------------------
 nx=nlci ; ny=nlcj

!ind2 = nbghostcellsfine_tot_x + 1
!ind3 = nbghostcellsfine_tot_y + 1
ind2 = 2 + nbghostcells
ind3 = ind2
nbghostcellsfine_tot_x=nbghostcells+1
nbghostcellsfine_tot_y=nbghostcells+1

irafx = Agrif_irhox()

CALL agrif_nemo_init  ! specific namelist part if needed

CALL agrif_declare_variable((/2,2/),(/ind2,ind3/),(/'x','y'/),(/1,1/),(/nx,ny/),glamt_id)
CALL agrif_declare_variable((/1,2/),(/ind2-1,ind3/),(/'x','y'/),(/1,1/),(/nx,ny/),glamu_id)
CALL agrif_declare_variable((/2,1/),(/ind2,ind3-1/),(/'x','y'/),(/1,1/),(/nx,ny/),glamv_id)
CALL agrif_declare_variable((/1,1/),(/ind2-1,ind3-1/),(/'x','y'/),(/1,1/),(/nx,ny/),glamf_id)

CALL agrif_declare_variable((/2,2/),(/ind2,ind3/),(/'x','y'/),(/1,1/),(/nx,ny/),gphit_id)
CALL agrif_declare_variable((/1,2/),(/ind2-1,ind3/),(/'x','y'/),(/1,1/),(/nx,ny/),gphiu_id)
CALL agrif_declare_variable((/2,1/),(/ind2,ind3-1/),(/'x','y'/),(/1,1/),(/nx,ny/),gphiv_id)
CALL agrif_declare_variable((/1,1/),(/ind2-1,ind3-1/),(/'x','y'/),(/1,1/),(/nx,ny/),gphif_id)

! Horizontal scale factors

CALL agrif_declare_variable((/2,2/),(/ind2,ind3/),(/'x','y'/),(/1,1/),(/nx,ny/),e1t_id)
CALL agrif_declare_variable((/1,2/),(/ind2-1,ind3/),(/'x','y'/),(/1,1/),(/nx,ny/),e1u_id)
CALL agrif_declare_variable((/2,1/),(/ind2,ind3-1/),(/'x','y'/),(/1,1/),(/nx,ny/),e1v_id)
CALL agrif_declare_variable((/1,1/),(/ind2-1,ind3-1/),(/'x','y'/),(/1,1/),(/nx,ny/),e1f_id)

CALL agrif_declare_variable((/2,2/),(/ind2,ind3/),(/'x','y'/),(/1,1/),(/nx,ny/),e2t_id)
CALL agrif_declare_variable((/1,2/),(/ind2-1,ind3/),(/'x','y'/),(/1,1/),(/nx,ny/),e2u_id)
CALL agrif_declare_variable((/2,1/),(/ind2,ind3-1/),(/'x','y'/),(/1,1/),(/nx,ny/),e2v_id)
CALL agrif_declare_variable((/1,1/),(/ind2-1,ind3-1/),(/'x','y'/),(/1,1/),(/nx,ny/),e2f_id)

! Bathymetry

CALL agrif_declare_variable((/2,2/),(/ind2,ind3/),(/'x','y'/),(/1,1/),(/nx,ny/),bathy_id)

! Vertical scale factors
CALL agrif_declare_variable((/2,2,0/),(/ind2,ind3,0/),(/'x','y','N'/),(/1,1,1/),(/nx,ny,jpk/),e3t_id)
CALL agrif_declare_variable((/2,2,0/),(/ind2,ind3,0/),(/'x','y','N'/),(/1,1,1/),(/nx,ny,jpk/),e3t_copy_id)
CALL agrif_declare_variable((/2,2,0/),(/ind2,ind3,0/),(/'x','y','N'/),(/1,1,1/),(/nx,ny,jpk+1/),e3t_connect_id)

CALL agrif_declare_variable((/1,2,0/),(/ind2-1,ind3,0/),(/'x','y','N'/),(/1,1,1/),(/nx,ny,jpk/),e3u_id)
CALL agrif_declare_variable((/2,1,0/),(/ind2,ind3-1,0/),(/'x','y','N'/),(/1,1,1/),(/nx,ny,jpk/),e3v_id)

! Bottom level

CALL agrif_declare_variable((/2,2/),(/ind2,ind3/),(/'x','y'/),(/1,1/),(/nx,ny/),bottom_level_id)

CALL Agrif_Set_bcinterp(glamt_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(glamt_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( glamt_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(glamu_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(glamu_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( glamu_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(glamv_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(glamv_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( glamv_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(glamf_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(glamf_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( glamf_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(gphit_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(gphit_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( gphit_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(gphiu_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(gphiu_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( gphiu_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(gphiv_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(gphiv_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( gphiv_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(gphif_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(gphif_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( gphif_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

!

CALL Agrif_Set_bcinterp(e1t_id,interp=AGRIF_ppm)
CALL Agrif_Set_interp(e1t_id,interp=AGRIF_ppm)
CALL Agrif_Set_bc( e1t_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(e1u_id, interp1=Agrif_linear, interp2=AGRIF_ppm)
CALL Agrif_Set_interp(e1u_id, interp1=Agrif_linear, interp2=AGRIF_ppm)
CALL Agrif_Set_bc( e1u_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(e1v_id,interp1=AGRIF_ppm, interp2=Agrif_linear)
CALL Agrif_Set_interp(e1v_id, interp1=AGRIF_ppm, interp2=Agrif_linear)
CALL Agrif_Set_bc( e1v_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(e1f_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(e1f_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( e1f_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(e2t_id,interp=AGRIF_ppm)
CALL Agrif_Set_interp(e2t_id,interp=AGRIF_ppm)
CALL Agrif_Set_bc( e2t_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(e2u_id,interp1=Agrif_linear, interp2=AGRIF_ppm)
CALL Agrif_Set_interp(e2u_id,interp1=Agrif_linear, interp2=AGRIF_ppm)
CALL Agrif_Set_bc( e2u_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(e2v_id,interp1=AGRIF_ppm, interp2=Agrif_linear)
CALL Agrif_Set_interp(e2v_id,interp1=AGRIF_ppm, interp2=Agrif_linear)
CALL Agrif_Set_bc( e2v_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(e2f_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(e2f_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( e2f_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

CALL Agrif_Set_bcinterp(bathy_id,interp=AGRIF_linear)
CALL Agrif_Set_interp(bathy_id,interp=AGRIF_linear)
CALL Agrif_Set_bc( bathy_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )

! Vertical scale factors
CALL Agrif_Set_bcinterp(e3t_id,interp=AGRIF_ppm)
CALL Agrif_Set_interp(e3t_id,interp=AGRIF_ppm)
CALL Agrif_Set_bc( e3t_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )
CALL Agrif_Set_Updatetype( e3t_id, update = AGRIF_Update_Average)

CALL Agrif_Set_bcinterp(e3t_copy_id,interp=AGRIF_constant)
CALL Agrif_Set_interp(e3t_copy_id,interp=AGRIF_constant)
CALL Agrif_Set_bc( e3t_copy_id, (/-npt_copy*irafx-1,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/))

CALL Agrif_Set_bcinterp(e3t_connect_id,interp=AGRIF_ppm)
CALL Agrif_Set_interp(e3t_connect_id,interp=AGRIF_ppm)
CALL Agrif_Set_bc( e3t_connect_id, (/-(npt_copy+npt_connect)*irafx-1,-npt_copy*irafx-2/))

CALL Agrif_Set_bcinterp(e3u_id, interp1=Agrif_linear, interp2=AGRIF_ppm)
CALL Agrif_Set_interp(e3u_id, interp1=Agrif_linear, interp2=AGRIF_ppm)
CALL Agrif_Set_bc( e3u_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )
CALL Agrif_Set_Updatetype(e3u_id,update1 = Agrif_Update_Copy, update2 = Agrif_Update_Average)

CALL Agrif_Set_bcinterp(e3v_id,interp1=AGRIF_ppm, interp2=Agrif_linear)
CALL Agrif_Set_interp(e3v_id, interp1=AGRIF_ppm, interp2=Agrif_linear)
CALL Agrif_Set_bc( e3v_id, (/0,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/) )
CALL Agrif_Set_Updatetype(e3v_id,update1 = Agrif_Update_Average, update2 = Agrif_Update_Copy)
   
! Bottom level
CALL Agrif_Set_bcinterp(bottom_level_id,interp=AGRIF_constant)
CALL Agrif_Set_interp(bottom_level_id,interp=AGRIF_constant)
CALL Agrif_Set_bc( bottom_level_id, (/-npt_copy*irafx-1,max(nbghostcellsfine_tot_x,nbghostcellsfine_tot_y)-1/))
CALL Agrif_Set_Updatetype( bottom_level_id, update = AGRIF_Update_Max)

end subroutine agrif_declare_var


subroutine agrif_init_lonlat()
use agrif_profiles
use agrif_util
external :: init_glamt, init_glamu, init_glamv, init_glamf
external :: init_gphit, init_gphiu, init_gphiv, init_gphif

call Agrif_Init_variable(glamt_id, procname = init_glamt)
call Agrif_Init_variable(glamu_id, procname = init_glamu)
call Agrif_Init_variable(glamv_id, procname = init_glamv)
call Agrif_Init_variable(glamf_id, procname = init_glamf)

call Agrif_Init_variable(gphit_id, procname = init_gphit)
call Agrif_Init_variable(gphiu_id, procname = init_gphiu)
call Agrif_Init_variable(gphiv_id, procname = init_gphiv)
call Agrif_Init_variable(gphif_id, procname = init_gphif)

end subroutine agrif_init_lonlat

subroutine agrif_init_scales()
use agrif_profiles
use agrif_util
external :: init_e1t, init_e1u, init_e1v, init_e1f
external :: init_e2t, init_e2u, init_e2v, init_e2f

call Agrif_Init_variable(e1t_id, procname = init_e1t)
call Agrif_Init_variable(e1u_id, procname = init_e1u)
call Agrif_Init_variable(e1v_id, procname = init_e1v)
call Agrif_Init_variable(e1f_id, procname = init_e1f)

call Agrif_Init_variable(e2t_id, procname = init_e2t)
call Agrif_Init_variable(e2u_id, procname = init_e2u)
call Agrif_Init_variable(e2v_id, procname = init_e2v)
call Agrif_Init_variable(e2f_id, procname = init_e2f)

end subroutine agrif_init_scales



   SUBROUTINE init_glamt( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = glamt(i1:i2,j1:j2)
      ELSE
         glamt(i1:i2,j1:j2)=ptab
      ENDIF
      !
   END SUBROUTINE init_glamt

    SUBROUTINE init_glamu( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = glamu(i1:i2,j1:j2)
      ELSE
         glamu(i1:i2,j1:j2)=ptab
      ENDIF
      !
   END SUBROUTINE init_glamu

   SUBROUTINE init_glamv( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = glamv(i1:i2,j1:j2)
      ELSE
         glamv(i1:i2,j1:j2)=ptab
      ENDIF
      !
   END SUBROUTINE init_glamv

   SUBROUTINE init_glamf( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = glamf(i1:i2,j1:j2)
      ELSE
         glamf(i1:i2,j1:j2)=ptab
      ENDIF
      !
   END SUBROUTINE init_glamf

   SUBROUTINE init_gphit( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = gphit(i1:i2,j1:j2)
      ELSE
         gphit(i1:i2,j1:j2)=ptab
      ENDIF
      !
   END SUBROUTINE init_gphit

    SUBROUTINE init_gphiu( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = gphiu(i1:i2,j1:j2)
      ELSE
         gphiu(i1:i2,j1:j2)=ptab
      ENDIF
      !
   END SUBROUTINE init_gphiu

    SUBROUTINE init_gphiv( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = gphiv(i1:i2,j1:j2)
      ELSE
         gphiv(i1:i2,j1:j2)=ptab
      ENDIF
      !
   END SUBROUTINE init_gphiv


   SUBROUTINE init_gphif( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = gphif(i1:i2,j1:j2)
      ELSE
         gphif(i1:i2,j1:j2)=ptab
      ENDIF
      !
   END SUBROUTINE init_gphif


   SUBROUTINE init_e1t( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = e1t(i1:i2,j1:j2)
      ELSE
         e1t(i1:i2,j1:j2)=ptab/Agrif_rhoy()
      ENDIF
      !
   END SUBROUTINE init_e1t

   SUBROUTINE init_e1u( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = e1u(i1:i2,j1:j2)
      ELSE
         e1u(i1:i2,j1:j2)=ptab/Agrif_rhoy()
      ENDIF
      !
   END SUBROUTINE init_e1u

   SUBROUTINE init_e1v( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = e1v(i1:i2,j1:j2)
      ELSE
         e1v(i1:i2,j1:j2)=ptab/Agrif_rhoy()
      ENDIF
      !
   END SUBROUTINE init_e1v

   SUBROUTINE init_e1f( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = e1f(i1:i2,j1:j2)
      ELSE
         e1f(i1:i2,j1:j2)=ptab/Agrif_rhoy()
      ENDIF
      !
   END SUBROUTINE init_e1f

  SUBROUTINE init_e2t( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = e2t(i1:i2,j1:j2)
      ELSE
         e2t(i1:i2,j1:j2)=ptab/Agrif_rhoy()
      ENDIF
      !
   END SUBROUTINE init_e2t

   SUBROUTINE init_e2u( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = e2u(i1:i2,j1:j2)
      ELSE
         e2u(i1:i2,j1:j2)=ptab/Agrif_rhoy()
      ENDIF
      !
   END SUBROUTINE init_e2u

   SUBROUTINE init_e2v( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = e2v(i1:i2,j1:j2)
      ELSE
         e2v(i1:i2,j1:j2)=ptab/Agrif_rhoy()
      ENDIF
      !
   END SUBROUTINE init_e2v

   SUBROUTINE init_e2f( ptab, i1, i2, j1, j2, before, nb,ndir)
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
         ptab(i1:i2,j1:j2) = e2f(i1:i2,j1:j2)
      ELSE
         e2f(i1:i2,j1:j2)=ptab/Agrif_rhoy()
      ENDIF
      !
   END SUBROUTINE init_e2f


SUBROUTINE agrif_nemo_init
USE agrif_parameters
USE in_out_manager
USE lib_mpp

   
   !!
   IMPLICIT NONE
   
   INTEGER ::   ios
   
   NAMELIST/namagrif/ nn_cln_update,ln_spc_dyn,rn_sponge_tra,rn_sponge_dyn,ln_chk_bathy,npt_connect, npt_copy

      REWIND( numnam_ref )              ! Namelist namagrif in reference namelist : nesting parameters
      READ  ( numnam_ref, namagrif, IOSTAT = ios, ERR = 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namagrif in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist namzgr in configuration namelist : nesting parameters
      READ  ( numnam_cfg, namagrif, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namagrif in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, namagrif )

      IF(lwp) THEN                     ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'agrif_nemo_init : nesting'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namagrif : set nesting parameters'
         WRITE(numout,*) '      npt_copy     = ', npt_copy
         WRITE(numout,*) '      npt_connect  = ', npt_connect
      ENDIF
      
END SUBROUTINE agrif_nemo_init
   

SUBROUTINE Agrif_detect( kg, ksizex )
      !!----------------------------------------------------------------------
      !!                      *** ROUTINE Agrif_detect ***
      !!----------------------------------------------------------------------
   INTEGER, DIMENSION(2) :: ksizex
   INTEGER, DIMENSION(ksizex(1),ksizex(2)) :: kg 
      !!----------------------------------------------------------------------
   !
   RETURN
   !
END SUBROUTINE Agrif_detect
SUBROUTINE agrif_before_regridding
END SUBROUTINE agrif_before_regridding

SUBROUTINE Agrif_InvLoc( indloc, nprocloc, i, indglob )
      !!----------------------------------------------------------------------
      !!                     *** ROUTINE Agrif_InvLoc ***
      !!----------------------------------------------------------------------
   USE dom_oce
   !!
   IMPLICIT NONE
   !
   INTEGER :: indglob, indloc, nprocloc, i
      !!----------------------------------------------------------------------
   !
   SELECT CASE( i )
   CASE(1)   ;   indglob = indloc + nimppt(nprocloc+1) - 1
   CASE(2)   ;   indglob = indloc + njmppt(nprocloc+1) - 1
   CASE DEFAULT
      indglob = indloc
   END SELECT
   !
END SUBROUTINE Agrif_InvLoc

SUBROUTINE Agrif_get_proc_info( imin, imax, jmin, jmax )
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_get_proc_info ***
      !!----------------------------------------------------------------------
   USE par_oce
   USE dom_oce 
   !!
   IMPLICIT NONE
   !
   INTEGER, INTENT(out) :: imin, imax
   INTEGER, INTENT(out) :: jmin, jmax
      !!----------------------------------------------------------------------
   !
   imin = nimppt(Agrif_Procrank+1)  ! ?????
   jmin = njmppt(Agrif_Procrank+1)  ! ?????
   imax = imin + jpi - 1
   jmax = jmin + jpj - 1
   ! 
END SUBROUTINE Agrif_get_proc_info

SUBROUTINE Agrif_estimate_parallel_cost(imin, imax,jmin, jmax, nbprocs, grid_cost)
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_estimate_parallel_cost ***
      !!----------------------------------------------------------------------
   USE par_oce
   !!
   IMPLICIT NONE
   !
   INTEGER,  INTENT(in)  :: imin, imax
   INTEGER,  INTENT(in)  :: jmin, jmax
   INTEGER,  INTENT(in)  :: nbprocs
   REAL(wp), INTENT(out) :: grid_cost
      !!----------------------------------------------------------------------
   !
   grid_cost = REAL(imax-imin+1,wp)*REAL(jmax-jmin+1,wp) / REAL(nbprocs,wp)
   !
END SUBROUTINE Agrif_estimate_parallel_cost

#endif
