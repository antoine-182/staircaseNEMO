MODULE agrif_oce
   !!======================================================================
   !!                       ***  MODULE agrif_oce  ***
   !! AGRIF :   define in memory AGRIF variables
   !!----------------------------------------------------------------------
   !! History :  2.0  ! 2007-12  (R. Benshila)  Original code
   !!----------------------------------------------------------------------
#if defined key_agrif
   !!----------------------------------------------------------------------
   !!   'key_agrif'                                              AGRIF zoom
   !!----------------------------------------------------------------------
   USE par_oce      ! ocean parameters
   USE dom_oce      ! domain parameters

   IMPLICIT NONE
   PRIVATE 

   PUBLIC agrif_oce_alloc ! routine called by nemo_init in nemogcm.F90
  
   !                                              !!* Namelist namagrif: AGRIF parameters
   LOGICAL , PUBLIC ::   ln_agrif_2way = .TRUE.    !: activate two way nesting 
   LOGICAL , PUBLIC ::   ln_spc_dyn    = .FALSE.   !: use zeros (.false.) or not (.true.) in
                                                   !: bdys dynamical fields interpolation
   REAL(wp), PUBLIC ::   rn_sponge_tra = 2800.     !: sponge coeff. for tracers
   REAL(wp), PUBLIC ::   rn_sponge_dyn = 2800.     !: sponge coeff. for dynamics
   REAL(wp), PUBLIC ::   rn_trelax_tra = 0.01      !: time relaxation parameter for tracers
   REAL(wp), PUBLIC ::   rn_trelax_dyn = 0.01      !: time relaxation parameter for momentum
   LOGICAL , PUBLIC ::   ln_chk_bathy  = .FALSE.   !: check of parent bathymetry 
   !
   INTEGER , PUBLIC, PARAMETER ::   nn_sponge_len = 2  !: Sponge width (in number of parent grid points)
   LOGICAL , PUBLIC :: spongedoneT = .FALSE.       !: tracer   sponge layer indicator
   LOGICAL , PUBLIC :: spongedoneU = .FALSE.       !: dynamics sponge layer indicator
   LOGICAL , PUBLIC :: lk_agrif_fstep = .TRUE.     !: if true: first step
   LOGICAL , PUBLIC :: lk_agrif_debug = .FALSE.    !: if true: print debugging info

   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tabspongedone_tsn
# if defined key_top
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tabspongedone_trn
# endif
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tabspongedone_u
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: tabspongedone_v
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: utint_stage
   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: vtint_stage
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: fspu, fspv !: sponge arrays
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: fspt, fspf !:   "      "

   ! Barotropic arrays used to store open boundary data during time-splitting loop:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  ubdy, vbdy, hbdy
   INTEGER , PUBLIC,              SAVE                 ::  Kbb_a, Kmm_a, Krhs_a   !: AGRIF module-specific copies of time-level indices

# if defined key_vertical
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: ht0_parent, hu0_parent, hv0_parent
   INTEGER,  PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: mbkt_parent, mbku_parent, mbkv_parent
# endif

   INTEGER, PUBLIC :: tsn_id                                                  ! AGRIF profile for tracers interpolation and update
   INTEGER, PUBLIC :: un_interp_id, vn_interp_id                              ! AGRIF profiles for interpolations
   INTEGER, PUBLIC :: un_update_id, vn_update_id                              ! AGRIF profiles for udpates
   INTEGER, PUBLIC :: tsn_sponge_id, un_sponge_id, vn_sponge_id               ! AGRIF profiles for sponge layers
# if defined key_top
   INTEGER, PUBLIC :: trn_id, trn_sponge_id
# endif  
   INTEGER, PUBLIC :: unb_id, vnb_id, ub2b_interp_id, vb2b_interp_id
   INTEGER, PUBLIC :: ub2b_update_id, vb2b_update_id
   INTEGER, PUBLIC :: e3t_id, e1u_id, e2v_id, sshn_id
   INTEGER, PUBLIC :: scales_t_id
   INTEGER, PUBLIC :: avt_id, avm_id, en_id                ! TKE related identificators
   INTEGER, PUBLIC :: mbkt_id, ht0_id
   INTEGER, PUBLIC :: kindic_agr
   
   !!----------------------------------------------------------------------
   !! NEMO/NST 4.0 , NEMO Consortium (2018)
   !! $Id: agrif_oce.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS 

   INTEGER FUNCTION agrif_oce_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION agrif_oce_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ierr
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( fspu(jpi,jpj), fspv(jpi,jpj),    	    &
         &      fspt(jpi,jpj), fspf(jpi,jpj),               &
         &      tabspongedone_tsn(jpi,jpj),                 &
         &      utint_stage(jpi,jpj), vtint_stage(jpi,jpj), &
# if defined key_top         
         &      tabspongedone_trn(jpi,jpj),           &
# endif   
# if defined key_vertical
         &      ht0_parent(jpi,jpj), mbkt_parent(jpi,jpj),  &
         &      hu0_parent(jpi,jpj), mbku_parent(jpi,jpj),  &
         &      hv0_parent(jpi,jpj), mbkv_parent(jpi,jpj),  &
# endif      
         &      tabspongedone_u  (jpi,jpj),           &
         &      tabspongedone_v  (jpi,jpj), STAT = ierr(1) )

      ALLOCATE( ubdy(jpi,jpj), vbdy(jpi,jpj), hbdy(jpi,jpj), STAT = ierr(2) )

      agrif_oce_alloc = MAXVAL(ierr)
      !
   END FUNCTION agrif_oce_alloc

#endif
   !!======================================================================
END MODULE agrif_oce
