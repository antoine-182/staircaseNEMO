MODULE traisf
   !!==============================================================================
   !!                       ***  MODULE  traisf  ***
   !! Ocean active tracers:  ice shelf boundary condition
   !!==============================================================================
   !! History :    4.0  !  2019-09  (P. Mathiot) original file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_isf       : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   USE isf_oce                                     ! Ice shelf variables
   USE dom_oce                                     ! ocean space domain variables
   USE isfutils, ONLY : debug                      ! debug option
   USE timing  , ONLY : timing_start, timing_stop  ! Timing
   USE in_out_manager                              ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_isf   ! routine called by step.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: trasbc.F90 10499 2019-01-10 15:12:24Z deazer $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_isf ( kt, Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_isf  ***
      !!
      !! ** Purpose :  Compute the temperature trend due to the ice shelf melting (qhoce + qhc)
      !!
      !! ** Action  : - update pts(:,:,:,:,Krhs) for cav, par and cpl case
      !!----------------------------------------------------------------------
      INTEGER                                  , INTENT(in   ) :: kt        ! ocean time step
      INTEGER                                  , INTENT(in   ) :: Kmm, Krhs ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) :: pts       ! active tracers and RHS of tracer equation
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_isf')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_isf : Ice shelf heat fluxes'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      ! cavity case
      IF ( ln_isfcav_mlt ) CALL tra_isf_mlt(misfkt_cav, misfkb_cav, rhisf_tbl_cav, rfrac_tbl_cav, risf_cav_tsc, risf_cav_tsc_b, pts(:,:,:,:,Krhs))
      !
      ! parametrisation case
      IF ( ln_isfpar_mlt ) CALL tra_isf_mlt(misfkt_par, misfkb_par, rhisf_tbl_par, rfrac_tbl_par, risf_par_tsc, risf_par_tsc_b, pts(:,:,:,:,Krhs))
      !
      ! ice sheet coupling case
      IF ( ln_isfcpl ) THEN
         !
         ! Dynamical stability at start up after change in under ice shelf cavity geometry is achieve by correcting the divergence.
         ! This is achieved by applying a volume flux in order to keep the horizontal divergence after remapping
         ! the same as at the end of the latest time step. So correction need to be apply at nit000 (euler time step) and
         ! half of it at nit000+1 (leap frog time step).
         ! in accordance to this, the heat content flux due to injected water need to be added in the temperature and salt trend
         ! at time step nit000 and nit000+1
         IF ( kt == nit000  ) CALL tra_isf_cpl(Kmm, risfcpl_tsc       , pts(:,:,:,:,Krhs))
         IF ( kt == nit000+1) CALL tra_isf_cpl(Kmm, risfcpl_tsc*0.5_wp, pts(:,:,:,:,Krhs))
         !
         ! ensure 0 trend due to unconservation of the ice shelf coupling
         IF ( ln_isfcpl_cons ) CALL tra_isf_cpl(Kmm, risfcpl_cons_tsc, pts(:,:,:,:,Krhs))
         !
      END IF
      !
      IF ( ln_isfdebug ) THEN
         CALL debug('tra_isf: pts(:,:,:,:,Krhs) T', pts(:,:,:,1,Krhs))
         CALL debug('tra_isf: pts(:,:,:,:,Krhs) S', pts(:,:,:,2,Krhs))
      END IF
      !
      IF( ln_timing )   CALL timing_stop('tra_isf')
      !
   END SUBROUTINE tra_isf
   !
   SUBROUTINE tra_isf_mlt(ktop, kbot, phtbl, pfrac, ptsc, ptsc_b, pts)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_isf_mlt  ***
      !!
      !! *** Purpose :  Compute the temperature trend due to the ice shelf melting (qhoce + qhc) for cav or par case
      !!
      !! *** Action :: Update pts(:,:,:,:,Krhs) with the surface boundary condition trend
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) :: pts
      !!----------------------------------------------------------------------
      INTEGER , DIMENSION(jpi,jpj)     , INTENT(in   ) :: ktop , kbot
      REAL(wp), DIMENSION(jpi,jpj)     , INTENT(in   ) :: phtbl, pfrac
      REAL(wp), DIMENSION(jpi,jpj,jpts), INTENT(in   ) :: ptsc , ptsc_b
      !!----------------------------------------------------------------------
      INTEGER                      :: ji,jj,jk  ! loop index
      INTEGER                      :: ikt, ikb  ! top and bottom level of the tbl
      REAL(wp), DIMENSION(jpi,jpj) :: ztc       ! total ice shelf tracer trend
      !!----------------------------------------------------------------------
      !
      ! compute 2d total trend due to isf
      ztc(:,:) = 0.5_wp * ( ptsc(:,:,jp_tem) + ptsc_b(:,:,jp_tem) ) / phtbl(:,:)
      !
      ! update pts(:,:,:,:,Krhs)
      DO_2D_11_11
         !
         ikt = ktop(ji,jj)
         ikb = kbot(ji,jj)
         !
         ! level fully include in the ice shelf boundary layer
         DO jk = ikt, ikb - 1
            pts(ji,jj,jk,jp_tem) = pts(ji,jj,jk,jp_tem) + ztc(ji,jj)
         END DO
         !
         ! level partially include in ice shelf boundary layer
         pts(ji,jj,ikb,jp_tem) = pts(ji,jj,ikb,jp_tem) + ztc(ji,jj) * pfrac(ji,jj)
         !
      END_2D
      !
   END SUBROUTINE tra_isf_mlt
   !
   SUBROUTINE tra_isf_cpl( Kmm, ptsc, ptsa )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_isf_cpl  ***
      !!
      !! *** Action :: Update pts(:,:,:,:,Krhs) with the ice shelf coupling trend
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) :: ptsa
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) :: Kmm   ! ocean time level index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) :: ptsc
      !!----------------------------------------------------------------------
      INTEGER :: jk
      !!----------------------------------------------------------------------
      !
      DO jk = 1,jpk
         ptsa(:,:,jk,jp_tem) =   &
            &  ptsa(:,:,jk,jp_tem) + ptsc(:,:,jk,jp_tem) * r1_e1e2t(:,:) / e3t(:,:,jk,Kmm)
         ptsa(:,:,jk,jp_sal) =   &
            &  ptsa(:,:,jk,jp_sal) + ptsc(:,:,jk,jp_sal) * r1_e1e2t(:,:) / e3t(:,:,jk,Kmm)
      END DO
      !
   END SUBROUTINE tra_isf_cpl
   !
END MODULE traisf
