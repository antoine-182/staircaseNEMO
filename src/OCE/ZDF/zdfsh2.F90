MODULE zdfsh2
   !!======================================================================
   !!                       ***  MODULE  zdfsh2  ***
   !! Ocean physics:  shear production term of TKE 
   !!=====================================================================
   !! History :   -   !  2014-10  (A. Barthelemy, G. Madec)  original code
   !!   NEMO     4.0  !  2017-04  (G. Madec)  remove u-,v-pts avm
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_sh2       : compute mixing the shear production term of TKE
   !!----------------------------------------------------------------------
   USE oce
   USE dom_oce        ! domain: ocean
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_sh2        ! called by zdftke, zdfglf, and zdfric
   
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: zdfsh2.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_sh2( Kbb, Kmm, p_avm, p_sh2  ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_sh2  ***
      !!
      !! ** Purpose :   Compute the shear production term of a TKE equation
      !!
      !! ** Method  : - a stable discretization of this term is linked to the
      !!                time-space discretization of the vertical diffusion
      !!                of the OGCM. NEMO uses C-grid, a leap-frog environment 
      !!                and an implicit computation of vertical mixing term,
      !!                so the shear production at w-point is given by:
      !!                   sh2 = mi[   mi(avm) * dk[ub]/e3ub * dk[un]/e3un   ] 
      !!                       + mj[   mj(avm) * dk[vb]/e3vb * dk[vn]/e3vn   ] 
      !!                NB: wet-point only horizontal averaging of shear
      !!
      !! ** Action  : - p_sh2 shear prod. term at w-point (inner domain only)
      !!                                                   *****
      !! References :   Bruchard, OM 2002
      !! ---------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   Kbb, Kmm             ! ocean time level indices
      REAL(wp), DIMENSION(:,:,:) , INTENT(in   ) ::   p_avm                ! vertical eddy viscosity (w-points)
      REAL(wp), DIMENSION(:,:,:) , INTENT(  out) ::   p_sh2                ! shear production of TKE (w-points)
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop arguments
      REAL(wp), DIMENSION(jpi,jpj) ::   zsh2u, zsh2v   ! 2D workspace
      !!--------------------------------------------------------------------
      !
      DO jk = 2, jpkm1
         DO_2D_10_10
            zsh2u(ji,jj) = ( p_avm(ji+1,jj,jk) + p_avm(ji,jj,jk) ) &
               &         * (   uu(ji,jj,jk-1,Kmm) -   uu(ji,jj,jk,Kmm) ) &
               &         * (   uu(ji,jj,jk-1,Kbb) -   uu(ji,jj,jk,Kbb) ) & 
               &         / ( e3uw(ji,jj,jk  ,Kmm) * e3uw(ji,jj,jk,Kbb) ) &
               &         * wumask(ji,jj,jk)
            zsh2v(ji,jj) = ( p_avm(ji,jj+1,jk) + p_avm(ji,jj,jk) ) &
               &         * (   vv(ji,jj,jk-1,Kmm) -   vv(ji,jj,jk,Kmm) ) &
               &         * (   vv(ji,jj,jk-1,Kbb) -   vv(ji,jj,jk,Kbb) ) &
               &         / ( e3vw(ji,jj,jk  ,Kmm) * e3vw(ji,jj,jk,Kbb) ) &
               &         * wvmask(ji,jj,jk)
         END_2D
         DO_2D_00_00
            p_sh2(ji,jj,jk) = 0.25 * (   ( zsh2u(ji-1,jj) + zsh2u(ji,jj) ) * ( 2. - umask(ji-1,jj,jk) * umask(ji,jj,jk) )   &
               &                       + ( zsh2v(ji,jj-1) + zsh2v(ji,jj) ) * ( 2. - vmask(ji,jj-1,jk) * vmask(ji,jj,jk) )   )
         END_2D
      END DO 
      !
   END SUBROUTINE zdf_sh2

   !!======================================================================
END MODULE zdfsh2
