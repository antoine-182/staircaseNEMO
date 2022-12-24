MODULE agrif_oce_interp
   !!======================================================================
   !!                   ***  MODULE  agrif_oce_interp  ***
   !! AGRIF: interpolation package for the ocean dynamics (OPA)
   !!======================================================================
   !! History :  2.0  !  2002-06  (L. Debreu)  Original cade
   !!            3.2  !  2009-04  (R. Benshila) 
   !!            3.6  !  2014-09  (R. Benshila) 
   !!----------------------------------------------------------------------
#if defined key_agrif
   !!----------------------------------------------------------------------
   !!   'key_agrif'                                              AGRIF zoom
   !!----------------------------------------------------------------------
   !!   Agrif_tra     :
   !!   Agrif_dyn     : 
   !!   Agrif_ssh     :
   !!   Agrif_dyn_ts  :
   !!   Agrif_dta_ts  :
   !!   Agrif_ssh_ts  :
   !!   Agrif_avm     : 
   !!   interpu       :
   !!   interpv       :
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce
   USE dom_oce      
   USE zdf_oce
   USE agrif_oce
   USE phycst
   USE dynspg_ts, ONLY: un_adv, vn_adv
   !
   USE in_out_manager
   USE agrif_oce_sponge
   USE lib_mpp
   USE vremap
 
   IMPLICIT NONE
   PRIVATE

   PUBLIC   Agrif_dyn, Agrif_ssh, Agrif_dyn_ts, Agrif_dyn_ts_flux, Agrif_ssh_ts, Agrif_dta_ts
   PUBLIC   Agrif_tra, Agrif_avm
   PUBLIC   interpun , interpvn
   PUBLIC   interptsn, interpsshn, interpavm
   PUBLIC   interpunb, interpvnb , interpub2b, interpvb2b
   PUBLIC   interpe3t
#if defined key_vertical
   PUBLIC   interpht0, interpmbkt
# endif
   INTEGER ::   bdy_tinterp = 0

   !!----------------------------------------------------------------------
   !! NEMO/NST 4.0 , NEMO Consortium (2018)
   !! $Id: agrif_oce_interp.F90 12377 2020-02-12 14:39:06Z acc $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE Agrif_tra
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_tra  ***
      !!----------------------------------------------------------------------
      !
      IF( Agrif_Root() )   RETURN
      !
      Agrif_SpecialValue    = 0._wp
      Agrif_UseSpecialValue = .TRUE.
      !
      CALL Agrif_Bc_variable( tsn_id, procname=interptsn )
      !
      Agrif_UseSpecialValue = .FALSE.
      !
   END SUBROUTINE Agrif_tra


   SUBROUTINE Agrif_dyn( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_DYN  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) ::   kt
      !
      INTEGER ::   ji, jj, jk       ! dummy loop indices
      INTEGER ::   ibdy1, jbdy1, ibdy2, jbdy2
      REAL(wp), DIMENSION(jpi,jpj) ::   zub, zvb
      !!----------------------------------------------------------------------  
      !
      IF( Agrif_Root() )   RETURN
      !
      Agrif_SpecialValue    = 0._wp
      Agrif_UseSpecialValue = ln_spc_dyn
      !
      CALL Agrif_Bc_variable( un_interp_id, procname=interpun )
      CALL Agrif_Bc_variable( vn_interp_id, procname=interpvn )
      !
      Agrif_UseSpecialValue = .FALSE.
      !
      ! --- West --- !
      ibdy1 = 2
      ibdy2 = 1+nbghostcells 
      !
      IF( .NOT.ln_dynspg_ts ) THEN  ! Store transport
         DO ji = mi0(ibdy1), mi1(ibdy2)
            uu_b(ji,:,Krhs_a) = 0._wp

            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  uu_b(ji,jj,Krhs_a) = uu_b(ji,jj,Krhs_a) + e3u(ji,jj,jk,Krhs_a) * uu(ji,jj,jk,Krhs_a) * umask(ji,jj,jk)
               END DO
            END DO

            DO jj = 1, jpj
               uu_b(ji,jj,Krhs_a) = uu_b(ji,jj,Krhs_a) * r1_hu(ji,jj,Krhs_a)
            END DO
         END DO
      ENDIF
      !
      DO ji = mi0(ibdy1), mi1(ibdy2)
         zub(ji,:) = 0._wp    ! Correct transport
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               zub(ji,jj) = zub(ji,jj) & 
                  & + e3u(ji,jj,jk,Krhs_a)  * uu(ji,jj,jk,Krhs_a)*umask(ji,jj,jk)
            END DO
         END DO
         DO jj=1,jpj
            zub(ji,jj) = zub(ji,jj) * r1_hu(ji,jj,Krhs_a)
         END DO
            
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               uu(ji,jj,jk,Krhs_a) = ( uu(ji,jj,jk,Krhs_a) + uu_b(ji,jj,Krhs_a)-zub(ji,jj)) * umask(ji,jj,jk)
            END DO
         END DO
      END DO
            
      IF( ln_dynspg_ts ) THEN       ! Set tangential velocities to time splitting estimate
         DO ji = mi0(ibdy1), mi1(ibdy2)
            zvb(ji,:) = 0._wp
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  zvb(ji,jj) = zvb(ji,jj) + e3v(ji,jj,jk,Krhs_a) * vv(ji,jj,jk,Krhs_a) * vmask(ji,jj,jk)
               END DO
            END DO
            DO jj = 1, jpj
               zvb(ji,jj) = zvb(ji,jj) * r1_hv(ji,jj,Krhs_a)
            END DO
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  vv(ji,jj,jk,Krhs_a) = ( vv(ji,jj,jk,Krhs_a) + vv_b(ji,jj,Krhs_a)-zvb(ji,jj))*vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      ! --- East --- !
      ibdy1 = jpiglo-1-nbghostcells
      ibdy2 = jpiglo-2 
      !
      IF( .NOT.ln_dynspg_ts ) THEN  ! Store transport
         DO ji = mi0(ibdy1), mi1(ibdy2)
            uu_b(ji,:,Krhs_a) = 0._wp
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  uu_b(ji,jj,Krhs_a) = uu_b(ji,jj,Krhs_a) & 
                      & + e3u(ji,jj,jk,Krhs_a) * uu(ji,jj,jk,Krhs_a) * umask(ji,jj,jk)
               END DO
            END DO
            DO jj = 1, jpj
               uu_b(ji,jj,Krhs_a) = uu_b(ji,jj,Krhs_a) * r1_hu(ji,jj,Krhs_a)
            END DO
         END DO
      ENDIF
      !
      DO ji = mi0(ibdy1), mi1(ibdy2)
         zub(ji,:) = 0._wp    ! Correct transport
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               zub(ji,jj) = zub(ji,jj) & 
                  & + e3u(ji,jj,jk,Krhs_a)  * uu(ji,jj,jk,Krhs_a) * umask(ji,jj,jk)
            END DO
         END DO
         DO jj=1,jpj
            zub(ji,jj) = zub(ji,jj) * r1_hu(ji,jj,Krhs_a)
         END DO
            
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               uu(ji,jj,jk,Krhs_a) = ( uu(ji,jj,jk,Krhs_a) & 
                 & + uu_b(ji,jj,Krhs_a)-zub(ji,jj))*umask(ji,jj,jk)
            END DO
         END DO
      END DO
            
      IF( ln_dynspg_ts ) THEN       ! Set tangential velocities to time splitting estimate
         ibdy1 = jpiglo-nbghostcells
         ibdy2 = jpiglo-1 
         DO ji = mi0(ibdy1), mi1(ibdy2)
            zvb(ji,:) = 0._wp
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  zvb(ji,jj) = zvb(ji,jj) &
                     & + e3v(ji,jj,jk,Krhs_a) * vv(ji,jj,jk,Krhs_a) * vmask(ji,jj,jk)
               END DO
            END DO
            DO jj = 1, jpj
               zvb(ji,jj) = zvb(ji,jj) * r1_hv(ji,jj,Krhs_a)
            END DO
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  vv(ji,jj,jk,Krhs_a) = ( vv(ji,jj,jk,Krhs_a) & 
                      & + vv_b(ji,jj,Krhs_a)-zvb(ji,jj)) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      ! --- South --- !
      jbdy1 = 2
      jbdy2 = 1+nbghostcells 
      !
      IF( .NOT.ln_dynspg_ts ) THEN  ! Store transport
         DO jj = mj0(jbdy1), mj1(jbdy2)
            vv_b(:,jj,Krhs_a) = 0._wp
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  vv_b(ji,jj,Krhs_a) = vv_b(ji,jj,Krhs_a) & 
                      & + e3v(ji,jj,jk,Krhs_a) * vv(ji,jj,jk,Krhs_a) * vmask(ji,jj,jk)
               END DO
            END DO
            DO ji=1,jpi
               vv_b(ji,jj,Krhs_a) = vv_b(ji,jj,Krhs_a) * r1_hv(ji,jj,Krhs_a)     
            END DO
         END DO
      ENDIF
      !
      DO jj = mj0(jbdy1), mj1(jbdy2)
         zvb(:,jj) = 0._wp    ! Correct transport
         DO jk=1,jpkm1
            DO ji=1,jpi
               zvb(ji,jj) = zvb(ji,jj) & 
                  & + e3v(ji,jj,jk,Krhs_a) * vv(ji,jj,jk,Krhs_a) * vmask(ji,jj,jk)
            END DO
         END DO
         DO ji = 1, jpi
            zvb(ji,jj) = zvb(ji,jj) * r1_hv(ji,jj,Krhs_a)
         END DO

         DO jk = 1, jpkm1
            DO ji = 1, jpi
               vv(ji,jj,jk,Krhs_a) = ( vv(ji,jj,jk,Krhs_a) & 
                 & + vv_b(ji,jj,Krhs_a) - zvb(ji,jj) ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
            
      IF( ln_dynspg_ts ) THEN       ! Set tangential velocities to time splitting estimate
         DO jj = mj0(jbdy1), mj1(jbdy2)
            zub(:,jj) = 0._wp
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  zub(ji,jj) = zub(ji,jj) & 
                     & + e3u(ji,jj,jk,Krhs_a) * uu(ji,jj,jk,Krhs_a) * umask(ji,jj,jk)
               END DO
            END DO
            DO ji = 1, jpi
               zub(ji,jj) = zub(ji,jj) * r1_hu(ji,jj,Krhs_a)
            END DO
               
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  uu(ji,jj,jk,Krhs_a) = ( uu(ji,jj,jk,Krhs_a) & 
                    & + uu_b(ji,jj,Krhs_a) - zub(ji,jj) ) * umask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      ! --- North --- !
      jbdy1 = jpjglo-1-nbghostcells
      jbdy2 = jpjglo-2 
      !
      IF( .NOT.ln_dynspg_ts ) THEN  ! Store transport
         DO jj = mj0(jbdy1), mj1(jbdy2)
            vv_b(:,jj,Krhs_a) = 0._wp
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  vv_b(ji,jj,Krhs_a) = vv_b(ji,jj,Krhs_a) & 
                      & + e3v(ji,jj,jk,Krhs_a) * vv(ji,jj,jk,Krhs_a) * vmask(ji,jj,jk)
               END DO
            END DO
            DO ji=1,jpi
               vv_b(ji,jj,Krhs_a) = vv_b(ji,jj,Krhs_a) * r1_hv(ji,jj,Krhs_a)
            END DO
         END DO
      ENDIF
      !
      DO jj = mj0(jbdy1), mj1(jbdy2)
         zvb(:,jj) = 0._wp    ! Correct transport
         DO jk=1,jpkm1
            DO ji=1,jpi
               zvb(ji,jj) = zvb(ji,jj) & 
                  & + e3v(ji,jj,jk,Krhs_a) * vv(ji,jj,jk,Krhs_a) * vmask(ji,jj,jk)
            END DO
         END DO
         DO ji = 1, jpi
            zvb(ji,jj) = zvb(ji,jj) * r1_hv(ji,jj,Krhs_a)
         END DO

         DO jk = 1, jpkm1
            DO ji = 1, jpi
               vv(ji,jj,jk,Krhs_a) = ( vv(ji,jj,jk,Krhs_a) & 
                 & + vv_b(ji,jj,Krhs_a) - zvb(ji,jj) ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
            
      IF( ln_dynspg_ts ) THEN       ! Set tangential velocities to time splitting estimate
         jbdy1 = jpjglo-nbghostcells
         jbdy2 = jpjglo-1
         DO jj = mj0(jbdy1), mj1(jbdy2)
            zub(:,jj) = 0._wp
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  zub(ji,jj) = zub(ji,jj) & 
                     & + e3u(ji,jj,jk,Krhs_a) * uu(ji,jj,jk,Krhs_a) * umask(ji,jj,jk)
               END DO
            END DO
            DO ji = 1, jpi
               zub(ji,jj) = zub(ji,jj) * r1_hu(ji,jj,Krhs_a)
            END DO
               
            DO jk = 1, jpkm1
               DO ji = 1, jpi
                  uu(ji,jj,jk,Krhs_a) = ( uu(ji,jj,jk,Krhs_a) & 
                    & + uu_b(ji,jj,Krhs_a) - zub(ji,jj) ) * umask(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE Agrif_dyn


   SUBROUTINE Agrif_dyn_ts( jn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_dyn_ts  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) ::   jn
      !!
      INTEGER :: ji, jj
      INTEGER :: istart, iend, jstart, jend
      !!----------------------------------------------------------------------  
      !
      IF( Agrif_Root() )   RETURN
      !
      !--- West ---!
      istart = 2
      iend   = nbghostcells+1
      DO ji = mi0(istart), mi1(iend)
         DO jj=1,jpj
            va_e(ji,jj) = vbdy(ji,jj) * hvr_e(ji,jj)
            ua_e(ji,jj) = ubdy(ji,jj) * hur_e(ji,jj)
         END DO
      END DO
      !
      !--- East ---!
      istart = jpiglo-nbghostcells
      iend   = jpiglo-1
      DO ji = mi0(istart), mi1(iend)
         DO jj=1,jpj
            va_e(ji,jj) = vbdy(ji,jj) * hvr_e(ji,jj)
         END DO
      END DO
      istart = jpiglo-nbghostcells-1
      iend   = jpiglo-2
      DO ji = mi0(istart), mi1(iend)
         DO jj=1,jpj
            ua_e(ji,jj) = ubdy(ji,jj) * hur_e(ji,jj)
         END DO
      END DO
      !
      !--- South ---!
      jstart = 2
      jend   = nbghostcells+1
      DO jj = mj0(jstart), mj1(jend)
         DO ji=1,jpi
            ua_e(ji,jj) = ubdy(ji,jj) * hur_e(ji,jj)
            va_e(ji,jj) = vbdy(ji,jj) * hvr_e(ji,jj)
         END DO
      END DO
      !
      !--- North ---!
      jstart = jpjglo-nbghostcells
      jend   = jpjglo-1
      DO jj = mj0(jstart), mj1(jend)
         DO ji=1,jpi
            ua_e(ji,jj) = ubdy(ji,jj) * hur_e(ji,jj)
         END DO
      END DO
      jstart = jpjglo-nbghostcells-1
      jend   = jpjglo-2
      DO jj = mj0(jstart), mj1(jend)
         DO ji=1,jpi
            va_e(ji,jj) = vbdy(ji,jj) * hvr_e(ji,jj)
         END DO
      END DO
      !
   END SUBROUTINE Agrif_dyn_ts

   SUBROUTINE Agrif_dyn_ts_flux( jn, zu, zv )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_dyn_ts_flux  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) ::   jn
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   zu, zv
      !!
      INTEGER :: ji, jj
      INTEGER :: istart, iend, jstart, jend
      !!----------------------------------------------------------------------  
      !
      IF( Agrif_Root() )   RETURN
      !
      !--- West ---!
      istart = 2
      iend   = nbghostcells+1
      DO ji = mi0(istart), mi1(iend)
         DO jj=1,jpj
            zv(ji,jj) = vbdy(ji,jj) * e1v(ji,jj)
            zu(ji,jj) = ubdy(ji,jj) * e2u(ji,jj)
         END DO
      END DO
      !
      !--- East ---!
      istart = jpiglo-nbghostcells
      iend   = jpiglo-1
      DO ji = mi0(istart), mi1(iend)
         DO jj=1,jpj
            zv(ji,jj) = vbdy(ji,jj) * e1v(ji,jj)
         END DO
      END DO
      istart = jpiglo-nbghostcells-1
      iend   = jpiglo-2
      DO ji = mi0(istart), mi1(iend)
         DO jj=1,jpj
            zu(ji,jj) = ubdy(ji,jj) * e2u(ji,jj)
         END DO
      END DO
      !
      !--- South ---!
      jstart = 2
      jend   = nbghostcells+1
      DO jj = mj0(jstart), mj1(jend)
         DO ji=1,jpi
            zu(ji,jj) = ubdy(ji,jj) * e2u(ji,jj)
            zv(ji,jj) = vbdy(ji,jj) * e1v(ji,jj)
         END DO
      END DO
      !
      !--- North ---!
      jstart = jpjglo-nbghostcells
      jend   = jpjglo-1
      DO jj = mj0(jstart), mj1(jend)
         DO ji=1,jpi
            zu(ji,jj) = ubdy(ji,jj) * e2u(ji,jj)
         END DO
      END DO
      jstart = jpjglo-nbghostcells-1
      jend   = jpjglo-2
      DO jj = mj0(jstart), mj1(jend)
         DO ji=1,jpi
            zv(ji,jj) = vbdy(ji,jj) * e1v(ji,jj)
         END DO
      END DO
      !
   END SUBROUTINE Agrif_dyn_ts_flux

   SUBROUTINE Agrif_dta_ts( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_dta_ts  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) ::   kt
      !!
      INTEGER :: ji, jj
      LOGICAL :: ll_int_cons
      !!----------------------------------------------------------------------  
      !
      IF( Agrif_Root() )   RETURN
      !
      ll_int_cons = ln_bt_fw ! Assume conservative temporal integration in the forward case only
      !
      ! Enforce volume conservation if no time refinement:  
      IF ( Agrif_rhot()==1 ) ll_int_cons=.TRUE.  
      !
      ! Interpolate barotropic fluxes
      Agrif_SpecialValue = 0._wp
      Agrif_UseSpecialValue = ln_spc_dyn
      !
      ! Set bdy time interpolation stage to 0 (latter incremented locally do deal with corners)
      utint_stage(:,:) = 0
      vtint_stage(:,:) = 0
      !
      IF( ll_int_cons ) THEN  ! Conservative interpolation
         ! order matters here !!!!!!
         CALL Agrif_Bc_variable( ub2b_interp_id, calledweight=1._wp, procname=interpub2b ) ! Time integrated
         CALL Agrif_Bc_variable( vb2b_interp_id, calledweight=1._wp, procname=interpvb2b ) 
         !
         bdy_tinterp = 1
         CALL Agrif_Bc_variable( unb_id        , calledweight=1._wp, procname=interpunb  ) ! After
         CALL Agrif_Bc_variable( vnb_id        , calledweight=1._wp, procname=interpvnb  )  
         !
         bdy_tinterp = 2
         CALL Agrif_Bc_variable( unb_id        , calledweight=0._wp, procname=interpunb  ) ! Before
         CALL Agrif_Bc_variable( vnb_id        , calledweight=0._wp, procname=interpvnb  )   
      ELSE ! Linear interpolation
         !
         ubdy(:,:) = 0._wp   ;   vbdy(:,:) = 0._wp 
         CALL Agrif_Bc_variable( unb_id, procname=interpunb )
         CALL Agrif_Bc_variable( vnb_id, procname=interpvnb )
      ENDIF
      Agrif_UseSpecialValue = .FALSE.
      ! 
   END SUBROUTINE Agrif_dta_ts


   SUBROUTINE Agrif_ssh( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_ssh  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) ::   kt
      !
      INTEGER  :: ji, jj
      INTEGER  :: istart, iend, jstart, jend
      !!----------------------------------------------------------------------  
      !
      IF( Agrif_Root() )   RETURN
      !      
      ! Linear time interpolation of sea level
      !
      Agrif_SpecialValue    = 0._wp
      Agrif_UseSpecialValue = .TRUE.
      CALL Agrif_Bc_variable(sshn_id, procname=interpsshn )
      Agrif_UseSpecialValue = .FALSE.
      !
      ! --- West --- !
      istart = 2
      iend   = 1 + nbghostcells
      DO ji = mi0(istart), mi1(iend)
         DO jj = 1, jpj
            ssh(ji,jj,Krhs_a) = hbdy(ji,jj)
         ENDDO
      ENDDO
      !
      ! --- East --- !
      istart = jpiglo - nbghostcells
      iend   = jpiglo - 1
      DO ji = mi0(istart), mi1(iend)
         DO jj = 1, jpj
            ssh(ji,jj,Krhs_a) = hbdy(ji,jj)
         ENDDO
      ENDDO
      !
      ! --- South --- !
      jstart = 2
      jend   = 1 + nbghostcells
      DO jj = mj0(jstart), mj1(jend)
         DO ji = 1, jpi
            ssh(ji,jj,Krhs_a) = hbdy(ji,jj)
         ENDDO
      ENDDO
      !
      ! --- North --- !
      jstart = jpjglo - nbghostcells
      jend   = jpjglo - 1
      DO jj = mj0(jstart), mj1(jend)
         DO ji = 1, jpi
            ssh(ji,jj,Krhs_a) = hbdy(ji,jj)
         ENDDO
      ENDDO
      !
   END SUBROUTINE Agrif_ssh


   SUBROUTINE Agrif_ssh_ts( jn )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_ssh_ts  ***
      !!----------------------------------------------------------------------  
      INTEGER, INTENT(in) ::   jn
      !!
      INTEGER :: ji, jj
      INTEGER  :: istart, iend, jstart, jend
      !!----------------------------------------------------------------------  
      !
      IF( Agrif_Root() )   RETURN
      !
      ! --- West --- !
      istart = 2
      iend   = 1+nbghostcells
      DO ji = mi0(istart), mi1(iend)
         DO jj = 1, jpj
            ssha_e(ji,jj) = hbdy(ji,jj)
         ENDDO
      ENDDO
      !
      ! --- East --- !
      istart = jpiglo - nbghostcells
      iend   = jpiglo - 1
      DO ji = mi0(istart), mi1(iend)
         DO jj = 1, jpj
            ssha_e(ji,jj) = hbdy(ji,jj)
         ENDDO
      ENDDO
      !
      ! --- South --- !
      jstart = 2
      jend   = 1+nbghostcells
      DO jj = mj0(jstart), mj1(jend)
         DO ji = 1, jpi
            ssha_e(ji,jj) = hbdy(ji,jj)
         ENDDO
      ENDDO
      !
      ! --- North --- !
      jstart = jpjglo - nbghostcells
      jend   = jpjglo - 1
      DO jj = mj0(jstart), mj1(jend)
         DO ji = 1, jpi
            ssha_e(ji,jj) = hbdy(ji,jj)
         ENDDO
      ENDDO
      !
   END SUBROUTINE Agrif_ssh_ts

   SUBROUTINE Agrif_avm
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE Agrif_avm  ***
      !!----------------------------------------------------------------------  
      REAL(wp) ::   zalpha
      !!----------------------------------------------------------------------  
      !
      IF( Agrif_Root() )   RETURN
      !
      zalpha = 1._wp ! JC: proper time interpolation impossible  
                     ! => use last available value from parent 
      !
      Agrif_SpecialValue    = 0.e0
      Agrif_UseSpecialValue = .TRUE.
      !
      CALL Agrif_Bc_variable( avm_id, calledweight=zalpha, procname=interpavm )       
      !
      Agrif_UseSpecialValue = .FALSE.
      !
   END SUBROUTINE Agrif_avm
   

   SUBROUTINE interptsn( ptab, i1, i2, j1, j2, k1, k2, n1, n2, before )
      !!----------------------------------------------------------------------
      !!                  *** ROUTINE interptsn ***
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,n1:n2), INTENT(inout) ::   ptab
      INTEGER                                     , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2, n1, n2
      LOGICAL                                     , INTENT(in   ) ::   before
      !
      INTEGER  ::   ji, jj, jk, jn  ! dummy loop indices
      INTEGER  ::   N_in, N_out
      ! vertical interpolation:
      REAL(wp) :: zhtot
      REAL(wp), DIMENSION(k1:k2,1:jpts) :: tabin
      REAL(wp), DIMENSION(k1:k2) :: h_in
      REAL(wp), DIMENSION(1:jpk) :: h_out
      !!----------------------------------------------------------------------

      IF( before ) THEN         
         DO jn = 1,jpts
            DO jk=k1,k2
               DO jj=j1,j2
                 DO ji=i1,i2
                       ptab(ji,jj,jk,jn) = ts(ji,jj,jk,jn,Kmm_a)
                 END DO
              END DO
           END DO
        END DO

# if defined key_vertical
        ! Interpolate thicknesses
        ! Warning: these are masked, hence extrapolated prior interpolation.
        DO jk=k1,k2
           DO jj=j1,j2
              DO ji=i1,i2
                  ptab(ji,jj,jk,jpts+1) = tmask(ji,jj,jk) * e3t(ji,jj,jk,Kmm_a)
              END DO
           END DO
        END DO

        ! Extrapolate thicknesses in partial bottom cells:
        ! Set them to Agrif_SpecialValue (0.). Correct bottom thicknesses are retrieved later on
        IF (ln_zps) THEN
           DO jj=j1,j2
              DO ji=i1,i2
                  jk = mbkt(ji,jj)
                  ptab(ji,jj,jk,jpts+1) = 0._wp
              END DO
           END DO           
        END IF
     
        ! Save ssh at last level:
        IF (.NOT.ln_linssh) THEN
           ptab(i1:i2,j1:j2,k2,jpts+1) = ssh(i1:i2,j1:j2,Kmm_a)*tmask(i1:i2,j1:j2,1) 
        ELSE
           ptab(i1:i2,j1:j2,k2,jpts+1) = 0._wp
        END IF      
# endif
      ELSE 

# if defined key_vertical 
         IF (ln_linssh) ptab(i1:i2,j1:j2,k2,n2) = 0._wp 
            
         DO jj=j1,j2
            DO ji=i1,i2
               ts(ji,jj,:,:,Krhs_a) = 0._wp
               N_in = mbkt_parent(ji,jj)
               zhtot = 0._wp
               DO jk=1,N_in !k2 = jpk of parent grid
                  IF (jk==N_in) THEN
                     h_in(jk) = ht0_parent(ji,jj) + ptab(ji,jj,k2,n2) - zhtot
                  ELSE
                     h_in(jk) = ptab(ji,jj,jk,n2)
                  ENDIF
                  zhtot = zhtot + h_in(jk)
                  tabin(jk,:) = ptab(ji,jj,jk,n1:n2-1)
               END DO
               N_out = 0
               DO jk=1,jpk ! jpk of child grid
                  IF (tmask(ji,jj,jk) == 0._wp) EXIT 
                  N_out = N_out + 1
                  h_out(jk) = e3t(ji,jj,jk,Krhs_a)
               ENDDO
               IF (N_in*N_out > 0) THEN
                  CALL reconstructandremap(tabin(1:N_in,1:jpts),h_in(1:N_in),ts(ji,jj,1:N_out,1:jpts,Krhs_a),h_out(1:N_out),N_in,N_out,jpts)
               ENDIF
            ENDDO
         ENDDO
# else
         !
         DO jn=1, jpts
            ts(i1:i2,j1:j2,1:jpk,jn,Krhs_a)=ptab(i1:i2,j1:j2,1:jpk,jn)*tmask(i1:i2,j1:j2,1:jpk) 
         END DO
# endif

      ENDIF
      !
   END SUBROUTINE interptsn

   SUBROUTINE interpsshn( ptab, i1, i2, j1, j2, before )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpsshn  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      !
      !!----------------------------------------------------------------------  
      !
      IF( before) THEN
         ptab(i1:i2,j1:j2) = ssh(i1:i2,j1:j2,Kmm_a)
      ELSE
         hbdy(i1:i2,j1:j2) = ptab(i1:i2,j1:j2) * tmask(i1:i2,j1:j2,1)
      ENDIF
      !
   END SUBROUTINE interpsshn

   SUBROUTINE interpun( ptab, i1, i2, j1, j2, k1, k2, m1, m2, before )
      !!----------------------------------------------------------------------
      !!                  *** ROUTINE interpun ***
      !!---------------------------------------------    
      !!
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,m1,m2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,m1:m2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      !!
      INTEGER :: ji,jj,jk
      REAL(wp) :: zrhoy, zhtot
      ! vertical interpolation:
      REAL(wp), DIMENSION(k1:k2) :: tabin, h_in
      REAL(wp), DIMENSION(1:jpk) :: h_out
      INTEGER  :: N_in, N_out
      REAL(wp) :: h_diff
      !!---------------------------------------------    
      !
      IF (before) THEN 
         DO jk=1,jpk
            DO jj=j1,j2
               DO ji=i1,i2
                  ptab(ji,jj,jk,1) = (e2u(ji,jj) * e3u(ji,jj,jk,Kmm_a) * uu(ji,jj,jk,Kmm_a)*umask(ji,jj,jk)) 
# if defined key_vertical
                  ! Interpolate thicknesses (masked for subsequent extrapolation)
                  ptab(ji,jj,jk,2) = umask(ji,jj,jk) * e2u(ji,jj) * e3u(ji,jj,jk,Kmm_a)
# endif
               END DO
            END DO
         END DO
# if defined key_vertical
         ! Extrapolate thicknesses in partial bottom cells:
         ! Set them to Agrif_SpecialValue (0.). Correct bottom thicknesses are retrieved later on
         IF (ln_zps) THEN
            DO jj=j1,j2
               DO ji=i1,i2
                  jk = mbku(ji,jj)
                  ptab(ji,jj,jk,2) = 0._wp
               END DO
            END DO           
         END IF
        ! Save ssh at last level:
        ptab(i1:i2,j1:j2,k2,2) = 0._wp
        IF (.NOT.ln_linssh) THEN
           ! This vertical sum below should be replaced by the sea-level at U-points (optimization):
           DO jk=1,jpk
              ptab(i1:i2,j1:j2,k2,2) = ptab(i1:i2,j1:j2,k2,2) + e3u(i1:i2,j1:j2,jk,Kmm_a) * umask(i1:i2,j1:j2,jk)
           END DO
           ptab(i1:i2,j1:j2,k2,2) = ptab(i1:i2,j1:j2,k2,2) - hu_0(i1:i2,j1:j2)
        END IF 
# endif
         !
      ELSE
         zrhoy = Agrif_rhoy()
# if defined key_vertical
! VERTICAL REFINEMENT BEGIN

         IF (ln_linssh) ptab(i1:i2,j1:j2,k2,2) = 0._wp 

         DO ji=i1,i2
            DO jj=j1,j2
               uu(ji,jj,:,Krhs_a) = 0._wp
               N_in = mbku_parent(ji,jj)
               zhtot = 0._wp
               DO jk=1,N_in
                  IF (jk==N_in) THEN
                     h_in(jk) = hu0_parent(ji,jj) + ptab(ji,jj,k2,2) - zhtot
                  ELSE
                     h_in(jk) = ptab(ji,jj,jk,2)/(e2u(ji,jj)*zrhoy) 
                  ENDIF
                  zhtot = zhtot + h_in(jk)
                  tabin(jk) = ptab(ji,jj,jk,1)/(e2u(ji,jj)*zrhoy*h_in(jk))
              ENDDO
                  
              N_out = 0
              DO jk=1,jpk
                 if (umask(ji,jj,jk) == 0) EXIT
                 N_out = N_out + 1
                 h_out(N_out) = e3u(ji,jj,jk,Krhs_a)
              ENDDO
              IF (N_in*N_out > 0) THEN
                 CALL reconstructandremap(tabin(1:N_in),h_in(1:N_in),uu(ji,jj,1:N_out,Krhs_a),h_out(1:N_out),N_in,N_out,1)
              ENDIF
            ENDDO
         ENDDO

# else
         DO jk = 1, jpkm1
            DO jj=j1,j2
               uu(i1:i2,jj,jk,Krhs_a) = ptab(i1:i2,jj,jk,1) / ( zrhoy * e2u(i1:i2,jj) * e3u(i1:i2,jj,jk,Krhs_a) )
            END DO
         END DO
# endif

      ENDIF
      ! 
   END SUBROUTINE interpun

   SUBROUTINE interpvn( ptab, i1, i2, j1, j2, k1, k2, m1, m2, before )
      !!----------------------------------------------------------------------
      !!                  *** ROUTINE interpvn ***
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: i1,i2,j1,j2,k1,k2,m1,m2
      REAL(wp), DIMENSION(i1:i2,j1:j2,k1:k2,m1:m2), INTENT(inout) :: ptab
      LOGICAL, INTENT(in) :: before
      !
      INTEGER :: ji,jj,jk
      REAL(wp) :: zrhox
      ! vertical interpolation:
      REAL(wp), DIMENSION(k1:k2) :: tabin, h_in
      REAL(wp), DIMENSION(1:jpk) :: h_out
      INTEGER  :: N_in, N_out
      REAL(wp) :: h_diff, zhtot
      !!---------------------------------------------    
      !      
      IF (before) THEN          
         DO jk=k1,k2
            DO jj=j1,j2
               DO ji=i1,i2
                  ptab(ji,jj,jk,1) = (e1v(ji,jj) * e3v(ji,jj,jk,Kmm_a) * vv(ji,jj,jk,Kmm_a)*vmask(ji,jj,jk))
# if defined key_vertical
                  ! Interpolate thicknesses (masked for subsequent extrapolation)
                  ptab(ji,jj,jk,2) = vmask(ji,jj,jk) * e1v(ji,jj) * e3v(ji,jj,jk,Kmm_a)
# endif
               END DO
            END DO
         END DO
# if defined key_vertical
         ! Extrapolate thicknesses in partial bottom cells:
         ! Set them to Agrif_SpecialValue (0.). Correct bottom thicknesses are retrieved later on
         IF (ln_zps) THEN
            DO jj=j1,j2
               DO ji=i1,i2
                  jk = mbkv(ji,jj)
                  ptab(ji,jj,jk,2) = 0._wp
               END DO
            END DO           
         END IF
        ! Save ssh at last level:
        ptab(i1:i2,j1:j2,k2,2) = 0._wp
        IF (.NOT.ln_linssh) THEN
           ! This vertical sum below should be replaced by the sea-level at V-points (optimization):
           DO jk=1,jpk
              ptab(i1:i2,j1:j2,k2,2) = ptab(i1:i2,j1:j2,k2,2) + e3v(i1:i2,j1:j2,jk,Kmm_a) * vmask(i1:i2,j1:j2,jk)
           END DO
           ptab(i1:i2,j1:j2,k2,2) = ptab(i1:i2,j1:j2,k2,2) - hv_0(i1:i2,j1:j2)
        END IF 
# endif
      ELSE       
         zrhox = Agrif_rhox()
# if defined key_vertical

         IF (ln_linssh) ptab(i1:i2,j1:j2,k2,2) = 0._wp 

         DO jj=j1,j2
            DO ji=i1,i2
               vv(ji,jj,:,Krhs_a) = 0._wp
               N_in = mbkv_parent(ji,jj)
               zhtot = 0._wp
               DO jk=1,N_in
                  IF (jk==N_in) THEN
                     h_in(jk) = hv0_parent(ji,jj) + ptab(ji,jj,k2,2) - zhtot
                  ELSE
                     h_in(jk) = ptab(ji,jj,jk,2)/(e1v(ji,jj)*zrhox) 
                  ENDIF
                  zhtot = zhtot + h_in(jk)
                  tabin(jk) = ptab(ji,jj,jk,1)/(e1v(ji,jj)*zrhox*h_in(jk))
              ENDDO
         
               N_out = 0
               DO jk=1,jpk
                  if (vmask(ji,jj,jk) == 0) EXIT
                  N_out = N_out + 1
                  h_out(N_out) = e3v(ji,jj,jk,Krhs_a)
               END DO
               IF (N_in*N_out > 0) THEN
                  call reconstructandremap(tabin(1:N_in),h_in(1:N_in),vv(ji,jj,1:N_out,Krhs_a),h_out(1:N_out),N_in,N_out,1)
               ENDIF
            END DO
         END DO
# else
         DO jk = 1, jpkm1
            vv(i1:i2,j1:j2,jk,Krhs_a) = ptab(i1:i2,j1:j2,jk,1) / ( zrhox * e1v(i1:i2,j1:j2) * e3v(i1:i2,j1:j2,jk,Krhs_a) )
         END DO
# endif
      ENDIF
      !        
   END SUBROUTINE interpvn

   SUBROUTINE interpunb( ptab, i1, i2, j1, j2, before)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpunb  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      !
      INTEGER  ::   ji, jj
      REAL(wp) ::   zrhoy, zrhot, zt0, zt1, ztcoeff
      !!----------------------------------------------------------------------  
      !
      IF( before ) THEN 
         ptab(i1:i2,j1:j2) = e2u(i1:i2,j1:j2) * hu(i1:i2,j1:j2,Kmm_a) * uu_b(i1:i2,j1:j2,Kmm_a)
      ELSE
         zrhoy = Agrif_Rhoy()
         zrhot = Agrif_rhot()
         ! Time indexes bounds for integration
         zt0 = REAL(Agrif_NbStepint()  , wp) / zrhot
         zt1 = REAL(Agrif_NbStepint()+1, wp) / zrhot      
         ! 
         DO ji = i1, i2
            DO jj = j1, j2
               IF ( utint_stage(ji,jj) < (bdy_tinterp + 1) ) THEN
                  IF    ( utint_stage(ji,jj) == 1  ) THEN
                     ztcoeff = zrhot * (  zt1**2._wp * (       zt1 - 1._wp)        &
                        &               - zt0**2._wp * (       zt0 - 1._wp)        )
                  ELSEIF( utint_stage(ji,jj) == 2  ) THEN
                     ztcoeff = zrhot * (  zt1        * (       zt1 - 1._wp)**2._wp &
                        &               - zt0        * (       zt0 - 1._wp)**2._wp )
                  ELSEIF( utint_stage(ji,jj) == 0  ) THEN                
                     ztcoeff = 1._wp
                  ELSE
                     ztcoeff = 0._wp
                  ENDIF
                  !   
                  ubdy(ji,jj) = ubdy(ji,jj) + ztcoeff * ptab(ji,jj)
                  !            
                  IF (( utint_stage(ji,jj) == 2 ).OR.( utint_stage(ji,jj) == 0 )) THEN
                     ubdy(ji,jj) = ubdy(ji,jj) / (zrhoy*e2u(ji,jj)) * umask(ji,jj,1)
                  ENDIF
                  !
                  utint_stage(ji,jj) = utint_stage(ji,jj) + 1
               ENDIF
            END DO
         END DO
      END IF
      ! 
   END SUBROUTINE interpunb


   SUBROUTINE interpvnb( ptab, i1, i2, j1, j2, before )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpvnb  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      !
      INTEGER  ::   ji, jj
      REAL(wp) ::   zrhox, zrhot, zt0, zt1, ztcoeff   
      !!----------------------------------------------------------------------  
      ! 
      IF( before ) THEN 
         ptab(i1:i2,j1:j2) = e1v(i1:i2,j1:j2) * hv(i1:i2,j1:j2,Kmm_a) * vv_b(i1:i2,j1:j2,Kmm_a)
      ELSE
         zrhox = Agrif_Rhox()
         zrhot = Agrif_rhot()
         ! Time indexes bounds for integration
         zt0 = REAL(Agrif_NbStepint()  , wp) / zrhot
         zt1 = REAL(Agrif_NbStepint()+1, wp) / zrhot 
         !     
         DO ji = i1, i2
            DO jj = j1, j2
               IF ( vtint_stage(ji,jj) < (bdy_tinterp + 1) ) THEN
                  IF    ( vtint_stage(ji,jj) == 1  ) THEN
                     ztcoeff = zrhot * (  zt1**2._wp * (       zt1 - 1._wp)        &
                        &               - zt0**2._wp * (       zt0 - 1._wp)        )
                  ELSEIF( vtint_stage(ji,jj) == 2  ) THEN
                     ztcoeff = zrhot * (  zt1        * (       zt1 - 1._wp)**2._wp &
                        &               - zt0        * (       zt0 - 1._wp)**2._wp )
                  ELSEIF( vtint_stage(ji,jj) == 0  ) THEN                
                     ztcoeff = 1._wp
                  ELSE
                     ztcoeff = 0._wp
                  ENDIF
                  !   
                  vbdy(ji,jj) = vbdy(ji,jj) + ztcoeff * ptab(ji,jj)
                  !            
                  IF (( vtint_stage(ji,jj) == 2 ).OR.( vtint_stage(ji,jj) == 0 )) THEN
                     vbdy(ji,jj) = vbdy(ji,jj) / (zrhox*e1v(ji,jj)) * vmask(ji,jj,1)
                  ENDIF
                  !
                  vtint_stage(ji,jj) = vtint_stage(ji,jj) + 1
               ENDIF
            END DO
         END DO          
      ENDIF
      !
   END SUBROUTINE interpvnb


   SUBROUTINE interpub2b( ptab, i1, i2, j1, j2, before )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpub2b  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      !
      INTEGER  ::   ji,jj
      REAL(wp) ::   zrhot, zt0, zt1, zat
      !!----------------------------------------------------------------------  
      IF( before ) THEN
         IF ( ln_bt_fw ) THEN
            ptab(i1:i2,j1:j2) = e2u(i1:i2,j1:j2) * ub2_b(i1:i2,j1:j2)
         ELSE
            ptab(i1:i2,j1:j2) = e2u(i1:i2,j1:j2) * un_adv(i1:i2,j1:j2)
         ENDIF
      ELSE
         zrhot = Agrif_rhot()
         ! Time indexes bounds for integration
         zt0 = REAL(Agrif_NbStepint()  , wp) / zrhot
         zt1 = REAL(Agrif_NbStepint()+1, wp) / zrhot
         ! Polynomial interpolation coefficients:
         zat = zrhot * (  zt1**2._wp * (-2._wp*zt1 + 3._wp)    &
            &           - zt0**2._wp * (-2._wp*zt0 + 3._wp)    ) 
         !
         ubdy(i1:i2,j1:j2) = zat * ptab(i1:i2,j1:j2) 
         !
         ! Update interpolation stage:
         utint_stage(i1:i2,j1:j2) = 1
      ENDIF
      ! 
   END SUBROUTINE interpub2b
   

   SUBROUTINE interpvb2b( ptab, i1, i2, j1, j2, before )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpvb2b  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      !
      INTEGER ::   ji,jj
      REAL(wp) ::   zrhot, zt0, zt1, zat
      !!----------------------------------------------------------------------  
      !
      IF( before ) THEN
         IF ( ln_bt_fw ) THEN
            ptab(i1:i2,j1:j2) = e1v(i1:i2,j1:j2) * vb2_b(i1:i2,j1:j2)
         ELSE
            ptab(i1:i2,j1:j2) = e1v(i1:i2,j1:j2) * vn_adv(i1:i2,j1:j2)
         ENDIF
      ELSE      
         zrhot = Agrif_rhot()
         ! Time indexes bounds for integration
         zt0 = REAL(Agrif_NbStepint()  , wp) / zrhot
         zt1 = REAL(Agrif_NbStepint()+1, wp) / zrhot
         ! Polynomial interpolation coefficients:
         zat = zrhot * (  zt1**2._wp * (-2._wp*zt1 + 3._wp)    &
            &           - zt0**2._wp * (-2._wp*zt0 + 3._wp)    ) 
         !
         vbdy(i1:i2,j1:j2) = zat * ptab(i1:i2,j1:j2)
         !
         ! update interpolation stage:
         vtint_stage(i1:i2,j1:j2) = 1
      ENDIF
      !      
   END SUBROUTINE interpvb2b


   SUBROUTINE interpe3t( ptab, i1, i2, j1, j2, k1, k2, before )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpe3t  ***
      !!----------------------------------------------------------------------  
      INTEGER                              , INTENT(in   ) :: i1, i2, j1, j2, k1, k2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2), INTENT(inout) :: ptab
      LOGICAL                              , INTENT(in   ) :: before
      !
      INTEGER :: ji, jj, jk
      !!----------------------------------------------------------------------  
      !    
      IF( before ) THEN
         ptab(i1:i2,j1:j2,k1:k2) = tmask(i1:i2,j1:j2,k1:k2) * e3t_0(i1:i2,j1:j2,k1:k2)
      ELSE
         !
         DO jk = k1, k2
            DO jj = j1, j2
               DO ji = i1, i2
                  IF( ABS( ptab(ji,jj,jk) - tmask(ji,jj,jk) * e3t_0(ji,jj,jk) ) > 1.D-2) THEN
                     WRITE(numout,*) ' Agrif error for e3t_0: parent , child, i, j, k ',  & 
                     &                 ptab(ji,jj,jk), tmask(ji,jj,jk) * e3t_0(ji,jj,jk), &
                     &                 ji+nimpp-1, jj+njmpp-1, jk
                     kindic_agr = kindic_agr + 1
                  ENDIF
               END DO
            END DO
         END DO
         !
      ENDIF
      ! 
   END SUBROUTINE interpe3t


   SUBROUTINE interpavm( ptab, i1, i2, j1, j2, k1, k2, m1, m2, before )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interavm  ***
      !!----------------------------------------------------------------------  
      INTEGER                                    , INTENT(in   ) ::   i1, i2, j1, j2, k1, k2, m1, m2
      REAL(wp),DIMENSION(i1:i2,j1:j2,k1:k2,m1:m2), INTENT(inout) ::   ptab
      LOGICAL                                    , INTENT(in   ) ::   before
      !
      INTEGER  :: ji, jj, jk
      INTEGER  :: N_in, N_out
      REAL(wp), DIMENSION(k1:k2) :: tabin, z_in
      REAL(wp), DIMENSION(1:jpk) :: z_out
      !!----------------------------------------------------------------------  
      !      
      IF (before) THEN         
         DO jk=k1,k2
            DO jj=j1,j2
              DO ji=i1,i2
                    ptab(ji,jj,jk,1) = avm_k(ji,jj,jk)
              END DO
           END DO
        END DO

# if defined key_vertical
        ! Interpolate thicknesses
        ! Warning: these are masked, hence extrapolated prior interpolation.
        DO jk=k1,k2
           DO jj=j1,j2
              DO ji=i1,i2
                  ptab(ji,jj,jk,2) = tmask(ji,jj,jk) * e3t(ji,jj,jk,Kmm_a)
              END DO
           END DO
        END DO

        ! Extrapolate thicknesses in partial bottom cells:
        ! Set them to Agrif_SpecialValue (0.). Correct bottom thicknesses are retrieved later on
        IF (ln_zps) THEN
           DO jj=j1,j2
              DO ji=i1,i2
                  jk = mbkt(ji,jj)
                  ptab(ji,jj,jk,2) = 0._wp
              END DO
           END DO           
        END IF
     
        ! Save ssh at last level:
        IF (.NOT.ln_linssh) THEN
           ptab(i1:i2,j1:j2,k2,2) = ssh(i1:i2,j1:j2,Kmm_a)*tmask(i1:i2,j1:j2,1) 
        ELSE
           ptab(i1:i2,j1:j2,k2,2) = 0._wp
        END IF      
# endif
      ELSE 
#ifdef key_vertical         
         IF (ln_linssh) ptab(i1:i2,j1:j2,k2,2) = 0._wp 
         avm_k(i1:i2,j1:j2,k1:k2) = 0._wp
            
         DO jj = j1, j2
            DO ji =i1, i2
               N_in = mbkt_parent(ji,jj)
               IF ( tmask(ji,jj,1) == 0._wp) N_in = 0
               z_in(N_in+1) = ht0_parent(ji,jj) + ptab(ji,jj,k2,2)
               DO jk = N_in, 1, -1  ! Parent vertical grid               
                     z_in(jk) = z_in(jk+1) - ptab(ji,jj,jk,2)
                    tabin(jk) = ptab(ji,jj,jk,1)
               END DO
               N_out = mbkt(ji,jj) 
               DO jk = 1, N_out        ! Child vertical grid
                  z_out(jk) = gdepw(ji,jj,jk,Kmm_a)
               ENDDO
               IF (N_in*N_out > 0) THEN
                  CALL remap_linear(tabin(1:N_in),z_in(1:N_in),avm_k(ji,jj,1:N_out),z_out(1:N_out),N_in,N_out,1)
               ENDIF
            ENDDO
         ENDDO
#else
         avm_k(i1:i2,j1:j2,k1:k2) = ptab (i1:i2,j1:j2,k1:k2,1)
#endif
      ENDIF
      !
   END SUBROUTINE interpavm

# if defined key_vertical
   SUBROUTINE interpmbkt( ptab, i1, i2, j1, j2, before )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpsshn  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      !
      !!----------------------------------------------------------------------  
      !
      IF( before) THEN
         ptab(i1:i2,j1:j2) = REAL(mbkt(i1:i2,j1:j2),wp)
      ELSE
         mbkt_parent(i1:i2,j1:j2) = NINT(ptab(i1:i2,j1:j2))
      ENDIF
      !
   END SUBROUTINE interpmbkt

   SUBROUTINE interpht0( ptab, i1, i2, j1, j2, before )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE interpsshn  ***
      !!----------------------------------------------------------------------  
      INTEGER                         , INTENT(in   ) ::   i1, i2, j1, j2
      REAL(wp), DIMENSION(i1:i2,j1:j2), INTENT(inout) ::   ptab
      LOGICAL                         , INTENT(in   ) ::   before
      !
      !!----------------------------------------------------------------------  
      !
      IF( before) THEN
         ptab(i1:i2,j1:j2) = ht_0(i1:i2,j1:j2)
      ELSE
         ht0_parent(i1:i2,j1:j2) = ptab(i1:i2,j1:j2)
      ENDIF
      !
   END SUBROUTINE interpht0
#endif

#else
   !!----------------------------------------------------------------------
   !!   Empty module                                          no AGRIF zoom
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE Agrif_OCE_Interp_empty
      WRITE(*,*)  'agrif_oce_interp : You should not have seen this print! error?'
   END SUBROUTINE Agrif_OCE_Interp_empty
#endif

   !!======================================================================
END MODULE agrif_oce_interp
