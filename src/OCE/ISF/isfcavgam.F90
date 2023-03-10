MODULE isfcavgam
   !!======================================================================
   !!                       ***  MODULE  isfgammats  ***
   !! Ice shelf gamma module :  compute exchange coeficient at the ice/ocean interface
   !!======================================================================
   !! History :  4.1  !  (P. Mathiot) original
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   isfcav_gammats       : compute exchange coeficient gamma 
   !!----------------------------------------------------------------------
   USE isf_oce
   USE isfutils, ONLY: debug
   USE isftbl  , ONLY: isf_tbl

   USE oce     , ONLY: uu, vv, rn2         ! ocean dynamics and tracers
   USE phycst  , ONLY: grav, vkarmn        ! physical constant
   USE eosbn2  , ONLY: eos_rab             ! equation of state
   USE zdfdrg  , ONLY: rCd0_top, r_ke0_top ! vertical physics: top/bottom drag coef.
   USE iom     , ONLY: iom_put             !
   USE lib_mpp , ONLY: ctl_stop

   USE dom_oce        ! ocean space and time domain
   USE in_out_manager ! I/O manager
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   PUBLIC   isfcav_gammats
   
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcisf.F90 10536 2019-01-16 19:21:09Z mathiot $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   !
   !!-----------------------------------------------------------------------------------------------------
   !!                              PUBLIC SUBROUTINES 
   !!-----------------------------------------------------------------------------------------------------
   !
   SUBROUTINE isfcav_gammats( Kmm, pttbl, pstbl, pqoce, pqfwf, pgt, pgs )
      !!----------------------------------------------------------------------
      !! ** Purpose    : compute the coefficient echange for heat and fwf flux
      !!
      !! ** Method     : select the gamma formulation
      !!                 3 method available (cst, vel and vel_stab)
      !!---------------------------------------------------------------------
      !!-------------------------- OUT -------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) :: pgt  , pgs      ! gamma t and gamma s 
      !!-------------------------- IN  -------------------------------------
      INTEGER                                     :: Kmm             ! ocean time level index
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pqoce, pqfwf    ! isf heat and fwf
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pttbl, pstbl    ! top boundary layer tracer
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj)                :: zutbl, zvtbl    ! top boundary layer velocity
      !!---------------------------------------------------------------------
      !
      !==========================================
      ! 1.: compute velocity in the tbl if needed
      !==========================================
      !
      SELECT CASE ( cn_gammablk )
      CASE ( 'spe'  ) 
         ! gamma is constant (specified in namelist)
         ! nothing to do
      CASE ('vel', 'vel_stab')
         ! compute velocity in tbl
         CALL isf_tbl(Kmm, uu(:,:,:,Kmm) ,zutbl(:,:),'U', miku, rhisf_tbl_cav)
         CALL isf_tbl(Kmm, vv(:,:,:,Kmm) ,zvtbl(:,:),'V', mikv, rhisf_tbl_cav)
         !
         ! mask velocity in tbl with ice shelf mask
         zutbl(:,:) = zutbl(:,:) * mskisf_cav(:,:)
         zvtbl(:,:) = zvtbl(:,:) * mskisf_cav(:,:)
         !
         ! output
         CALL iom_put('utbl',zutbl(:,:))
         CALL iom_put('vtbl',zvtbl(:,:))
      CASE DEFAULT
         CALL ctl_stop('STOP','method to compute gamma (cn_gammablk) is unknown (should not see this)')
      END SELECT
      ! 
      !==========================================
      ! 2.: compute gamma
      !==========================================
      !
      SELECT CASE ( cn_gammablk )
      CASE ( 'spe'  ) ! gamma is constant (specified in namelist)
         pgt(:,:) = rn_gammat0
         pgs(:,:) = rn_gammas0
      CASE ( 'vel' ) ! gamma is proportional to u*
         CALL gammats_vel      (                   zutbl, zvtbl, rCd0_top, r_ke0_top,               pgt, pgs )
      CASE ( 'vel_stab' ) ! gamma depends of stability of boundary layer and u*
         CALL gammats_vel_stab (Kmm, pttbl, pstbl, zutbl, zvtbl, rCd0_top, r_ke0_top, pqoce, pqfwf, pgt, pgs )
      CASE DEFAULT
         CALL ctl_stop('STOP','method to compute gamma (cn_gammablk) is unknown (should not see this)')
      END SELECT
      !
      !==========================================
      ! 3.: output and debug
      !==========================================
      !
      CALL iom_put('isfgammat', pgt(:,:))
      CALL iom_put('isfgammas', pgs(:,:))
      !
      IF (ln_isfdebug) THEN
         CALL debug( 'isfcav_gam pgt:', pgt(:,:) )
         CALL debug( 'isfcav_gam pgs:', pgs(:,:) )
      END IF
      !
   END SUBROUTINE isfcav_gammats
   !
   !!-----------------------------------------------------------------------------------------------------
   !!                              PRIVATE SUBROUTINES 
   !!-----------------------------------------------------------------------------------------------------
   !
   SUBROUTINE gammats_vel( putbl, pvtbl, pCd, pke2, &   ! <<== in
      &                                  pgt, pgs   )   ! ==>> out gammats [m/s]
      !!----------------------------------------------------------------------
      !! ** Purpose    : compute the coefficient echange coefficient 
      !!
      !! ** Method     : gamma is velocity dependent ( gt= gt0 * Ustar )
      !!
      !! ** Reference  : Asay-Davis et al., Geosci. Model Dev., 9, 2471-2497, 2016
      !!---------------------------------------------------------------------
      !!-------------------------- OUT -------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) :: pgt, pgs     ! gammat and gammas [m/s]
      !!-------------------------- IN  -------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: putbl, pvtbl ! velocity in the losch top boundary layer
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pCd          ! drag coefficient
      REAL(wp),                     INTENT(in   ) :: pke2         ! background velocity
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj) :: zustar
      !!---------------------------------------------------------------------
      !
      ! compute ustar (AD15 eq. 27)
      zustar(:,:) = SQRT( pCd(:,:) * ( putbl(:,:) * putbl(:,:) + pvtbl(:,:) * pvtbl(:,:) + pke2 ) ) * mskisf_cav(:,:)
      !
      ! Compute gammats
      pgt(:,:) = zustar(:,:) * rn_gammat0
      pgs(:,:) = zustar(:,:) * rn_gammas0
      !
      ! output ustar
      CALL iom_put('isfustar',zustar(:,:))
      !
   END SUBROUTINE gammats_vel

   SUBROUTINE gammats_vel_stab( Kmm, pttbl, pstbl, putbl, pvtbl, pCd, pke2, pqoce, pqfwf, &  ! <<== in
      &                                                                     pgt  , pgs    )  ! ==>> out gammats [m/s]
      !!----------------------------------------------------------------------
      !! ** Purpose    : compute the coefficient echange coefficient 
      !!
      !! ** Method     : gamma is velocity dependent and stability dependent
      !!
      !! ** Reference  : Holland and Jenkins, 1999, JPO, p1787-1800
      !!---------------------------------------------------------------------
      !!-------------------------- OUT -------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) :: pgt, pgs     ! gammat and gammas
      !!-------------------------- IN  -------------------------------------
      INTEGER                                     :: Kmm            ! ocean time level index
      REAL(wp),                     INTENT(in   ) :: pke2           ! background velocity squared
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pqoce, pqfwf   ! surface heat flux and fwf flux
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pCd            ! drag coeficient
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: putbl, pvtbl   ! velocity in the losch top boundary layer
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) :: pttbl, pstbl   ! tracer   in the losch top boundary layer
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj                     ! loop index
      INTEGER  :: ikt                        ! local integer
      REAL(wp) :: zdku, zdkv                 ! U, V shear 
      REAL(wp) :: zPr, zSc, zRc              ! Prandtl, Scmidth and Richardson number 
      REAL(wp) :: zmob, zmols                ! Monin Obukov length, coriolis factor at T point
      REAL(wp) :: zbuofdep, zhnu             ! Bouyancy length scale, sublayer tickness
      REAL(wp) :: zhmax                      ! limitation of mol
      REAL(wp) :: zetastar                   ! stability parameter
      REAL(wp) :: zgmolet, zgmoles, zgturb   ! contribution of modelecular sublayer and turbulence 
      REAL(wp) :: zcoef                      ! temporary coef
      REAL(wp) :: zdep
      REAL(wp) :: zeps = 1.0e-20_wp    
      REAL(wp), PARAMETER :: zxsiN = 0.052_wp   ! dimensionless constant
      REAL(wp), PARAMETER :: znu   = 1.95e-6_wp ! kinamatic viscosity of sea water (m2.s-1)
      REAL(wp), DIMENSION(2) :: zts, zab
      REAL(wp), DIMENSION(jpi,jpj) :: zustar    ! friction velocity
      !!---------------------------------------------------------------------
      !
      ! compute ustar
      zustar(:,:) = SQRT( pCd * ( putbl(:,:) * putbl(:,:) + pvtbl(:,:) * pvtbl(:,:) + pke2 ) )
      !
      ! output ustar
      CALL iom_put('isfustar',zustar(:,:))
      !
      ! compute Pr and Sc number (eq ??)
      zPr =   13.8_wp
      zSc = 2432.0_wp
      !
      ! compute gamma mole (eq ??)
      zgmolet = 12.5_wp * zPr ** (2.0/3.0) - 6.0_wp
      zgmoles = 12.5_wp * zSc ** (2.0/3.0) - 6.0_wp
      !
      ! compute gamma
      DO ji = 2, jpi
         DO jj = 2, jpj
            ikt = mikt(ji,jj)

            IF( zustar(ji,jj) == 0._wp ) THEN           ! only for kt = 1 I think
               pgt = rn_gammat0
               pgs = rn_gammas0
            ELSE
               ! compute Rc number (as done in zdfric.F90)
!!gm better to do it like in the new zdfric.F90   i.e. avm weighted Ri computation
               zcoef = 0.5_wp / e3w(ji,jj,ikt+1,Kmm)
               !                                            ! shear of horizontal velocity
               zdku = zcoef * (  uu(ji-1,jj  ,ikt  ,Kmm) + uu(ji,jj,ikt  ,Kmm)  &
                  &             -uu(ji-1,jj  ,ikt+1,Kmm) - uu(ji,jj,ikt+1,Kmm)  )
               zdkv = zcoef * (  vv(ji  ,jj-1,ikt  ,Kmm) + vv(ji,jj,ikt  ,Kmm)  &
                  &             -vv(ji  ,jj-1,ikt+1,Kmm) - vv(ji,jj,ikt+1,Kmm)  )
               !                                            ! richardson number (minimum value set to zero)
               zRc = MAX(rn2(ji,jj,ikt+1), 0._wp) / MAX( zdku*zdku + zdkv*zdkv, zeps )

               ! compute bouyancy 
               zts(jp_tem) = pttbl(ji,jj)
               zts(jp_sal) = pstbl(ji,jj)
               zdep        = gdepw(ji,jj,ikt,Kmm)
               !
               CALL eos_rab( zts, zdep, zab, Kmm )
               !
               ! compute length scale (Eq ??)
               zbuofdep = grav * ( zab(jp_tem) * pqoce(ji,jj) - zab(jp_sal) * pqfwf(ji,jj) )
               !
               ! compute Monin Obukov Length
               ! Maximum boundary layer depth (Eq ??)
               zhmax = gdept(ji,jj,mbkt(ji,jj),Kmm) - gdepw(ji,jj,mikt(ji,jj),Kmm) -0.001_wp
               !
               ! Compute Monin obukhov length scale at the surface and Ekman depth: (Eq ??)
               zmob   = zustar(ji,jj) ** 3 / (vkarmn * (zbuofdep + zeps))
               zmols  = SIGN(1._wp, zmob) * MIN(ABS(zmob), zhmax) * tmask(ji,jj,ikt)
               !
               ! compute eta* (stability parameter) (Eq ??)
               zetastar = 1._wp / ( SQRT(1._wp + MAX(zxsiN * zustar(ji,jj) / ( ABS(ff_f(ji,jj)) * zmols * zRc ), 0._wp)))
               !
               ! compute the sublayer thickness (Eq ??)
               zhnu = 5 * znu / zustar(ji,jj)
               !
               ! compute gamma turb (Eq ??)
               zgturb = 1._wp / vkarmn * LOG(zustar(ji,jj) * zxsiN * zetastar * zetastar / ( ABS(ff_f(ji,jj)) * zhnu )) &
               &      + 1._wp / ( 2 * zxsiN * zetastar ) - 1._wp / vkarmn
               !
               ! compute gammats
               pgt(ji,jj) = zustar(ji,jj) / (zgturb + zgmolet)
               pgs(ji,jj) = zustar(ji,jj) / (zgturb + zgmoles)
            END IF
         END DO
      END DO

   END SUBROUTINE gammats_vel_stab

END MODULE isfcavgam
