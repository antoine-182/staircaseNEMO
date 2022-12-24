MODULE domcfg
   !!==============================================================================
   !!                       ***  MODULE domcfg   ***
   !! Ocean initialization : domain configuration initialization
   !!==============================================================================
   !! History :  1.0  ! 2003-09  (G. Madec)  Original code
   !!            3.2  ! 2009-07  (R. Benshila) Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dom_cfg        : initialize the domain configuration
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dom_cfg    ! called by opa.F90

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LODYC-IPSL  (2009)
   !! $Id: domcfg.F90 6140 2015-12-21 11:35:23Z timgraham $ 
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dom_cfg
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_cfg  ***
      !!                    
      !! ** Purpose :   set the domain configuration
      !!
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_cfg : set the ocean configuration'
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   ocean model configuration used :   cp_cfg = ', cp_cfg, ' jp_cfg = ', jp_cfg
         !
         WRITE(numout,*) '   global domain lateral boundaries'
         !
         IF( jperio == 0 )   WRITE(numout,*) '      jperio= 0, closed'
         IF( jperio == 1 )   WRITE(numout,*) '      jperio= 1, cyclic east-west'
         IF( jperio == 2 )   WRITE(numout,*) '      jperio= 2, equatorial symmetric'
         IF( jperio == 3 )   WRITE(numout,*) '      jperio= 3, north fold with T-point pivot'
         IF( jperio == 4 )   WRITE(numout,*) '      jperio= 4, cyclic east-west and north fold with T-point pivot'
         IF( jperio == 5 )   WRITE(numout,*) '      jperio= 5, north fold with F-point pivot'
         IF( jperio == 6 )   WRITE(numout,*) '      jperio= 6, cyclic east-west and north fold with F-point pivot'
      ENDIF
      !
      IF( jperio <  0 .OR. jperio > 6 )   CALL ctl_stop( 'jperio is out of range' )
      !
      CALL dom_glo                   ! global domain versus zoom and/or local domain
      !
   END SUBROUTINE dom_cfg

   SUBROUTINE dom_glo
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_glo  ***
      !!
      !! ** Purpose :   initialization of global domain <--> local domain indices
      !!
      !! ** Method  :   
      !!
      !! ** Action  : - mig , mjg : local  domain indices ==> global domain indices
      !!              - mi0 , mi1 : global domain indices ==> local  domain indices
      !!              - mj0,, mj1   (global point not in the local domain ==> mi0>mi1 and/or mj0>mj1)
      !!----------------------------------------------------------------------
      INTEGER ::   ji, jj   ! dummy loop argument
      !!----------------------------------------------------------------------
      !
      DO ji = 1, jpi                 ! local domain indices ==> global domain indices
        mig(ji) = ji + nimpp - 1
      END DO
      DO jj = 1, jpj
        mjg(jj) = jj + njmpp - 1
      END DO
      !                              ! global domain indices ==> local domain indices
      !                                   ! (return (m.0,m.1)=(1,0) if data domain gridpoint is to the west/south of the 
      !                                   ! local domain, or (m.0,m.1)=(jp.+1,jp.) to the east/north of local domain. 
      DO ji = 1, jpiglo
        mi0(ji) = MAX( 1 , MIN( ji - nimpp + 1, jpi+1 ) )
        mi1(ji) = MAX( 0 , MIN( ji - nimpp + 1, jpi   ) )
      END DO
      DO jj = 1, jpjglo
        mj0(jj) = MAX( 1 , MIN( jj - njmpp + 1, jpj+1 ) )
        mj1(jj) = MAX( 0 , MIN( jj - njmpp + 1, jpj   ) )
      END DO
      IF(lwp) THEN                   ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dom_glo : domain: global <<==>> local '
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   global domain:   jpiglo = ', jpiglo, ' jpjglo = ', jpjglo, ' jpkglo = ', jpkglo
         WRITE(numout,*) '   local  domain:   jpi    = ', jpi   , ' jpj    = ', jpj   , ' jpk    = ', jpk
         WRITE(numout,*)
         WRITE(numout,*) '   conversion from local to global domain indices (and vise versa) done'
         IF( nn_print >= 1 ) THEN
            WRITE(numout,*)
            WRITE(numout,*) '          conversion local  ==> global i-index domain (mig)'
            WRITE(numout,25)              (mig(ji),ji = 1,jpi)
            WRITE(numout,*)
            WRITE(numout,*) '          conversion global ==> local  i-index domain'
            WRITE(numout,*) '             starting index (mi0)'
            WRITE(numout,25)              (mi0(ji),ji = 1,jpiglo)
            WRITE(numout,*) '             ending index (mi1)'
            WRITE(numout,25)              (mi1(ji),ji = 1,jpiglo)
            WRITE(numout,*)
            WRITE(numout,*) '          conversion local  ==> global j-index domain (mjg)'
            WRITE(numout,25)              (mjg(jj),jj = 1,jpj)
            WRITE(numout,*)
            WRITE(numout,*) '          conversion global ==> local  j-index domain'
            WRITE(numout,*) '             starting index (mj0)'
            WRITE(numout,25)              (mj0(jj),jj = 1,jpjglo)
            WRITE(numout,*) '             ending index (mj1)'
            WRITE(numout,25)              (mj1(jj),jj = 1,jpjglo)
         ENDIF
      ENDIF
 25   FORMAT( 100(10x,19i4,/) )
      !
   END SUBROUTINE dom_glo
   !!======================================================================
END MODULE domcfg
