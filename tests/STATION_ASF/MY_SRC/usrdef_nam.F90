MODULE usrdef_nam
   !!======================================================================
   !!                     ***  MODULE usrdef_nam   ***
   !!
   !!                     ===  STATION_ASF configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni, G. Madec)  Original code
   !! History :  4.x  ! 2019-10  (L. Brodeau) for STATION_ASF (C1D meets SAS)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp, njmpp             ! ocean space and time domain
   USE dom_oce  , ONLY: ln_zco, ln_zps, ln_sco   ! flag of type of coordinate
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called in nemogcm.F90 module

   !                              !!* namusr_def namelist *!!
   REAL(wp), PUBLIC::   rn_dept1   ! Depth (m) at which the SST is taken/measured == depth of first T point!

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 11536 2019-09-11 13:54:18Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_nam( cd_cfg, kk_cfg, kpi, kpj, kpk, kperio )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read user defined namelist and define the domain size
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here STATION_ASF configuration
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(out) ::   cd_cfg          ! configuration name
      INTEGER         , INTENT(out) ::   kk_cfg          ! configuration resolution
      INTEGER         , INTENT(out) ::   kpi, kpj, kpk   ! global domain sizes 
      INTEGER         , INTENT(out) ::   kperio          ! lateral global domain b.c. 
      !
      INTEGER ::   ios   ! Local integer
      !!
      NAMELIST/namusr_def/ rn_dept1
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )
      !
      cd_cfg = 'STATION_ASF'               ! name & resolution (not used)
      kk_cfg = 0

      ! Global Domain size: STATION_ASF domain is 3 x 3 grid-points x 75 or vertical levels
      kpi = 3
      kpj = 3
      kpk = 1
      !
      !                             ! Set the lateral boundary condition of the global domain
      kperio =  7                   ! C1D configuration : 3x3 basin with cyclic Est-West and Norht-South condition
      !
      !                             ! control print
      IF(lwp) THEN
         WRITE(numout,*) '   '
         WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namusr_def : STATION_ASF case'
         WRITE(numout,*) '         Detpth of first T-point (bulk SST): rn_dept1 = ', rn_dept1
         WRITE(numout,*) '         jpiglo, jpjglo  = ', kpi, kpj
         WRITE(numout,*) '      number of model levels                              kpk = ', kpk
         WRITE(numout,*) '   '
         WRITE(numout,*) '   Lateral b.c. of the domain set to       jperio = ', kperio
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
