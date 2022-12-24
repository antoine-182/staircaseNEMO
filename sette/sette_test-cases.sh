#!/bin/bash
############################################################
# Author : Simona Flavoni for NEMO
# Contact: sflod@locean-ipsl.upmc.fr
#
# sette_test-cases.sh   : principal script of SET TEsts for NEMO (SETTE)
#                       : this script : compiles, run and tests TEST_CASES
#
#                       : TO DO: test if nitend is equal to end of run.stat
# ----------------------------------------------------------------------
# NEMO/SETTE , NEMO Consortium (2018)
# Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
# ----------------------------------------------------------------------
#
#############################################################
#set -vx
set -o posix
#set -u
#set -e
# ===========
# DESCRIPTION
# ===========
#
# Variables to be checked by user:
#
# COMPILER          : name of compiler as defined in NEMOGCM/ARCH directory 
# BATCH_COMMAND_PAR :  name of the command for submitting parallel batch jobs
# BATCH_COMMAND_SEQ :  name of the command for submitting sequential batch jobs  
# INTERACT_FLAG     : flag to run in interactive mode "yes"
#                           to run in batch mode "no"
# MPIRUN_FLAG       : flag to run in parallel (MPI) "yes"
#                           to run in sequential mode (NB_PROC = 1) "no"
# USING_XIOS        : flag to control the activation of key_iomput
#                      "yes" to compile using key_iomput and link to the external XIOS library
#                      "no"  to compile without key_iomput and link to the old IOIPSL library
# USING_MPMD        : flag to control the use of stand-alone IO servers
#                     requires USING_XIOS="yes"
#                      "yes" to run in MPMD (detached) mode with stand-alone IO servers
#                      "no"  to run in SPMD (attached) mode without separate IO servers 
# NUM_XIOSERVERS    : number of stand-alone IO servers to employ
#                     set to zero if USING_MPMD="no"
#
# Principal script is sette_test-cases.sh, that calls 
#
#  makenemo  : to create successive exectuables in ${CONFIG_NAME}/BLD/bin/nemo.exe 
#              and links to nemo in ${CONFIG_NAME}/EXP00)
#
#  param.cfg : sets and loads following directories:
#
#   FORCING_DIR         : is the directory for forcing files (tarfile)
#   INPUT_DIR           : is the directory for input files storing 
#   TMPDIR              : is the temporary directory (if needed)
#   NEMO_VALIDATION_DIR : is the validation directory
#
#   (NOTE: this file is the same for all configrations to be tested with sette_test-cases.sh)
#
#   all_functions.sh : loads functions used by sette (note: new functions can be added here)
#   set_namelist     : function declared in all_functions that sets namelist parameters 
#   post_test_tidyup : creates validation storage directory and copies required output files 
#                      (run.stat and ocean.output) in it after execution of test.
#
#  VALIDATION tree is:
#
#   NEMO_VALIDATION_DIR/WCONFIG_NAME/WCOMPILER_NAME/TEST_NAME/REVISION_NUMBER(or DATE)
#
#  prepare_exe_dir.sh : defines and creates directory where the test is executed
#                       execution directory takes name of TEST_NAME defined for every test 
#                       in sette_test-cases.sh. (each test in executed in its own directory)
#
#  set_valid_dir       : rename ocean.output/run.stat and tracer.stat to avoid checking them in the report 
#
#  clean_valid_dir    : rename ocean.output/run.stat and tracer.stat to avoid checking them in the report 
#                       ( not doing it could lead to false positive )
#
#  prepare_job.sh     : to generate the script run_job.sh
#
#  fcm_job.sh         : run in batch (INTERACT_FLAG="no") or interactive (INTERACT_FLAG="yes")
#                        see sette_test-cases.sh and BATCH_TEMPLATE directory
#
#  NOTE: jobs requiring initial or forcing data need to have an input_CONFIG.cfg in which 
#        can be found paths to the input tar file)
#  NOTE: if job is not launched for any reason you have the executable ready in ${EXE_DIR} 
#        directory
#  NOTE: the changed namelists are left in ${EXE_DIR} directory whereas original namelists 
#        remain in ${NEW_CONF}/EXP00
# 
#  NOTE: a log file, output.sette, is created in ${SETTE_DIR} with the echoes of 
#        executed commands
#
#  NOTE: if sette_test-cases.sh is stopped in output.sette there is written the last command 
#        executed by sette_test-cases.sh
#
# example use: ./sette_test-cases.sh 
#########################################################################################
#
# LOAD param value
SETTE_DIR=$(cd $(dirname "$0"); pwd)
MAIN_DIR=$(dirname $SETTE_DIR)
. ./param.cfg

export BATCH_COMMAND_PAR=${BATCH_CMD}
export BATCH_COMMAND_SEQ=${BATCH_CMD}
export INTERACT_FLAG="no"
export MPIRUN_FLAG="yes"
#
export DEL_KEYS="key_iomput"
if [ ${USING_XIOS} == "yes" ] 
 then 
   export DEL_KEYS=""
fi
#
export ADD_KEYS=""
if [ ${ADD_NOSIGNEDZERO} == "yes" ]
 then
   export ADD_KEYS="key_nosignedzero"
fi
#
# Settings which control the use of stand alone servers (only relevant if using xios)
#
export NUM_XIOSERVERS=4
export JOB_PREFIX=${JOB_PREFIX_MPMD}
#
if [ ${USING_MPMD} == "no" ] 
 then
   export NUM_XIOSERVERS=0
   export JOB_PREFIX=${JOB_PREFIX_NOMPMD}
fi
#
#
if [ ${USING_MPMD} == "yes" ] && [ ${USING_XIOS} == "no" ]
 then
   echo "Incompatible choices. MPMD mode requires the XIOS server"
   exit
fi

# Directory to run the tests
CONFIG_DIR0=${MAIN_DIR}/cfgs
TOOLS_DIR=${MAIN_DIR}/tools

CMP_NAM=${1:-$COMPILER}
# Copy job_batch_COMPILER file for specific compiler into job_batch_template
cd ${SETTE_DIR}
cp BATCH_TEMPLATE/${JOB_PREFIX}-${COMPILER} job_batch_template || exit
# Description of configuration tested:
# OVERFLOW       : TEST s-coordinates : (tracers) Advection schemes: FCT2, FCT4, ubs 
#                                     & (dynamics) advection schemes: flux form (ubs, centered), vector form (een)
#                       zps-coordinates : (tracers) Advection schemes: FCT2, FCT4, ubs
#                                     & (dynamics) advection schemes: flux form (ubs, centered), vector form (een, and een + Hollingsworth correction)
# LOCK_EXCHANGE  : 
# VORTEX         : 
# ICE_AGRIF      : 
# ISOMIP         : 
# WAD


for config in ${TEST_CONFIGS}
do

# ---------
#  OVERFLOW
# ---------
if [ ${config} == "OVERFLOW" ] ;  then
    ## Restartability tests for OVERFLOW
    export TEST_NAME="LONG"
    cd ${MAIN_DIR}
    . ./makenemo -m ${CMP_NAM} -n OVERFLOW_ST -a OVERFLOW -j 8 add_key "${ADD_KEYS}" del_key "${DEL_KEYS}"
    cd ${SETTE_DIR}
    . ./param.cfg
    . ./all_functions.sh
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=1
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
    cd ${EXE_DIR}  
    set_namelist namelist_cfg cn_exp \"OVF_LONG\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_itend 120
    set_namelist namelist_cfg nn_stock 60
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_EMPTY.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}

    cd ${SETTE_DIR}
    export TEST_NAME="SHORT"
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"OVF_SHORT\"
    set_namelist namelist_cfg nn_it000 61
    set_namelist namelist_cfg nn_itend 120
    set_namelist namelist_cfg nn_stock 60
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    set_namelist namelist_cfg ln_rstart .true.
    set_namelist namelist_cfg nn_rstctl 2
    set_namelist namelist_cfg cn_ocerst_in \"OVF_LONG_00000060_restart\"
    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    ln -sf ../LONG/OVF_LONG_00000060_restart.nc .

    cd ${SETTE_DIR}
    . ./prepare_job.sh input_EMPTY.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}


    ## Test for all advection, vert. coordinates, vector form, flux form: test runability and complete all time steps
    ## Needed namelist-xxxx for every type of run tested
    cd ${CONFIG_DIR}/${NEW_CONF}/EXP00

    for file in $(echo `ls namelist_*_cfg `) ; do
        TEST_NAME=`echo $file | sed -e "s/namelist_//" | sed -e "s/_cfg//"`
        TEST_NAME="EXP-${TEST_NAME}"
        if [ ! -d ${CONFIG_DIR}/${NEW_CONF}/${TEST_NAME} ] ; then mkdir ${CONFIG_DIR}/${NEW_CONF}/${TEST_NAME} ; fi
        export TEST_NAME="${TEST_NAME}"
         ##
        cd ${SETTE_DIR}
        . ./param.cfg
        . ./all_functions.sh
        . ./prepare_exe_dir.sh
        set_valid_dir
        clean_valid_dir
        JOB_FILE=${EXE_DIR}/run_job.sh
        NPROC=1
        if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
        cd ${EXE_DIR}
        rm namelist_*_*_*_*
        cp -pL ${CONFIG_DIR}/${NEW_CONF}/EXP00/$file namelist_cfg
        if [ ${USING_MPMD} == "yes" ] ; then
           set_xio_using_server iodef.xml true
        else
           set_xio_using_server iodef.xml false
        fi
        cd ${SETTE_DIR}
        . ./prepare_job.sh input_EMPTY.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
        cd ${SETTE_DIR}
        . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}
       ##
     done
fi

# --------------
#  LOCK_EXCHANGE
# --------------
if [ ${config} == "LOCK_EXCHANGE" ] ;  then
    ## Restartability tests for LOCK_EXCHANGE
    export TEST_NAME="LONG"
    cd ${MAIN_DIR}
    . ./makenemo -m ${CMP_NAM} -n LOCK_EXCHANGE_ST -a LOCK_EXCHANGE -j 8 add_key "${ADD_KEYS}" del_key "${DEL_KEYS}"
    cd ${SETTE_DIR}
    . ./param.cfg
    . ./all_functions.sh
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=1
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"LOCK_LONG\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_stock 60
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    set_namelist namelist_cfg nn_itend 120
    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_EMPTY.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}

    cd ${SETTE_DIR}
    export TEST_NAME="SHORT"
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"LOCK_SHORT\"
    set_namelist namelist_cfg nn_it000 61
    set_namelist namelist_cfg nn_itend 120
    set_namelist namelist_cfg nn_stock 60
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    set_namelist namelist_cfg ln_rstart .true.
    set_namelist namelist_cfg nn_rstctl 2
    set_namelist namelist_cfg cn_ocerst_in \"LOCK_LONG_00000060_restart\"
    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    ln -sf ../LONG/LOCK_LONG_00000060_restart.nc .

    cd ${SETTE_DIR}
    . ./prepare_job.sh input_EMPTY.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}

    ## Test for all advection, vector form, flux form: test runability and complete all time steps
    ## Needed namelist-xxxx for every type of run tested
    cd ${CONFIG_DIR}/${NEW_CONF}/EXP00

    for file in $(echo `ls namelist_*_cfg `) ; do
        echo ''
        TEST_NAME=`echo $file | sed -e "s/namelist_//" | sed -e "s/_cfg//"`
        TEST_NAME="EXP-${TEST_NAME}"
        `mkdir ${CONFIG_DIR}/${NEW_CONF}/${TEST_NAME}`
        export TEST_NAME="${TEST_NAME}"
        ##  
        cd ${SETTE_DIR}
        . ./param.cfg
        . ./all_functions.sh
        . ./prepare_exe_dir.sh
        set_valid_dir
        clean_valid_dir
        JOB_FILE=${EXE_DIR}/run_job.sh
        NPROC=1
        if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
        cd ${EXE_DIR}
        rm namelist_*_*_*_*
        cp -pL ${CONFIG_DIR}/${NEW_CONF}/EXP00/$file namelist_cfg
        if [ ${USING_MPMD} == "yes" ] ; then
           set_xio_using_server iodef.xml true
        else
           set_xio_using_server iodef.xml false
        fi
        cd ${SETTE_DIR}
        . ./prepare_job.sh input_EMPTY.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
        cd ${SETTE_DIR}
        . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}
        ##
        echo ''
   done
fi

# ---------
# VORTEX
# ---------
if [ ${config} == "VORTEX" ] ;  then 
## Restartability tests for VORTEX
    export TEST_NAME="LONG"
    cd ${MAIN_DIR}
    . ./makenemo -m ${CMP_NAM} -n VORTEX_ST -a VORTEX -j 8  add_key "${ADD_KEYS}" del_key "${DEL_KEYS}"
    cd ${SETTE_DIR}
    . ./param.cfg
    . ./all_functions.sh
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=6
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi 
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"VORTEX_LONG\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_itend 240
    set_namelist namelist_cfg nn_stock 120
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.

    set_namelist 1_namelist_cfg cn_exp \"VORTEX_LONG\"
    set_namelist 1_namelist_cfg nn_it000 1
    set_namelist 1_namelist_cfg nn_itend 720
    set_namelist 1_namelist_cfg nn_stock 360
    set_namelist 1_namelist_cfg sn_cfctl%l_config .true.
    set_namelist 1_namelist_cfg sn_cfctl%l_runstat .true.
	
    if [ ${USING_MPMD} == "yes" ] ; then
        set_xio_using_server iodef.xml true
    else
        set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_EMPTY.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    
    cd ${SETTE_DIR}
    export TEST_NAME="SHORT"
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"VORTEX_SHORT\"
    set_namelist namelist_cfg nn_it000 121
    set_namelist namelist_cfg nn_itend 240
    set_namelist namelist_cfg nn_stock 120
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    set_namelist namelist_cfg ln_rstart .true.
    set_namelist namelist_cfg nn_rstctl 2
    set_namelist namelist_cfg cn_ocerst_in \"VORTEX_LONG_00000120_restart\"
    
    set_namelist 1_namelist_cfg cn_exp \"VORTEX_SHORT\"
    set_namelist 1_namelist_cfg nn_it000 361
    set_namelist 1_namelist_cfg nn_itend 720
    set_namelist 1_namelist_cfg sn_cfctl%l_config .true.
    set_namelist 1_namelist_cfg sn_cfctl%l_runstat .true.
    set_namelist 1_namelist_cfg nn_stock 360
    set_namelist 1_namelist_cfg ln_rstart .true.
    set_namelist 1_namelist_cfg nn_rstctl 2
    set_namelist 1_namelist_cfg cn_ocerst_in \"VORTEX_LONG_00000360_restart\"
      
    if [ ${USING_MPMD} == "yes" ] ; then
        set_xio_using_server iodef.xml true
    else
        set_xio_using_server iodef.xml false
    fi
    if [ $NPROC -eq 1 ] ;  then
        ln -sf ../LONG/VORTEX_LONG_00000120_restart.nc .
        ln -sf ../LONG/1_VORTEX_LONG_00000360_restart.nc .
    else
        for (( i=1; i<=$NPROC; i++)) ; do
            L_NPROC=$(( $i - 1 ))
            L_NPROC=`printf "%04d\n" ${L_NPROC}`
            ln -sf ../LONG/VORTEX_LONG_00000120_restart_${L_NPROC}.nc .
            ln -sf ../LONG/1_VORTEX_LONG_00000360_restart_${L_NPROC}.nc .
        done
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_EMPTY.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}

## Reproducibility tests for VORTEX
    export TEST_NAME="REPRO_2_3"
    cd ${MAIN_DIR}
    cd ${SETTE_DIR}
    . ./param.cfg
    . ./all_functions.sh
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=6
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"VORTEX_23\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_itend 240
    set_namelist namelist_cfg nn_stock 120
    set_namelist namelist_cfg jpni 2
    set_namelist namelist_cfg jpnj 3
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
   
    set_namelist 1_namelist_cfg cn_exp \"VORTEX_23\"
    set_namelist 1_namelist_cfg nn_it000 1
    set_namelist 1_namelist_cfg nn_itend 720
    set_namelist 1_namelist_cfg nn_stock 360
    set_namelist 1_namelist_cfg jpni 2
    set_namelist 1_namelist_cfg jpnj 3
    set_namelist 1_namelist_cfg sn_cfctl%l_config .true.
    set_namelist 1_namelist_cfg sn_cfctl%l_runstat .true.

    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_EMPTY.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}

    cd ${SETTE_DIR}
    export TEST_NAME="REPRO_3_2"
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=6
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
    cd ${EXE_DIR}

    set_namelist namelist_cfg cn_exp \"VORTEX_32\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_itend 240
    set_namelist namelist_cfg nn_stock 120
    set_namelist namelist_cfg jpni 3
    set_namelist namelist_cfg jpnj 2
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.

    set_namelist 1_namelist_cfg cn_exp \"VORTEX_32\"
    set_namelist 1_namelist_cfg nn_it000 1
    set_namelist 1_namelist_cfg nn_itend 720
    set_namelist 1_namelist_cfg nn_stock 360
    set_namelist 1_namelist_cfg jpni 3
    set_namelist 1_namelist_cfg jpnj 2
    set_namelist 1_namelist_cfg sn_cfctl%l_config .true.
    set_namelist 1_namelist_cfg sn_cfctl%l_runstat .true.

    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_EMPTY.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}

fi


# ---------
# ICE_AGRIF
# ---------
if [ ${config} == "ICE_AGRIF" ] ;  then 
## Restartability tests for ICE_AGRIF
    export TEST_NAME="LONG"
    cd ${MAIN_DIR}
    . ./makenemo -m ${CMP_NAM} -n ICE_AGRIF_ST -a ICE_AGRIF -j 8  add_key "${ADD_KEYS}" del_key "${DEL_KEYS}"
    cd ${SETTE_DIR}
    . ./param.cfg
    . ./all_functions.sh
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=6
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi 
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"ICE_AGRIF_LONG\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_itend 200
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    set_namelist namelist_cfg nn_stock 100
    
    set_namelist 1_namelist_cfg cn_exp \"ICE_AGRIF_LONG\"
    set_namelist 1_namelist_cfg nn_it000 1
    set_namelist 1_namelist_cfg nn_itend 600
    set_namelist 1_namelist_cfg sn_cfctl%l_config .true.
    set_namelist 1_namelist_cfg sn_cfctl%l_runstat .true.
    set_namelist 1_namelist_cfg nn_stock 300
    
    if [ ${USING_MPMD} == "yes" ] ; then
        set_xio_using_server iodef.xml true
    else
        set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_ICE_AGRIF.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    
    cd ${SETTE_DIR}
    export TEST_NAME="SHORT"
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"ICE_AGRIF_SHORT\"
    set_namelist namelist_cfg nn_it000 101
    set_namelist namelist_cfg nn_itend 200
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    set_namelist namelist_cfg nn_stock 100
    set_namelist namelist_cfg ln_rstart .true.
    set_namelist namelist_cfg nn_rstctl 2
    set_namelist namelist_cfg cn_ocerst_in \"ICE_AGRIF_LONG_00000100_restart\"
    set_namelist namelist_ice_cfg cn_icerst_in \"ICE_AGRIF_LONG_00000100_restart_ice\"
    
    set_namelist 1_namelist_cfg cn_exp \"ICE_AGRIF_SHORT\"
    set_namelist 1_namelist_cfg nn_it000 301
    set_namelist 1_namelist_cfg nn_itend 600
    set_namelist 1_namelist_cfg sn_cfctl%l_config .true.
    set_namelist 1_namelist_cfg sn_cfctl%l_runstat .true.
    set_namelist 1_namelist_cfg nn_stock 300
    set_namelist 1_namelist_cfg ln_rstart .true.
    set_namelist 1_namelist_cfg nn_rstctl 2
    set_namelist 1_namelist_cfg cn_ocerst_in \"ICE_AGRIF_LONG_00000300_restart\"
    set_namelist 1_namelist_ice_cfg cn_icerst_in \"ICE_AGRIF_LONG_00000300_restart_ice\"
    
    
    if [ ${USING_MPMD} == "yes" ] ; then
        set_xio_using_server iodef.xml true
    else
        set_xio_using_server iodef.xml false
    fi
    if [ $NPROC -eq 1 ] ;  then
        ln -sf ../LONG/ICE_AGRIF_LONG_00000100_restart.nc .
        ln -sf ../LONG/ICE_AGRIF_LONG_00000100_restart_ice.nc .
        ln -sf ../LONG/1_ICE_AGRIF_LONG_00000300_restart.nc .
        ln -sf ../LONG/1_ICE_AGRIF_LONG_00000300_restart_ice.nc .
    else
        for (( i=1; i<=$NPROC; i++)) ; do
            L_NPROC=$(( $i - 1 ))
            L_NPROC=`printf "%04d\n" ${L_NPROC}`
            ln -sf ../LONG/ICE_AGRIF_LONG_00000100_restart_${L_NPROC}.nc .
            ln -sf ../LONG/ICE_AGRIF_LONG_00000100_restart_ice_${L_NPROC}.nc .
            ln -sf ../LONG/1_ICE_AGRIF_LONG_00000300_restart_${L_NPROC}.nc .
            ln -sf ../LONG/1_ICE_AGRIF_LONG_00000300_restart_ice_${L_NPROC}.nc .
        done
    fi

    cd ${SETTE_DIR}
    . ./prepare_job.sh input_ICE_AGRIF.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}


## Reproducibility tests for ICE_AGRIF
    export TEST_NAME="REPRO_2_3"
    cd ${MAIN_DIR}
    cd ${SETTE_DIR}
    . ./param.cfg
    . ./all_functions.sh
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=6
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"ICE_AGRIF_23\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_itend 200
    set_namelist namelist_cfg nn_stock 100
    set_namelist namelist_cfg jpni 2
    set_namelist namelist_cfg jpnj 3
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    
    set_namelist 1_namelist_cfg cn_exp \"ICE_AGRIF_23\"
    set_namelist 1_namelist_cfg nn_it000 1
    set_namelist 1_namelist_cfg nn_itend 600
    set_namelist 1_namelist_cfg nn_stock 300
    set_namelist 1_namelist_cfg jpni 2
    set_namelist 1_namelist_cfg jpnj 3
    set_namelist 1_namelist_cfg sn_cfctl%l_config .true.
    set_namelist 1_namelist_cfg sn_cfctl%l_runstat .true.

    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_ICE_AGRIF.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}

    cd ${SETTE_DIR}
    export TEST_NAME="REPRO_3_2"
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=6
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
    cd ${EXE_DIR}

    set_namelist namelist_cfg cn_exp \"ICE_AGRIF_32\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_itend 200
    set_namelist namelist_cfg nn_stock 100
    set_namelist namelist_cfg jpni 3
    set_namelist namelist_cfg jpnj 2
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    
    set_namelist 1_namelist_cfg cn_exp \"ICE_AGRIF_32\"
    set_namelist 1_namelist_cfg nn_it000 1
    set_namelist 1_namelist_cfg nn_itend 600
    set_namelist 1_namelist_cfg nn_stock 300
    set_namelist 1_namelist_cfg jpni 3
    set_namelist 1_namelist_cfg jpnj 2
    set_namelist 1_namelist_cfg sn_cfctl%l_config .true.
    set_namelist 1_namelist_cfg sn_cfctl%l_runstat .true.
    
    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_ICE_AGRIF.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}

fi

# ------
# ISOMIP
# ------
if [ ${config} == "ISOMIP" ] ;  then
## Restartability tests
    export TEST_NAME="LONG"
    cd ${MAIN_DIR}
    . ./makenemo -m ${CMP_NAM} -n ISOMIP_ST -a ISOMIP -j 8 add_key "${ADD_KEYS}" del_key "${DEL_KEYS}"
    cd ${SETTE_DIR}
    . ./param.cfg
    . ./all_functions.sh
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=15
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"ISOMIP_LONG\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_itend 96
    set_namelist namelist_cfg nn_stock 48
    set_namelist namelist_cfg jpni 5
    set_namelist namelist_cfg jpnj 3
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_ISOMIP.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}

    cd ${SETTE_DIR}
    export TEST_NAME="SHORT"
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"ISOMIP_SHORT\"
    set_namelist namelist_cfg nn_it000 49
    set_namelist namelist_cfg nn_itend 96
    set_namelist namelist_cfg nn_stock 48
    set_namelist namelist_cfg ln_rstart .true.
    set_namelist namelist_cfg nn_rstctl 2
    set_namelist namelist_cfg jpni 5
    set_namelist namelist_cfg jpnj 3
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    set_namelist namelist_cfg cn_ocerst_in \"ISOMIP_LONG_00000048_restart\"
    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    for (( i=1; i<=$NPROC; i++)) ; do
        L_NPROC=$(( $i - 1 ))
        L_NPROC=`printf "%04d\n" ${L_NPROC}`
        ln -sf ../LONG/ISOMIP_LONG_00000048_restart_${L_NPROC}.nc .
    done

    cd ${SETTE_DIR}
    . ./prepare_job.sh input_ISOMIP.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}

## Reproducibility tests
    export TEST_NAME="REPRO_7_3"
    cd ${MAIN_DIR}
    cd ${SETTE_DIR}
    . ./param.cfg
    . ./all_functions.sh
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=21
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"ISOMIP_73\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_itend 48
    set_namelist namelist_cfg jpni 7
    set_namelist namelist_cfg jpnj 3
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_ISOMIP.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}

    cd ${SETTE_DIR}
    export TEST_NAME="REPRO_8_4"
    . ./prepare_exe_dir.sh
    set_valid_dir
    clean_valid_dir
    JOB_FILE=${EXE_DIR}/run_job.sh
    NPROC=32
    if [ -f ${JOB_FILE} ] ; then \rm ${JOB_FILE} ; fi
    cd ${EXE_DIR}
    set_namelist namelist_cfg cn_exp \"ISOMIP_84\"
    set_namelist namelist_cfg nn_it000 1
    set_namelist namelist_cfg nn_itend 48
    set_namelist namelist_cfg jpni 8
    set_namelist namelist_cfg jpnj 4
    set_namelist namelist_cfg sn_cfctl%l_config .true.
    set_namelist namelist_cfg sn_cfctl%l_runstat .true.
    if [ ${USING_MPMD} == "yes" ] ; then
       set_xio_using_server iodef.xml true
    else
       set_xio_using_server iodef.xml false
    fi
    cd ${SETTE_DIR}
    . ./prepare_job.sh input_ISOMIP.cfg $NPROC ${TEST_NAME} ${MPIRUN_FLAG} ${JOB_FILE} ${NUM_XIOSERVERS} ${NEMO_VALID}
    cd ${SETTE_DIR}
    . ./fcm_job.sh $NPROC ${JOB_FILE} ${INTERACT_FLAG} ${MPIRUN_FLAG}

fi

#----
done
