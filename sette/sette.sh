#!/bin/sh
# initialise user dependent variable
SETTE_DIR=$(cd $(dirname "$0"); pwd)
MAIN_DIR=$(dirname $SETTE_DIR)

# Parse command-line arguments
#
#  -t "test configurations": select active test configurations
#                            (overrides selection made using
#                            environment variable SETTE_TEST_CONFIGS)

if [ $# -gt 0 ]; then
  while getopts :ht: option; do 
     case $option in
        t) export SETTE_TEST_CONFIGS=$OPTARG
           echo ""
           echo "Configuration $SETTE_TEST_CONFIGS will be tested if they are available"
           echo "" ;;
        h | *) echo 'sette.sh allow no arguments (all configuration will be tested or -t "CFG1_to_test CFG2_to_test ..."'; exit 42 ;;
     esac
  done
  shift $((OPTIND - 1))
fi

. ./param.cfg

# run sette on reference configuration
./sette_reference-configurations.sh
if [[ $? != 0 ]]; then
   echo ""
   echo "--------------------------------------------------------------"
   echo "./sette_cfg-ref.sh didn't finish properly, need investigations"
   echo "--------------------------------------------------------------"
   echo ""
   exit 42
fi

# run sette on test cases
./sette_test-cases.sh
if [[ $? != 0 ]]; then
   echo ""
   echo "-----------------------------------------------------------------"
   echo "./sette_test-cases.sh didn't finish properly, need investigations"
   echo "-----------------------------------------------------------------"
   echo ""
   exit 42
fi

# run sette report
echo ""
echo "-------------------------------------------------------------"
echo "./sette_rpt.sh (script will wait all nemo_sette run are done)"
echo "-------------------------------------------------------------"
echo ""
NRUN=999
NIT=0
while [[ $NRUN -ne 0 && $nit -le 1080 ]]; do
   nit=$((nit+1))
   NRUN=$( ${BATCH_STAT} | grep nemo_sette | wc -l ) 
   if [[ $NRUN -ne 0 ]]; then 
      printf "%-3d %s\r" $NRUN 'nemo_sette run still in queue or running ...';
   else
      printf "%-50s\n" " "
      ./sette_rpt.sh
      exit
   fi
   sleep 10
done
printf "\n"
echo ""
echo "Something wrong happened, it tooks more than 3 hours to run all the sette tests"
echo ""
