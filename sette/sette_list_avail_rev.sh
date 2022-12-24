#!/bin/bash -f
# set -vx
#
lst_rev () {
    # get the list of revision available for a configuration
    if [ ! -d $1 ] ; then 
       CFGLST=-9999
    else
       CFGLST=`ls $1 | sort -u -r ` 
    fi
    # config name
    CONFIG=$2
    # list of all revision available
    ALLLST=${@:3}
    # number of revision total and for CONFIG
    nrevall=`echo $ALLLST | wc -w`
    nrevcfg=`echo $CFGLST | wc -w`
    # display
    echo ""
    printf "%-27s : " $CONFIG
    irev=1
    irevcfg=1
    while [[ $irev -le $nrevall ]] ; do
       rev=`echo $ALLLST | cut -d\  -f ${irev}`
       cfgrev=`echo $CFGLST | cut -d\  -f ${irevcfg}`
       if [ $cfgrev == $rev ] ; then
          printf "%-6s  " $rev
          irevcfg=$((irevcfg+1))
       else
          printf "%-5s  " "***** " 
       fi
       irev=$((irev+1))
    done
}

  SETTE_DIR=$(cd $(dirname "$0"); pwd)
  MAIN_DIR=$(dirname $SETTE_DIR)
  . ./param.cfg

  mach=${COMPILER}
  NEMO_VALID=${NEMO_VALIDATION_DIR}
 
 # list of all revision available
 DIRLST=`find ${NEMO_VALID} -maxdepth 3 -mindepth 3 -type d | sed -e 's/.*\/W.*\///' | sort -u -r`

 # display header
 echo ""
 echo " Compiler used is : $COMPILER"
 echo ""
 printf " List of all avail. rev. is : "
 for dir in `echo $DIRLST`; do printf "%5s  " $dir ; done
 printf "\n"

 # start checking configuration revision
 echo ""
 echo "   !---- check revision available for each configuration ----!   "
 for CONFIG in WGYRE_PISCES_ST WORCA2_ICE_PISCES_ST WORCA2_OFF_PISCES_ST WAMM12_ST WORCA2_SAS_ICE_ST WAGRIF_DEMO_ST WSPITZ12_ST WISOMIP_ST WOVERFLOW_ST WLOCK_EXCHANGE_ST WVORTEX_ST WICE_AGRIF_ST 
 do
    DIR=${NEMO_VALIDATION_DIR}/${CONFIG}/${COMPILER}
    lst_rev $DIR $CONFIG $DIRLST
 done
 printf "\n"
 printf "\n"
#
exit
