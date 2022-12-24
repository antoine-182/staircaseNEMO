#!/bin/bash

# NEMO directory where to fetch compiled STATION_ASF nemo.exe + setup:
NEMO_DIR="${HOME}/NEMO/NEMOvdev_r11085_ASINTER-05_Brodeau_Advanced_Bulk"

# Directory where to run the simulation:
WORK_DIR="${HOME}/tmp/STATION_ASF"


# FORC_DIR => Directory containing sea-surface + atmospheric forcings
#             (get it there https://drive.google.com/file/d/1MxNvjhRHmMrL54y6RX7WIaM9-LGl--ZP/):
if [ `hostname` = "merlat"        ]; then
    FORC_DIR="/MEDIA/data/STATION_ASF/input_data_STATION_ASF_2016-2018"
elif [ `hostname` = "luitel"        ]; then
    FORC_DIR="/data/gcm_setup/STATION_ASF/input_data_STATION_ASF_2016-2018"
elif [ `hostname` = "ige-meom-cal1" ]; then
    FORC_DIR="/mnt/meom/workdir/brodeau/STATION_ASF/input_data_STATION_ASF_2016-2018"
elif [ `hostname` = "salvelinus" ]; then
    FORC_DIR="/opt/data/STATION_ASF/input_data_STATION_ASF_2016-2018"
else
    echo "Boo!"; exit
fi
#======================
mkdir -p ${WORK_DIR}

NEMO_EXE="${NEMO_DIR}/tests/STATION_ASF/BLD/bin/nemo.exe"
if [ ! -f ${NEMO_EXE} ]; then echo " Mhhh, no compiled nemo.exe found into ${NEMO_DIR}/tests/STATION_ASF/BLD/bin !"; exit; fi

NEMO_EXPREF="${NEMO_DIR}/tests/STATION_ASF/EXPREF"
if [ ! -d ${NEMO_EXPREF} ]; then echo " Mhhh, no EXPREF directory ${NEMO_EXPREF} !"; exit; fi

rsync -avP ${NEMO_EXE}          ${WORK_DIR}/

for ff in "context_nemo.xml" "domain_def_nemo.xml" "field_def_nemo-oce.xml" "file_def_nemo-oce.xml" "grid_def_nemo.xml" "iodef.xml" "namelist_ref"; do
    if [ ! -f ${NEMO_EXPREF}/${ff} ]; then echo " Mhhh, ${ff} not found into ${NEMO_EXPREF} !"; exit; fi
    rsync -avPL ${NEMO_EXPREF}/${ff} ${WORK_DIR}/
done

# Copy forcing to work directory:
rsync -avP ${FORC_DIR}/Station_PAPA_50N-145W*.nc ${WORK_DIR}/

for CASE in "ECMWF-noskin" "COARE3p6-noskin" "ECMWF" "COARE3p6" "NCAR"; do

    echo ; echo
    echo "============================="
    echo " Going for ${CASE} experiment"
    echo "============================="
    echo

    scase=`echo "${CASE}" | tr '[:upper:]' '[:lower:]'`

    rm -f ${WORK_DIR}/namelist_cfg
    rsync -avPL ${NEMO_EXPREF}/namelist_${scase}_cfg ${WORK_DIR}/namelist_cfg

    cd ${WORK_DIR}/
    echo
    echo "Launching NEMO !"
    ./nemo.exe 1> out_nemo.out 2>err_nemo.err
    echo "Done!"
    echo

    # Moving output files:
    mkdir -p output
    mv -f STATION_ASF-${CASE}_*_grid*.nc output/
    
    # Saving logs:
    mkdir -p ${CASE}_log
    mv -f *.out *.err ocean.output output.namelist.dyn ${CASE}_log/

done
