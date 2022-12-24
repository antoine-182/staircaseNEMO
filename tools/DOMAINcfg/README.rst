================================
= HOW TO COMPILE
================================
The DOMAINcfg tool can be compiled using the maketools script in the NEMOGCM/TOOLS directory as follows:
:::::::::::::::::::::::::::::::::
./maketools -m $ARCH -n DOMAINcfg
:::::::::::::::::::::::::::::::::
where $ARCH indicates the arch file to be used from the directory NEMOGCM/ARCH. 

DOMAINcfg compiled will create "make_domain_cfg.exe" executable script (the main routine of this tool is make_domain_cfg.f90)

================================
= HOW TO RUN
================================
just run :
:::::::::::::::::::::::::::::::::::::::::::::::::
mpirun -np 1 ./make_domain_cfg.exe
:::::::::::::::::::::::::::::::::::::::::::::::::
NOTA: it can be run in multiproc mode, but in output there will be domain_cfg_00xx.nc files

================================
= HOW TO USE
================================
1) copy in DOMAINcfg  directory namelist_cfg all settings (that you had in 3.6_stable) of the configuration for which you want prepare domain_cfg.nc file
IMPORTANT : keep the namelist_ref committed inchanged. !!!

NEW OPTION:
* ln_e3_dep
till nemo_v3.6_stable e3 were done like an analytical derivative of depth function
now  e3=dk[depth] in discret sens

If you want to create same e3[tuvw] like 3.6 you've to use "ln_e3_dep=.false."

 ln_e3_dep   = .true.    ! =T : e3=dk[depth] in discret sens.
   !                       !      ===>>> will become the only possibility in v4.0
   !                       ! =F : e3 analytical derivative of depth function
   !                       !      only there for backward compatibility test with v3.6
   !

* ln_domclo
The definition of closed sea are now integrated into ./make_domain_cfg.exe. The need on an external python script and hard coded indices has been removed.
ln_domclo=T will generate the masks needed to mask or correct 
the fwf unbalance over closed sea. If non defined closed sea are detected, a mask containing this will also be generated.
The detection of the lake is done using a lat/lon seed and a flood filling algorithm (the default namelist should work for any resolution).
The definition of the target area is done by defining a lat/lon position center, then all wet points at a defined distance of it will be selected.
Options are available to select only coastal points if needed.
For each lake, name, lon_src, lat_src, lon_trg, lat_trg, river mouth area, correction scheme, radius trg, id trg need to be defined.
* lon/lat_src/trg are the seed location of the closed sea (src) and its river outflow (trg).
* river mouth area is to defined is you spread the closed sea emp correction locally (circle with raduis 'radius trg' in m) 
  or along the coast point into a circle with raduis 'radius trg' (in m) 
  or globally (open sea mask)
* correction scheme defined how you want to spread the closed sea correction 
  * glo : as emp globally
  * rnf : as a runoff for the net precip and on emp for net evaporation case locally or along the coast
  * emp : as emp locally or along the coast
* id trg is a integer used to defined group of closed sea and their river outflow (for example all the great lake and the St Laurent outflow).

2) copy in DOMAINcfg directory same input files (of related configuration) required in v3.6_stable.

DOMAINcfg package is EXACTLY what does exist in NEMO version 3.6 to define a model domain (both domain related namelist and initialization).
DOMAINcfg tool creates a netcdf file "domain_cfg.nc" containing all the ocean domain informations required to define an ocean configuration,
these files are :
          
    domain size
    domain characteristics (periodic)
    horizontal mesh
    Coriolis parameter
    depth and vertical scale factors

FOR EXAMPLE 
- for AMM12 : 
            coordinates.nc
            bathy_meter.nc
            bathy_level.nc
            amm12_rivers.nc
            coordinates.bdy.nc
            amm12_restart_oce.nc
         directories:   bdydta/
                        fluxes/
- for ORCA2 : 
            coordinates.nc
            bathy_meter.nc
            bathy_level.nc
            domain_def.xml
            field_def.xml
            iodef.xml
