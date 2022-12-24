#!/bin/bash
#SBATCH --job-name=a_COR              # nom du job
#SBATCH -A bho@cpu
###############SBATCH --qos=qos_cpu-dev             # test pas tres lourd 

#SBATCH --nodes=1                  # nombre de noeuds (1node = 2x20 proc)
#SBATCH --ntasks-per-node=1         # nombre de taches MPI par noeud
#SBATCH --hint=nomultithread       # 1 processus MPI par coeur physique (pas d'hyperthreading)
#SBATCH --exclusive

#SBATCH --time=00:45:00             # temps d execution maximum demande (HH:MM:SS)
#SBATCH --output=AM98_%j.out          # nom du fichier de sortie
#SBATCH --error=AM98_%j.out           # nom du fichier d'erreur (ici en commun avec la sortie)

# echo des commandes lancées
set -x

# on se place dans le repertoire de soumission
dirwork="EXP_ens_ns_divr" 

cd ../$dirwork/

# ref48 : 7x40 procs pour (0°), 10x40 (30°)
# ref4   : 40 procs 2h < 1/16° 45°
# nettoyage des modules charges en interactif et herites par defaut
module purge

# some cleaning
#rm -f AM98*_restart_*.nc
#rm -f mesh_mask*.nc
#rm -f output.init*.nc
#rm -f output.abort*.nc

 
# chargement des modules
intel_version=19.0.4
module load intel-compilers/${intel_version}
module load intel-mpi/${intel_version}
module load hdf5/1.10.5-mpi
module load netcdf/4.7.2-mpi
module load netcdf-fortran/4.5.2-mpi
module load cdo/1.9.7.1
module load nco/4.8.1
module load ncview/2.1.7-mpi
module load subversion/1.9.7
 
# execution
#./nemo 
# mpirun -n 8 ./nemo           # mal adapte a jeanzay
srun --qos=qos_cpu-dev ./nemo
#srun --qos=qos_cpu-t3 ./nemo
