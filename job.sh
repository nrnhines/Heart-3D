#!/bin/bash
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=68
#SBATCH --account=HBP_CDP2_it_0
#SBATCH --time=2:00:00
#SBATCH --mem=86000
#SBATCH --partition=knl_usr_prod
##SBATCH --qos=knl_qos_bprod
##SBATCH --qos=qos_prio
#SBATCH --mail-type=END
#SBATCH --mail-user=michael.hines@yale.edu
#SBATCH --constraint=cache

module load env-knl
module load intel/pe-xe-2017--binary
module load intelmpi/2017--binary
module load python/2.7.12

cd $SLURM_SUBMIT_DIR
mpirun nrniv -mpi -c 'stoptime=1000' parinit.hoc
