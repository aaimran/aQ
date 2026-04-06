# Default values
nprocs=40
fname="test"



module purge
module load compiler-rt/2024.0.0 ifort/2024.0.0 mpi/2021.13 R

mpirun -np $nprocs ../../build/./waveqlab3d ${fname}.in | tee ${fname}.out


