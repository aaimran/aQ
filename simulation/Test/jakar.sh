# Default values
nprocs=40


module purge
module load compiler-rt/2024.0.0 ifort/2024.0.0 mpi/2021.13 R

fname="test-3t_100m_m8"
mpirun -np $nprocs ../../build/./waveqlab3d ${fname}.in | tee ${fname}.out


fname="test-3t_100m_m8cg"
mpirun -np $nprocs ../../build/./waveqlab3d ${fname}.in | tee ${fname}.out