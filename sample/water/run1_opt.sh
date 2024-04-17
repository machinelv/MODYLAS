export OMP_NUM_THREADS=1
mpirun -np 8 ./modylas water_opt > output
