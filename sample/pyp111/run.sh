export OMP_NUM_THREADS=1
mpirun -np 8 ./modylas pyp111 > output
