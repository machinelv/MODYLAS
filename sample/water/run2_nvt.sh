cp water.mddef_for_nvt water_nvt.mddef
ln -s water_opt.mdff water_nvt.mdff
ln -s water_opt.restart.bin water_nvt.mdxyz.bin

export OMP_NUM_THREADS=1
mpirun -np 8 ./modylas water_nvt > output_nvt
