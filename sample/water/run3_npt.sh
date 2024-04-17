cp water.mddef_for_npt water_npt.mddef
ln -s water_opt.mdff water_npt.mdff
ln -s water_nvt.restart.bin water_npt.mdxyz.bin

export OMP_NUM_THREADS=1
mpirun -np 8 ./modylas water_npt > output_npt
