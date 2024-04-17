bin_PROGRAMS = modylas-text2bin modylas-mdtrj2xyz
modylas_SOURCES = version.f fmodules.f main_f90.f parse_f.f app_f90.f \
				  md_charmm_f90.f md_fmm_f90.f md_direct_f90.f \
				  MAP_module.f CMAP.f md_pmewald_f90.f ffte.f pzfft3dv.f \
				  param.h nve_integrate.f k_energy.f cell_edge.f \
				  mt19937ar.f mpitool.f domain_div.f assign2subcell.f \
				  nvt_integrate.f shake_rattle_roll.f \
				  npt_a_integrate.f p_calc.f nhchain.f nhc_a.f cell_convert.f \
				  npt_z_integrate.f nhc_z.f \
				  npt_pr_integrate.f nhc_pr.f \
				  array.f hash.f parse.f parse_input.f mdff.f \
				  diagonal.f matrix_inverse.f \
				  comm.f comm_2.f comm_3.f comm_fmm.f comm_qgrid.f \
				  opt_integrate.f  position_constrain.f \
                          calc_center_of_mass.f shake_rattle_com.f \
                          md_analyze_enectric_field.f md_separate_force.f \
                          md_localpres.f md_dlvo.f
modylas_text2bin_SOURCES = version.f fmodules.f MAP_module.f \
                           modylas-text2bin.f \
                           mpitool.f array.f hash.f parse.f \
                           parse_f.f parse_input.f mdff.f
modylas_mdtrj2xyz_SOURCES = version.f fmodules.f MAP_module.f \
                           modylas-mdtrj2xyz.f \
                           mpitool.f array.f hash.f parse.f \
                           parse_f.f parse_input.f mdff.f
if HAVE_MPI
bin_PROGRAMS += modylas
endif
F77LINK = $(FCLINK)
parse_input.o: parse.o
parse.o: hash.o
hash.o: array.o
fmodules.o: fmodules.f
	$(FCCOMPILE) -c -o $@ $<
%.o: %.f fmodules.o
	$(FCCOMPILE) -c -o $@ $<
%.obj: %.f fmodules.o
	$(FCCOMPILE) -c -o $@ `$(CYGPATH_W) '$<'`
html:
	doxygen
pdf: html
	cd latex; make; mv refman.pdf ../
clean:
	${RM} -f *.o *~ *.bak *.mod *.MOD
