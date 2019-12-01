
#  Written by C. Froese Fischer and Y. Qiu
#  Modified and extended by O. Zatsarinny
#  
#  Modified by C. Froese Fischer
#  January, 2006

.SUFFIXES: .f90 .mod

MODDIR = .

#FC = gfortran
#FC_FLAGS = -g
#LAPXK_DIR = ??

install:  spmchf


My_LIBS= spline_param.o spline_grid.o spline_galerkin.o spline_hl.o \
         spline_integrals.o spline_moments.o \
 	 mchf_atomic_state.o orbitals.o memory_use.o block_param.o \
         mchf_inout.o angular.o av_energy.o  mchf_param.o cl_arguments.o\
         gauss.o vbsplvd.o mrm.o bxv.o hlm.o hlc.o azl.o\
         mrk_cell.o rk.o density.o \
         yval.o minty.o bhwf.o rk_moments.o vinty.o moments.o \
         routines.o bhl.o bvmv.o quadr.o  el_nl.o bvalu2.o \
         lapack_interface.o  plot_bsw.o 

My_APPS= spmchf.o read_mchf_param.o get_case.o get_atom.o intgrl.o \
         get_spline_param.o define_spline.o mkgrid.o get_estimates.o \
         read_bsw.o unpacki.o update.o ecore.o solve_mchf.o \
         dvdson_param.o diag.o updatc.o \
         scf_matrix.o grad.o summry.o apply_orthogonality.o rotate.o \
         scf_nr.o  some_nr.o iniest2.o eig_out.o pack.o prepV.o \
         refine_grid.o  lagrange.o read_cl_param.o

spmchf: $(My_LIBS) $(My_APPS) $(dvdson)
	$(FC) -o spmchf $(FC_LD)  $(My_LIBS) $(My_APPS) $(dvdson) \
        -L $(LAPACK_DIR) -llapack -lblas
	
.f90.o:
	$(FC) -c $(FC_FLAGS) $< -I $(MODDIR) -o $@

clean:
	-rm -f *.o core  *.mod


