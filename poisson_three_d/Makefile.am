noinst_PROGRAMS= \
poisson_3d \
poisson_3d_no_bpf \
two_d_unsteady_heat \
two_d_unsteady_heat_t_adapt
#------------------------------------------------------------------------
poisson_3d_SOURCES = poisson_3d.cc

poisson_3d_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------

poisson_3d_no_bpf_SOURCES = poisson_3d_no_bpf.cc

poisson_3d_no_bpf_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------

two_d_unsteady_heat_SOURCES = two_d_unsteady_heat.cc

two_d_unsteady_heat_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------

two_d_unsteady_heat_t_adapt_SOURCES = two_d_unsteady_heat_t_adapt.cc

two_d_unsteady_heat_t_adapt_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -lunsteady_heat -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#------------------------------------------------------------------------
# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES = -I@includedir@  

# Extra actions for clean and distclean
distclean-local:
	touch Validation
	rm -r -f Validation
	rm -f *~
clean-local:
	touch Validation
	rm -r -f Validation
	rm -f *~

# Print top_builddir
print-top_builddir: 
	echo $(top_builddir)
