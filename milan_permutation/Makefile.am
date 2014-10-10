noinst_PROGRAMS= \
poisson_3d \
poisson_3d_no_bpf \
mesh_from_tetgen_poisson

#------------------------------------------------------------------------
poisson_3d_SOURCES = poisson_3d.cc
poisson_3d_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------

poisson_3d_no_bpf_SOURCES = poisson_3d_no_bpf.cc
poisson_3d_no_bpf_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------

mesh_from_tetgen_poisson_SOURCES = mesh_from_tetgen_poisson.cc
mesh_from_tetgen_poisson_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)




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
