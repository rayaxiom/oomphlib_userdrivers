noinst_PROGRAMS= \
sq_lgr \
rayleigh_channel \
rayleigh_traction_channel \
three_d_fp \
unstructured_three_d_fluid

#------------------------------------------------------------------------
sq_lgr_SOURCES = sq_lgr.cc

sq_lgr_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
rayleigh_channel_SOURCES = rayleigh_channel.cc

rayleigh_channel_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
rayleigh_traction_channel_SOURCES = rayleigh_traction_channel.cc

rayleigh_traction_channel_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
three_d_fp_SOURCES = three_d_fp.cc

three_d_fp_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
unstructured_three_d_fluid_SOURCES = unstructured_three_d_fluid.cc

unstructured_three_d_fluid_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
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
