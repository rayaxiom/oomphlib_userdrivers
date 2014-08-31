noinst_PROGRAMS= \
lin_elas_qele \
lin_elas_lele

#------------------------------------------------------------------------
lin_elas_qele_SOURCES = lin_elas_qele.cc

lin_elas_qele_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -llinear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------


lin_elas_lele_SOURCES = lin_elas_lele.cc

lin_elas_lele_LDADD = -L@libdir@ -lnavier_stokes -lpoisson -llinear_elasticity -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
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
