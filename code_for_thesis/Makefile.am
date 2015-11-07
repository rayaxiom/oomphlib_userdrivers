noinst_PROGRAMS= \
two_d_ns \
driven_cavity

#------------------------------------------------------------------------
two_d_ns_SOURCES = two_d_ns.cc

two_d_ns_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
#------------------------------------------------------------------------
driven_cavity_SOURCES = driven_cavity.cc

driven_cavity_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
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
