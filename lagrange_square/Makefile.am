noinst_PROGRAMS= \
sq_lgr

#------------------------------------------------------------------------
sq_lgr_SOURCES = sq_lgr.cc

sq_lgr_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
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
	echo $(abs_top_builddir)
