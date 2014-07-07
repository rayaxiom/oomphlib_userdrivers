noinst_PROGRAMS= \
qcirc2d \
qcirc2d_working_but_not_right

#------------------------------------------------------------------------
qcirc2d_SOURCES = qcirc2d.cc

qcirc2d_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
qcirc2d_working_but_not_right_SOURCES = qcirc2d_working_but_not_right.cc

qcirc2d_working_but_not_right_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
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
