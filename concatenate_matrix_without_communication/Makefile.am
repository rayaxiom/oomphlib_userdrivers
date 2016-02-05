noinst_PROGRAMS= \
mat_cat \
test_parallel

#------------------------------------------------------------------------
mat_cat_SOURCES = mat_cat.cc

mat_cat_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------

#------------------------------------------------------------------------
test_parallel_SOURCES = test_parallel.cc

test_parallel_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
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
