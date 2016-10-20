noinst_PROGRAMS= \
unstructured_three_d_fluid \
unstructured_original \
unstructured_lgr_trac \
unstructured_lgr_pressureimposed

#------------------------------------------------------------------------
unstructured_three_d_fluid_SOURCES = unstructured_three_d_fluid.cc

unstructured_three_d_fluid_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------
unstructured_lgr_trac_SOURCES = unstructured_lgr_trac.cc

unstructured_lgr_trac_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


#------------------------------------------------------------------------
unstructured_original_SOURCES = unstructured_original.cc

unstructured_original_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------
unstructured_lgr_pressureimposed_SOURCES = unstructured_lgr_pressureimposed.cc

unstructured_lgr_pressureimposed_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

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
