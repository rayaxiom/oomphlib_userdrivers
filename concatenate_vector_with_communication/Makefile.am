noinst_PROGRAMS= \
veccat \
veccat_without_comm \
vecsplit_without_comm 

#------------------------------------------------------------------------
veccat_SOURCES = veccat.cc

veccat_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
#------------------------------------------------------------------------
veccat_without_comm_SOURCES = veccat_without_comm.cc

veccat_without_comm_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------
vecsplit_without_comm_SOURCES = vecsplit_without_comm.cc

vecsplit_without_comm_LDADD = -L@libdir@ -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


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
