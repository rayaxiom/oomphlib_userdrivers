noinst_PROGRAMS= \
annular_wedge_threed \
aw3d_mod \
annular_wedge_threed_bk

#------------------------------------------------------------------------
annular_wedge_threed_SOURCES = annular_wedge_threed.cc

annular_wedge_threed_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

aw3d_mod_SOURCES = aw3d_mod.cc

aw3d_mod_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


annular_wedge_threed_bk_SOURCES = annular_wedge_threed_bk.cc

annular_wedge_threed_bk_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

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
