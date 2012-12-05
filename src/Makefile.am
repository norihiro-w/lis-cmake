SUBDIRS = matrix vector matvec precon solver esolver system precision
if ENABLE_FORTRAN
  SUBDIRS += fortran
endif
if ENABLE_SAAMG
  SUBDIRS += fortran
  SUBDIRS += fortran/amg
endif

lib_LTLIBRARIES  = liblis@LIBSUFFIX@.la
liblis@LIBSUFFIX@_la_LIBADD = matrix/libmatrix.la \
                   vector/libvector.la \
                   matvec/libmatvec.la \
                   precon/libprecon.la \
                   solver/libsolver.la \
                   esolver/libesolver.la \
                   system/libsystem.la \
                   precision/libprecision.la 
if ENABLE_FORTRAN
  liblis@LIBSUFFIX@_la_LIBADD += fortran/libfortran.la
endif
if ENABLE_SAAMG
  liblis@LIBSUFFIX@_la_LIBADD += fortran/amg/libsaamg.la
  liblis@LIBSUFFIX@_la_LIBADD += fortran/libfortran.la
endif
liblis@LIBSUFFIX@_la_SOURCES = 
