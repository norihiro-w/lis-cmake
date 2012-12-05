EXTRA_DIST = precision.inc
noinst_LTLIBRARIES = libsaamg.la

if ENABLE_MPI
  libsaamg_la_SOURCES = lis_m_data_type_AMGCG.F90 \
  lis_m_sort.F90 \
  lis_m_queue.F90 \
  lis_m_hash.F90 \
  lis_m_count_time_mod.F90 \
  lis_m_solver_SR2.F90 \
  lis_m_solver_Gnum.F90 \
  lis_m_aggregate.F90 \
  lis_m_data_creation_AMGCG.F90 \
  lis_m_solver_AMGCG.F90 \
  lis_m_finit.F90
else
  libsaamg_la_SOURCES = lis_s_data_type_AMGCG.F90 \
  lis_s_sort.F90 \
  lis_s_queue.F90 \
  lis_s_aggregate.F90 \
  lis_s_data_creation_AMGCG.F90 \
  lis_s_solver_AMGCG.F90 \
  lis_s_finit.F90
endif

CLEANFILES = *.mod
AM_CPPFLAGS = -I$(top_srcdir)/include
FCCOMPILE = $(FC) $(AMDEFS) $(AM_FCFLAGS) $(FCFLAGS)
LTFCCOMPILE = $(LIBTOOL) --mode=compile $(FC) $(AMDEFS) $(AM_FCFLAGS) $(FCFLAGS)
FCLD = $(FC)
FCLINK = $(LIBTOOL) --mode=link $(FCLD) $(AM_LDFLAGS) $(LDFLAGS) -o $@
PPFCCOMPILE = $(FC) $(AMDEFS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
        $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)
LTPPFCCOMPILE = $(LIBTOOL) --mode=compile $(FC) $(AMDEFS) $(DEFS) \
        $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
        $(AM_FCFLAGS) $(FCFLAGS)
PPF77COMPILE = $(F77) $(AMDEFS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FFLAGS) $(FFLAGS)
LTPPF77COMPILE = $(LIBTOOL) --mode=compile --tag=F77 $(F77) $(AMDEFS) $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_FFLAGS) $(FFLAGS)
COMPILE = $(CC) $(AMDEFS) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) \
	$(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
LTCOMPILE = $(LIBTOOL) --mode=compile --tag=CC $(CC) $(AMDEFS) $(DEFS) \
	$(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \
	$(AM_CFLAGS) $(CFLAGS)

.NOTPARALLEL:
.NO_PARALLEL:
