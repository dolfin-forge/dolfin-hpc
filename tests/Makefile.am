SUBDIRS = \
  common \
  ufl \
  log \
  la \
  elements \
  fem \
  function \
  io \
  main \
  mesh \
  parameter

TESTS = \
  common/exe \
  ufl/exe \
  log/exe \
  la/exe \
  elements/exe \
  fem/exe \
  function/exe \
  io/exe \
  main/exe \
  mesh/exe \
  parameter/exe

dolfin_tests.h:	Makefile
	$(SED) -e 's|#abs_top_srcdir#|@abs_top_srcdir@|g' \
		@top_srcdir@/tests/dolfin_tests.h.in > dolfin_tests.h

BUILT_SOURCES = dolfin_tests.h
DISTCLEANFILES= dolfin_tests.h

EXTRA_DIST= dolfin_tests.h.in

