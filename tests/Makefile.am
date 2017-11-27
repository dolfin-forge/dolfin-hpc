SUBDIRS = common \
	  ufl \
	  la \
	  elements \
	  fem \
	  function \
	  io \
	  mesh

TESTS =  common/exe \
	ufl/exe \
	la/exe \
	elements/exe \
	fem/exe \
	function/exe \
	io/exe \
	mesh/exe

tests.h:	Makefile
	$(SED) -e 's|#abs_top_srcdir#|@abs_top_srcdir@|g' \
		@top_srcdir@/tests/tests.h.in > tests.h

BUILT_SOURCES = tests.h
DISTCLEANFILES= tests.h

EXTRA_DIST= tests.h.in
