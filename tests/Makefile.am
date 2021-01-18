SUBDIRS = common \
          elements \
          fem \
          function \
          io \
          la \
          log \
          main \
          mesh \
          parameter

TESTS = common/exe \
        elements/exe \
        fem/exe \
        function/exe \
        io/exe \
        la/exe \
        log/exe \
        main/exe \
        mesh/exe \
        parameter/exe

dolfin_tests.h:	Makefile
	$(SED) -e 's|#abs_top_srcdir#|@abs_top_srcdir@|g' \
		@top_srcdir@/tests/dolfin_tests.h.in > dolfin_tests.h

BUILT_SOURCES = dolfin_tests.h
DISTCLEANFILES= dolfin_tests.h

EXTRA_DIST= dolfin_tests.h.in

