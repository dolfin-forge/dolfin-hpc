SUBDIRS = common \
          fem \
          function \
          io \
          la \
          log \
          main \
          mesh \
          parameter

TESTS = common/exe \
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

EXTRA_DIST= dolfin_tests.h.in\
	    elements/element_library.inc\
	    elements/ffc_Brezzi_Douglas_Marini_1_2d.h\
	    elements/ffc_Discontinuous_Lagrange_0_1d.h\
	    elements/ffc_Discontinuous_Lagrange_0_2d.h\
	    elements/ffc_Discontinuous_Lagrange_0_2dVector.h\
	    elements/ffc_Discontinuous_Lagrange_0_3d.h\
	    elements/ffc_Discontinuous_Lagrange_0_3dVector.h\
	    elements/ffc_Discontinuous_Lagrange_1_1d.h\
	    elements/ffc_Discontinuous_Lagrange_1_2d.h\
	    elements/ffc_Discontinuous_Lagrange_1_2dVector.h\
	    elements/ffc_Discontinuous_Lagrange_1_3d.h\
	    elements/ffc_Discontinuous_Lagrange_1_3dVector.h\
	    elements/ffc_Discontinuous_Lagrange_2_1d.h\
	    elements/ffc_Discontinuous_Lagrange_2_2d.h\
	    elements/ffc_Discontinuous_Lagrange_2_2dVector.h\
	    elements/ffc_Discontinuous_Lagrange_2_3d.h\
	    elements/ffc_Discontinuous_Lagrange_2_3dVector.h\
	    elements/ffc_Lagrange_1_1d.h\
	    elements/ffc_Lagrange_1_2d.h\
	    elements/ffc_Lagrange_1_2dVector.h\
	    elements/ffc_Lagrange_1_3d.h\
	    elements/ffc_Lagrange_1_3dVector.h\
	    elements/ffc_Lagrange_2_1d.h\
	    elements/ffc_Lagrange_2_2d.h\
	    elements/ffc_Lagrange_2_2dVector.h\
	    elements/ffc_Lagrange_2_3d.h\
	    elements/ffc_Lagrange_2_3dVector.h
