LDFLAGS = $(XML_LDFLAGS)
SUBDIRS = common elements fem function graph io\
	la log main math mesh mf nls ode parameter\
	pde quadrature

lib_LTLIBRARIES = libdolfin.la
libdolfin_la_SOURCES = 
libdolfin_la_LIBADD = \
		common/libdolfin_common.la\
		elements/libdolfin_elements.la\
		fem/libdolfin_fem.la\
		function/libdolfin_function.la\
		graph/libdolfin_graph.la\
		io/libdolfin_io.la\
		la/libdolfin_la.la\
		log/libdolfin_log.la\
		main/libdolfin_main.la\
		math/libdolfin_math.la\
		mesh/libdolfin_mesh.la\
		mf/libdolfin_mf.la\
		nls/libdolfin_nls.la\
		ode/libdolfin_ode.la\
		parameter/libdolfin_parameter.la\
		pde/libdolfin_pde.la\
		quadrature/libdolfin_quadrature.la

pkginclude_HEADERS = config/dolfin_config.h

