LDFLAGS = $(XML_LDFLAGS)
SUBDIRS = common elements fem function graph io\
	la log main math mesh parameter quadrature 

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = dolfin.pc

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
		parameter/libdolfin_parameter.la\
		quadrature/libdolfin_quadrature.la
		parameter/libdolfin_parameter.la

