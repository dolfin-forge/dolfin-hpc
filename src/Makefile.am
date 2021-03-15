LDFLAGS =
SUBDIRS = adaptivity common evolution fem function io\
	        la log main math mesh parameter quadrature insitu

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = dolfin.pc

lib_LTLIBRARIES = libdolfin.la
libdolfin_la_SOURCES =
libdolfin_la_LIBADD = adaptivity/libdolfin_adaptivity.la\
                      common/libdolfin_common.la\
                      evolution/libdolfin_evolution.la\
                      fem/libdolfin_fem.la\
                      function/libdolfin_function.la\
                      io/libdolfin_io.la\
                      insitu/libdolfin_insitu.la\
                      la/libdolfin_la.la\
                      log/libdolfin_log.la\
                      main/libdolfin_main.la\
                      math/libdolfin_math.la\
                      mesh/libdolfin_mesh.la\
                      parameter/libdolfin_parameter.la\
                      quadrature/libdolfin_quadrature.la

