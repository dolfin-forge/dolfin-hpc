LDFLAGS = $(XML_LDFLAGS)
SUBDIRS = common elements fem function graph io la

lib_LIBRARIES = libdolfin.a
libdolfin_a_SOURCES = 
libdolfin_a_LIBADD = \
		common/libdolfin_common.a\
		elements/libdolfin_elements.a\
		fem/libdolfin_fem.a\
		function/libdolfin_function.a

pkginclude_HEADERS = config/dolfin_config.h