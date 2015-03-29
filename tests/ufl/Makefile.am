check_PROGRAMS = exe
exe_SOURCES = test.cpp
exe_CPPFLAGS= -I@top_srcdir@/include @CHECK_CFLAGS@ $(CPPFLAGS)
exe_LDFLAGS = $(top_builddir)/src/.libs/libdolfin.a $(LDFLAGS) $(LIBS) @CHECK_LIBS@


#SUBDIRS = Algebra Argument array Cell CellSurfaceArea CellVolume Circumradius Class Coefficient Conditional Differentiation Domain ElementList EnrichedElement Equation Expression FacetArea FacetNormal Family FiniteElement FiniteElementBase Form GeometricQuantity Index Indexed IndexSum Integral List MixedElement Object QuadratureScheme repr RestrictedElement Space SpatialCoordinate TensorElement Tensors tuple Tuple type Variable VectorElement
