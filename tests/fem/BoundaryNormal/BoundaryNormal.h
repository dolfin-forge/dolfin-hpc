#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/Elements.h>
#include <dolfin/fem/NodeNormal.h>
#include <dolfin/mesh/Mesh.h>

#include "../../elements/element_library.inc"

using namespace dolfin;

//-----------------------------------------------------------------------------
void test( std::string file )
{
  using namespace ElementLibrary;

  Mesh mesh( file );

  std::string const element_str =   "VectorElement(FiniteElement('Lagrange', "
                                  + mesh.type().str() + ", 1), dim="
                                  + std::to_string( mesh.geometry_dimension() )
                                  + ")";
  std::string const dofmap_str  = "FFC dofmap for " + element_str;

  ufc::finite_element * element = create_finite_element( element_str.c_str() );
  ufc::dofmap *         dofmap  = create_dof_map( dofmap_str.c_str() );

  FiniteElementSpace Vh( mesh, element, *dofmap, true );

  NodeNormal nn( mesh );
  nn.init( Vh );
  nn.compute();
}
//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_NodeNormal )
{
  test( mesh_file( "cylinder.bin" ) );
  test( mesh_file( "aneurysm.bin" ) );
  test( mesh_file( "sphere.bin" ) );
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
