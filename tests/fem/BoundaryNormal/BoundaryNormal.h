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

  ufc::finite_element * element =
    create_finite_element( "FiniteElement('Lagrange', triangle, 1)" );

  ufc::dofmap * dofmap = create_dof_map(
    "FFC dofmap for FiniteElement('Lagrange', triangle, 1)" );

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
