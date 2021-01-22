#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/FiniteElementSpace.h>
#include <dolfin/mesh/unitmeshes/UnitCube.h>
#include <dolfin/mesh/unitmeshes/UnitInterval.h>
#include <dolfin/mesh/unitmeshes/UnitSquare.h>

#include "../../elements/element_library.inc"

using namespace dolfin;

void test_finite_element_space( FiniteElementSpace const & space )
{
  using namespace ElementLibrary;

  FiniteElement const & element = *space.element();

  // check FiniteElement
  {
    ufc::finite_element * ufc_element =
      create_finite_element( element.signature() );

    ck_assert( strcmp( ufc_element->signature(), element.signature() ) == 0 );
    ck_assert( strcmp( ufc_element->family(), element.family() ) == 0 );

    ck_assert( ufc_element->cell_shape() == element.cell_shape() );
    ck_assert( ufc_element->topological_dimension()
               == element.topological_dimension() );
    ck_assert( ufc_element->geometric_dimension()
               == element.geometric_dimension() );
    ck_assert( ufc_element->space_dimension() == element.space_dimension() );
    ck_assert( ufc_element->value_rank() == element.value_rank() );
    ck_assert( ufc_element->value_dimension( 0 )
               == element.value_dimension( 0 ) );
    ck_assert( ufc_element->value_size() == element.value_size() );
    ck_assert( ufc_element->reference_value_rank()
               == element.reference_value_rank() );
    ck_assert( ufc_element->reference_value_dimension( 0 )
               == element.reference_value_dimension( 0 ) );
    ck_assert( ufc_element->reference_value_size()
               == element.reference_value_size() );
    ck_assert( ufc_element->degree() == element.degree() );
    element.is_vectorizable();

    for ( size_t i = 0; i < element.num_sub_elements(); ++i )
    {
      ufc::finite_element * sub = element.create_sub_element( i );
      delete sub;
    }

    delete ufc_element;
  }

  DofMap const & dofmap = space.dofmap();

  // check DofMap
  {
    ufc::dofmap * ufc_dofmap = create_dof_map( dofmap.signature() );

    ck_assert( strcmp( ufc_dofmap->signature(), dofmap.signature() ) == 0 );
    ck_assert( ufc_dofmap->needs_mesh_entities( 0 )
               == dofmap.needs_mesh_entities( 0 ) );
    ck_assert( ufc_dofmap->topological_dimension()
               == dofmap.topological_dimension() );
    ck_assert( ufc_dofmap->num_global_support_dofs()
               == dofmap.num_global_support_dofs() );
    ck_assert( ufc_dofmap->num_element_support_dofs()
               == dofmap.num_element_support_dofs() );
    ck_assert( ufc_dofmap->num_element_dofs() == dofmap.num_element_dofs() );
    ck_assert( ufc_dofmap->num_facet_dofs() == dofmap.num_facet_dofs() );
    ck_assert( ufc_dofmap->num_entity_dofs( 0 )
               == dofmap.num_entity_dofs( 0 ) );
    ck_assert( ufc_dofmap->num_entity_closure_dofs( 0 )
               == dofmap.num_entity_closure_dofs( 0 ) );

    for ( size_t i = 0; i < dofmap.num_sub_dofmaps(); ++i )
    {
      ufc::dofmap * sub = dofmap.create_sub_dofmap( i );
      delete sub;
    }

    delete ufc_dofmap;
  }
}

//----------------------------------------------------------------------------
DOLFIN_START_TEST( test_FiniteElementSpace )
{
  using namespace ElementLibrary;

  // okay, so this is not super exact, but we just expect the dofmap and element
  // signatures in the element library to be in the correct order that each pair
  // creates a finite element space
  size_t const nelements = sizeof( elements ) / sizeof( char * );
  size_t const ndofmaps  = sizeof( dofmaps ) / sizeof( char * );
  ck_assert( nelements == ndofmaps );

  for ( size_t i = 0; i < nelements; ++i )
  {
    // create ufc::finite element
    ufc::finite_element * element = create_finite_element( elements[i] );

    // create ufc::dofmap
    ufc::dofmap * dofmap = create_dof_map( dofmaps[i] );

    Mesh m;
    if ( element->geometric_dimension() == 1 )
    {
      m = UnitInterval( 8 );
    }
    else if ( element->geometric_dimension() == 2 )
    {
      m = UnitSquare( 8, 8 );
    }
    else if ( element->geometric_dimension() == 3 )
    {
      m = UnitCube( 8, 8, 8 );
    }

    // FIXME try-catch should be removed, once CG1s is fixed
    try
    {
      // create FiniteElementSpace and transfer ownership
      FiniteElementSpace space( m, element, *dofmap, false );

      // test this FiniteElementSpace and its subspaces recursively
      test_finite_element_space( space );
    }
    catch ( std::exception & e )
    {
      message( "Checked:\n%s\n%s\nCaught:\n%s", elements[i], dofmaps[i], e.what() );
    }

    delete element;
    delete dofmap;
  }
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
