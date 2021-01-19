#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/DofMap.h>

#include <dolfin/mesh/UnitInterval.h>
#include <dolfin/mesh/UnitSquare.h>
#include <dolfin/mesh/UnitCube.h>

#include "../../elements/element_library.inc"

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_DofMap )
{
  using namespace ElementLibrary;

  for ( char const * dofmap_str : dofmaps )
  {
    ufc::dofmap * ufc_dofmap = create_dof_map( dofmap_str );
    size_t dim = ufc_dofmap->topological_dimension();

    Mesh m;
    if ( dim == 1 )
    {
      m = UnitInterval( 10 );
    }
    else if ( dim == 2 )
    {
      m = UnitSquare( 10, 10 );
    }
    else if ( dim == 3 )
    {
      m = UnitCube( 10, 10, 10 );
    }

    DofMap dofmap( m, *ufc_dofmap, false );
    ck_assert( strcmp( ufc_dofmap->signature(), dofmap.signature() ) == 0 );
    ck_assert( ufc_dofmap->needs_mesh_entities(0) == dofmap. needs_mesh_entities(0) );
    ck_assert( ufc_dofmap->topological_dimension() == dofmap.topological_dimension() );
    ck_assert( ufc_dofmap->num_global_support_dofs() == dofmap.num_global_support_dofs() );
    ck_assert( ufc_dofmap->num_element_support_dofs() == dofmap.num_element_support_dofs() );
    ck_assert( ufc_dofmap->num_element_dofs() == dofmap.num_element_dofs() );
    ck_assert( ufc_dofmap->num_facet_dofs() == dofmap.num_facet_dofs() );
    ck_assert( ufc_dofmap->num_entity_dofs(0) == dofmap.num_entity_dofs(0) );
    ck_assert( ufc_dofmap->num_entity_closure_dofs(0) == dofmap.num_entity_closure_dofs(0) );

    for ( size_t i = 0; i < dofmap.num_sub_dofmaps; ++i )
    {
     ufc::dofmap * sub = dofmap.create_sub_dofmap( i );
     delete sub;
    }
    delete ufc_dofmap;
  }
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
