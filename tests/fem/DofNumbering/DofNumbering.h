#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/DofNumbering.h>
#include <dolfin/fem/UFCCell.h>
#include <dolfin/mesh/celltypes/CellType.h>
#include <dolfin/mesh/entities/iterators/CellIterator.h>
#include <dolfin/mesh/EuclideanSpace.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Space.h>
#include <dolfin/mesh/unitmeshes/UnitCube.h>

#include <set>

#include "../../elements/element_library.inc"

using namespace dolfin;

//-----------------------------------------------------------------------------
DOLFIN_START_TEST( test_DofNumbering )
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
    ufc::dofmap *         dofmap  = create_dof_map( dofmaps[i] );
    ufc::finite_element * element = create_finite_element( elements[i] );
    CellType *            cell    = CellType::create( element->cell_shape() );

    Mesh refcell( *cell, EuclideanSpace( cell->dim() ) );
    cell->create_reference_cell( refcell );

    // FIXME try-catch should be removed, once CG1s is fixed
    try
    {
      DofNumbering * numbering = DofNumbering::create( refcell, *dofmap );

      numbering->init( refcell, *dofmap );
      numbering->build();
      bool           shared     = numbering->is_shared( 0 );
      bool           ghost      = numbering->is_ghost( 0 );
      size_t         offset     = numbering->offset();
      size_t         size       = numbering->size();
      size_t const * block      = numbering->block();
      size_t         block_size = numbering->block_size();

      delete numbering;
    }
    catch ( std::exception & e )
    {
      message( "Checked:\n%s\n%s\nCaught:\n%s", elements[i], dofmaps[i], e.what() );
    }

    delete dofmap;
    delete element;
    delete cell;
  }
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

// Verify that DG0vNumbering::tabulate_dofs returns dofs[k] = vs*global_cell_index + k
// and that all indices are unique across all owned cells.
DOLFIN_START_TEST( test_DofNumbering_DG0v )
{
  UnitCube mesh( 2, 2, 2 );

  ffc_discontinuous_lagrange_0_3dvector_dofmap_1 ufc_dofmap;
  DofNumbering * numbering = DofNumbering::create( mesh, ufc_dofmap );
  ck_assert_msg( numbering != nullptr, "DG0vNumbering not created" );
  numbering->init( mesh, ufc_dofmap );
  numbering->build();

  static const size_t vs = 3;
  std::set< size_t > seen_dofs;
  size_t nchecked = 0;

  for ( CellIterator c( mesh ); !c.end(); ++c )
  {
    if ( c->is_ghost() ) continue;
    UFCCell ufc_cell( *c );
    ufc_cell.update( *c );

    size_t dofs[vs];
    numbering->tabulate_dofs( dofs, ufc_cell, *c );

    size_t const global_idx = c->global_index();
    for ( size_t k = 0; k < vs; ++k )
    {
      ck_assert_msg( dofs[k] == vs * global_idx + k,
        "DG0v dofs[%zu]=%zu expected %zu (global_idx=%zu)",
        k, dofs[k], vs * global_idx + k, global_idx );
      ck_assert_msg( seen_dofs.count( dofs[k] ) == 0, "DG0v dof index collision" );
      seen_dofs.insert( dofs[k] );
    }
    ++nchecked;
  }

  ck_assert_msg( nchecked > 0, "No owned cells found" );
  delete numbering;
}
DOLFIN_END_TEST
//-----------------------------------------------------------------------------

#endif
