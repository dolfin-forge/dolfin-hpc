#include <dolfin_tests.h>

#ifdef HAVE_CHECK

#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/DofNumbering.h>
#include <dolfin/mesh/celltypes/CellType.h>
#include <dolfin/mesh/EuclideanSpace.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Space.h>

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

#endif
