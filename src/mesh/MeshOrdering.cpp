// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-01-30
// Last changed: 2008-06-17

#include <dolfin/log/log.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshOrdering.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
void MeshOrdering::order( Mesh & mesh )
{
	message( 1, "MeshOrdering: Ordering mesh entities..." );

	// Get cell type
	const CellType & cell_type = mesh.type();

	// Iterate over all cells and order the mesh entities locally
	// NOTE: ensure that boundary meshes and empty meshes are not reordered
	if ( mesh.topology().connectivity( cell_type.dim() ) )
	{
		uint const num_cells = mesh.topology().size( cell_type.dim() );
		for ( uint i = 0; i < num_cells; ++i )
		{
			cell_type.order_entities( mesh.topology(), i );
		}
	}

	// mesh._ordered = true;
}
//-----------------------------------------------------------------------------

}
