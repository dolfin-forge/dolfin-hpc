// Copyright (C) 2006-2007 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin.h>
#include <dolfin/mesh/MeshSmoothing.h>

using namespace dolfin;

int main()
{
	for ( unsigned int j = 0; j < 1; ++j )
	{
		// Create mesh of unit square
		UnitSquare mesh( 15, 15 );
		// Mesh mesh( "model.stl" );

		File f_mesh( ( j == 0 ) ? "simple.pvd" : "rivara.pvd" );

		f_mesh << mesh;

		// Uniform refinement
		mesh.refine();
		f_mesh << mesh;

		// Refine mesh close to x = (0.5, 0.5)
		Point p( 0.5, 0.5 );
		for ( unsigned int i = 0; i < 1; ++i )
		{
			// Mark cells for refinement
			MeshValues< bool, Cell > cell_markers( mesh );
			cell_markers = false;
			for ( CellIterator c( mesh ); !c.end(); ++c )
			{
				if ( c->midpoint().dist( p ) < 0.1 )
					cell_markers( *c ) = true;
			}

			// Refine mesh
			if ( j == 0 )
				mesh.refine();
			else
				RivaraRefinement::refine( mesh, cell_markers );

			// // Smooth mesh
			MeshSmoothing::smooth( mesh );

			f_mesh << mesh;
		}
	}
	return 0;
}
