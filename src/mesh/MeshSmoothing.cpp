// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/common/constants.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/MeshOrdering.h>
#include <dolfin/mesh/MeshSmoothing.h>
#include <dolfin/mesh/Vertex.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
void MeshSmoothing::smooth( Mesh & mesh )
{
	// Make sure we have cell-facet connectivity
	mesh.init( mesh.topology().dim(), mesh.topology().dim() - 1 );

	// Make sure the mesh is ordered
	MeshOrdering::order( mesh );
  // mesh.topology().finalize();

	// Mark vertices on the boundary so we may skip them
	BoundaryMesh boundary = mesh.exterior_boundary();
	MeshValues< uint, Vertex > vertex_map( boundary, 0.0 );
	MeshValues< bool, Vertex > on_boundary( mesh, 0 );
	on_boundary = false;
	for ( VertexIterator v( boundary ); !v.end(); ++v )
		on_boundary( ( vertex_map )( *v ), true );

	/*
	   FIXME: new algorithm
	   1) Iterate over shared vertices
	      - Compute neighbors and center off mass, localy
	- Exchange neighbors and center off mass
	- Compute new vertex coordinate on owner
	      - Exchange new coordinates
	   2) Iterate over all shared vertices
	      - Skip shared entities
	*/

	// Iterate over all vertices
	const uint    d = mesh.geometry().dim();
	Array< real > xx( d );
	for ( VertexIterator v( mesh ); !v.end(); ++v )
	{
		// Skip vertices on the boundary
		// FIXME: we also skip shared entities
		if ( on_boundary( *v ) )
			continue;

    if( mesh.is_distributed() && mesh.distdata()[0].is_shared( v->index() ) )
      continue;

		// Get coordinates of vertex
		real *      x = v->x();
		const Point p = v->point();

		// Compute center of mass of neighboring vertices
		for ( uint i = 0; i < d; i++ )
			xx[i] = 0.0;
		uint num_neighbors = 0;
		for ( VertexIterator vn( *v ); !vn.end(); ++vn )
		{
      std::cout << vn->x()[0] << " / " << vn->x()[1]  << std::endl;
			// Skip the vertex itself
			if ( v->index() == vn->index() )
				continue;
			num_neighbors += 1;

			// Compute center of mass
			const real * xn = vn->x();
			for ( uint i = 0; i < d; i++ )
				xx[i] += xn[i];
		}

		for ( uint i = 0; i < d; i++ )
			xx[i] /= static_cast< real >( num_neighbors );

		// Compute closest distance to boundary of star
		real rmin = 0.0;
		for ( CellIterator c( *v ); !c.end(); ++c )
		{
			// Get local number of vertex relative to facet
			const uint local_vertex = c->index( *v );

			// Get normal of corresponding facet
			Point n = c->normal( local_vertex );

			// Get first vertex in facet
			Facet f( mesh, c->entities( mesh.topology().dim() - 1 )[local_vertex] );
			VertexIterator fv( f );

			// Compute length of projection of v - fv onto normal
			const real r = std::abs( n.dot( p - fv->point() ) );
			if ( rmin == 0.0 )
				rmin = r;
			else
				rmin = std::min( rmin, r );
		}

		// Move vertex at most a distance rmin / 2
		real r = 0.0;
		for ( uint i = 0; i < d; i++ )
		{
			const real dx = xx[i] - x[i];
			r += dx * dx;
		}
		r = std::sqrt( r );
		if ( r < DOLFIN_EPS )
			continue;
		rmin = std::min( 0.5 * rmin, r );
		for ( uint i = 0; i < d; i++ )
			x[i] += rmin * ( xx[i] - x[i] ) / r;
	}
}

} // namespace dolfin
//-----------------------------------------------------------------------------
