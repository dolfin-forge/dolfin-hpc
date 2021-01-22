// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/Box.h>

#include <dolfin/mesh/EuclideanSpace.h>
#include <dolfin/mesh/MeshEditor.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

Box::Box( real   a,
          real   b,
          real   c,
          real   d,
          real   e,
          real   f,
          size_t nx,
          size_t ny,
          size_t nz )
  : Mesh( *CellType::create( CellType::tetrahedron ), EuclideanSpace( 3 ) )
{

  if ( nx < 1 || ny < 1 || nz < 1 )
  {
    error( "Size of box must be at least 1 in each dimension." );
  }

  rename( "mesh", "Mesh of the cuboid (a,b) x (c,d) x (e,f)" );

  // Open mesh for editing
  MeshEditor editor( *this, this->type(), this->space() );

  // Create vertices
  editor.init_vertices( ( nx + 1 ) * ( ny + 1 ) * ( nz + 1 ) );
  size_t vertex = 0;
  real   x[3]   = { 0.0 };
  for ( size_t iz = 0; iz <= nz; ++iz )
  {
    x[2] =
      e + ( static_cast< real >( iz ) ) * ( f - e ) / static_cast< real >( nz );
    for ( size_t iy = 0; iy <= ny; ++iy )
    {
      x[1] =
        c
        + ( static_cast< real >( iy ) ) * ( d - c ) / static_cast< real >( ny );
      for ( size_t ix = 0; ix <= nx; ++ix )
      {
        x[0] = a
               + ( static_cast< real >( ix ) ) * ( b - a )
                   / static_cast< real >( nx );
        editor.add_vertex( vertex++, x );
      }
    }
  }

  // Create tetrahedra
  editor.init_cells( 6 * nx * ny * nz );
  size_t cell = 0;
  for ( size_t iz = 0; iz < nz; iz++ )
  {
    for ( size_t iy = 0; iy < ny; iy++ )
    {
      for ( size_t ix = 0; ix < nx; ix++ )
      {
        size_t const v0 = iz * ( nx + 1 ) * ( ny + 1 ) + iy * ( nx + 1 ) + ix;
        size_t const v1 = v0 + 1;
        size_t const v2 = v0 + ( nx + 1 );
        size_t const v3 = v1 + ( nx + 1 );
        size_t const v4 = v0 + ( nx + 1 ) * ( ny + 1 );
        size_t const v5 = v1 + ( nx + 1 ) * ( ny + 1 );
        size_t const v6 = v2 + ( nx + 1 ) * ( ny + 1 );
        size_t const v7 = v3 + ( nx + 1 ) * ( ny + 1 );

        size_t const connectivity[24] = { v0, v1, v3, v7, v0, v1, v7, v5,
                                          v0, v5, v7, v4, v0, v3, v2, v7,
                                          v0, v6, v4, v7, v0, v2, v6, v7 };

        editor.add_cell( cell++, &connectivity[0] );
        editor.add_cell( cell++, &connectivity[4] );
        editor.add_cell( cell++, &connectivity[8] );
        editor.add_cell( cell++, &connectivity[12] );
        editor.add_cell( cell++, &connectivity[16] );
        editor.add_cell( cell++, &connectivity[20] );
      }
    }
  }

  // Close mesh editor
  editor.close();
}

//-----------------------------------------------------------------------------

} // namespace dolfin
