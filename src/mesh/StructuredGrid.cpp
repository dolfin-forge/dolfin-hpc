// Copyright (C) 2015 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/StructuredGrid.h>

#include <dolfin/mesh/EuclideanSpace.h>
#include <dolfin/mesh/MeshEditor.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
StructuredGrid::StructuredGrid( CellType const & type, size_t N )
  : Mesh( type, EuclideanSpace( type.space_dim() ) )
  , bbox_( type.dim() )
  , n_( N )
{
  init( type );
}

//-----------------------------------------------------------------------------
StructuredGrid::StructuredGrid( CellType const & type,
                                size_t           N,
                                BoundingBox      bbox )
  : Mesh( type, EuclideanSpace( type.space_dim() ) )
  , bbox_( bbox )
  , n_( N )

{
  init( type );
}

//-----------------------------------------------------------------------------
void StructuredGrid::init( CellType const & type )
{
  if ( n_ == 0 )
  {
    error( "StructuredGrid::init : provided number of cells per axis is zero" );
  }

  type.disp();

  //
  size_t const tdim = type.dim();
  MeshEditor   editor( *this, type );
  // Number of cells in each direction
  size_t * n = new size_t[tdim];
  for ( size_t i = 0; i < tdim; ++i )
  {
    n[i] = n_;
  } // isotropic

  // Number of vertices in each direction
  size_t m[Space::MAX_DIMENSION] = { 1 };
  size_t num_verts               = 1;
  size_t num_bricks              = 1;
  for ( size_t i = 0; i < tdim; ++i )
  {
    num_verts *= ( n[i] + 1 );
    num_bricks *= n[i];
    m[i] = num_verts;
  }

  //
  size_t * i = new size_t[tdim];
  real *   h = new real[tdim];
  real *   x = new real[tdim];

  // Create vertices
  size_t vertex = 0;
  for ( size_t d = 0; d < tdim; ++d )
  {
    i[d] = 0;
    h[d] = 1.0 / static_cast< real >( n[d] );
    x[d] = bbox_[0][d];
  }
  editor.init_vertices( num_verts );
  while ( vertex < num_verts )
  {
    editor.add_vertex( vertex++, x );
    x[0] = bbox_[0][0] + ( vertex % m[0] ) * h[0];
    for ( size_t d = 1; d < tdim; ++d )
    {
      if ( vertex % m[d - 1] == 0 )
      {
        i[d] = ( i[d] + 1 ) % ( n[d] + 1 );
        x[d] = bbox_[0][d] + i[d] * h[d];
      }
    }
  }

  // Create cellsnnnnn
  size_t       cell = 0;
  size_t const v0   = 0;
  size_t const v1   = 1;
  size_t const v2   = ( n[0] + 1 );
  size_t const v3   = ( n[0] + 1 ) + 1;
  size_t const v4   = v0 + ( n[0] + 1 ) * ( ( tdim > 1 ) ? ( n[1] + 1 ) : 0 );
  size_t const v5   = v1 + ( n[0] + 1 ) * ( ( tdim > 1 ) ? ( n[1] + 1 ) : 0 );
  size_t const v6   = v2 + ( n[0] + 1 ) * ( ( tdim > 1 ) ? ( n[1] + 1 ) : 0 );
  size_t const v7   = v3 + ( n[0] + 1 ) * ( ( tdim > 1 ) ? ( n[1] + 1 ) : 0 );
  switch ( type.cellType() )
  {
    case CellType::point:
    {
      editor.init_cells( num_bricks );
      for ( size_t k = 0; k < n[2]; ++k )
      {
        for ( size_t j = 0; j < n[1]; ++j )
        {
          for ( size_t i = 0; i < n[0]; ++i )
          {
            size_t cv = k * ( ( tdim > 1 ) ? ( n[1] + 1 ) : 0 ) * ( n[0] + 1 )
                        + j * ( n[0] + 1 ) + i;
            editor.add_cell( cell++, &cv );
          }
        }
      }
    }
    break;
    case CellType::interval:
    {
      size_t cv[2] = { v0, v1 };
      editor.init_cells( num_bricks );
      for ( size_t i = 0; i < n[0]; ++i )
      {
        editor.add_cell( cell++, &cv[0] );
        ++cv[0];
        ++cv[1];
      }
    }
    break;
    case CellType::quadrilateral:
    {
      size_t const co[4] = { v0, v1, v3, v2 };
      size_t       cv[4] = { 0 };
      editor.init_cells( num_bricks );
      for ( size_t j = 0; j < n[1]; ++j )
      {
        for ( size_t i = 0; i < n[0]; ++i )
        {
          cv[0] = j * ( n[0] + 1 ) + i;
          cv[1] = cv[0] + co[1];
          cv[2] = cv[0] + co[2];
          cv[3] = cv[0] + co[3];
          editor.add_cell( cell++, &cv[0] );
        }
      }
    }
    break;
    case CellType::hexahedron:
    {
      size_t const co[8] = { v0, v1, v3, v2, v4, v5, v7, v6 };
      size_t       cv[8] = { 0 };
      editor.init_cells( num_bricks );
      for ( size_t k = 0; k < n[2]; ++k )
      {
        for ( size_t j = 0; j < n[1]; ++j )
        {
          for ( size_t i = 0; i < n[0]; ++i )
          {
            cv[0] = k * ( n[1] + 1 ) * ( n[0] + 1 ) + j * ( n[0] + 1 ) + i;
            cv[1] = cv[0] + co[1];
            cv[2] = cv[0] + co[2];
            cv[3] = cv[0] + co[3];
            cv[4] = cv[0] + co[4];
            cv[5] = cv[0] + co[5];
            cv[6] = cv[0] + co[6];
            cv[7] = cv[0] + co[7];
            editor.add_cell( cell++, &cv[0] );
          }
        }
      }
    }
    break;
    case CellType::triangle:
    {
      size_t const co[6] = { v0, v1, v3, v0, v2, v3 };
      size_t       cv[6] = { 0 };
      editor.init_cells( 2 * num_bricks );
      for ( size_t j = 0; j < n[1]; ++j )
      {
        for ( size_t i = 0; i < n[0]; ++i )
        {
          cv[0] = j * ( n[0] + 1 ) + i;
          cv[1] = cv[0] + co[1];
          cv[2] = cv[0] + co[2];
          cv[3] = cv[0] + co[3];
          cv[4] = cv[0] + co[4];
          cv[5] = cv[0] + co[5];
          editor.add_cell( cell++, &cv[0] );
          editor.add_cell( cell++, &cv[3] );
        }
      }
    }
    break;
    case CellType::tetrahedron:
    {
      size_t const co[24] = { v0, v1, v3, v7, v0, v1, v7, v5, v0, v5, v7, v4,
                              v0, v3, v2, v7, v0, v6, v4, v7, v0, v2, v6, v7 };
      size_t       cv[24] = { 0 };
      editor.init_cells( 6 * num_bricks );
      for ( size_t k = 0; k < n[2]; ++k )
      {
        for ( size_t j = 0; j < n[1]; ++j )
        {
          for ( size_t i = 0; i < n[0]; ++i )
          {
            cv[0] = k * ( n[1] + 1 ) * ( n[0] + 1 ) + j * ( n[0] + 1 ) + i;
            for ( size_t c = 1; c < 24; ++c )
            {
              cv[c] = cv[0] + co[c];
            }
            editor.add_cell( cell++, &cv[0] );
            editor.add_cell( cell++, &cv[4] );
            editor.add_cell( cell++, &cv[8] );
            editor.add_cell( cell++, &cv[12] );
            editor.add_cell( cell++, &cv[16] );
            editor.add_cell( cell++, &cv[20] );
          }
        }
      }
    }
    break;
    default:
      error( "StructuredGrid::init : unsupported cell type" );
      break;
  }

  delete[] x;
  delete[] h;
  delete[] i;
  delete[] n;

  //
  editor.close();
}

//-----------------------------------------------------------------------------

} /* namespace dolfin */
