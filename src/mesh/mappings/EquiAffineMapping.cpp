// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU GPL Version 2.

#include <dolfin/mesh/mappings/EquiAffineMapping.h>

#include <dolfin/common/constants.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/celltypes/CellType.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Edge.h>
#include <dolfin/mesh/entities/Face.h>
#include <dolfin/mesh/entities/Vertex.h>

#include <algorithm>

namespace dolfin
{

#define RM( row, col, nrow ) ( ( row ) + ( ( nrow ) * ( col ) ) )

//-----------------------------------------------------------------------------
EquiAffineMapping::EquiAffineMapping( Mesh const & mesh )
  : Mapping()
  , J( nullptr )
  , K( nullptr )
  , gdim_( mesh.geometry_dimension() )
{
  dolfin_assert( gdim_ <= d_ );
  det = 0.0;
  J   = new real[d_ * d_];
  std::fill( &J[0], &J[d_ * d_], 0.0 );
  K = new real[d_ * d_];
  std::fill( &K[0], &K[d_ * d_], 0.0 );

  std::fill( &p0[0], &p0[d_], 0.0 );
  std::fill( &p1[0], &p1[d_], 0.0 );
  std::fill( &p2[0], &p2[d_], 0.0 );
  std::fill( &p3[0], &p3[d_], 0.0 );
}
//-----------------------------------------------------------------------------
EquiAffineMapping::~EquiAffineMapping()
{
  delete[] J;
  delete[] K;
}
//-----------------------------------------------------------------------------
void EquiAffineMapping::update( Cell const & cell )
{
  switch ( cell.type() )
  {
    case CellType::triangle:
      updateTriangle( cell );
      break;
    case CellType::tetrahedron:
      updateTetrahedron( cell );
      break;
    default:
      error( "Unknown cell type for EquiAffine map." );
      break;
  }
}
//-----------------------------------------------------------------------------
void EquiAffineMapping::map_from_reference_cell( real const * xref,
                                                 real *       x ) const
{
  //  J.mult(P, p);
  for ( size_t i = 0; i < gdim_; ++i )
  {
    x[i] = p0[i];
    for ( size_t j = 0; j < gdim_; ++j )
    {
      x[i] += J[RM( i, j, d_ )] * xref[j];
    }
  }
}
//-----------------------------------------------------------------------------
void EquiAffineMapping::map_to_reference_cell( real const * x,
                                               real *       xref ) const
{
  //  K.mult(p, P);
  for ( size_t i = 0; i < gdim_; i++ )
  {
    xref[i] = 0.0;
    for ( size_t j = 0; j < gdim_; j++ )
    {
      xref[i] += K[RM( i, j, d_ )] * ( x[j] - p0[j] );
    }
  }
}
//-----------------------------------------------------------------------------
void EquiAffineMapping::updateTriangle( Cell const & cell )
{

  dolfin_assert( cell.dim() == 2 );

  // Get coordinates
  Array< size_t > const & vertices = cell.entities( 0 );
  MeshGeometry const &    geom     = cell.mesh().geometry();
  std::fill( &p0[0], &p0[d_], 0.0 );
  std::fill( &p1[0], &p1[d_], 0.0 );
  std::fill( &p2[0], &p2[d_], 0.0 );
  std::copy(
    &geom.x( vertices[0] )[0], &geom.x( vertices[0] )[0] + gdim_, &p0[0] );
  std::copy(
    &geom.x( vertices[1] )[0], &geom.x( vertices[1] )[0] + gdim_, &p1[0] );
  std::copy(
    &geom.x( vertices[2] )[0], &geom.x( vertices[2] )[0] + gdim_, &p2[0] );

  // Compute Jacobian of map
  real a = 1.0 / sqrt( 3.0 );

  J[RM( 0, 0, d_ )] = -a * p0[0] + 2 * a * p1[0] - a * p2[0];
  J[RM( 1, 0, d_ )] = -a * p0[1] + 2 * a * p1[1] - a * p2[1];
  J[RM( 2, 0, d_ )] = +0.0;
  J[RM( 0, 1, d_ )] = -p0[0] + p2[0];
  J[RM( 1, 1, d_ )] = -p0[1] + p2[1];
  J[RM( 2, 1, d_ )] = +0.0;
  J[RM( 0, 2, d_ )] = +0.0;
  J[RM( 1, 2, d_ )] = +0.0;
  J[RM( 2, 2, d_ )] = +0.0;

  // Compute determinant
  det = J[RM( 0, 0, d_ )] * J[RM( 1, 1, d_ )]
        - J[RM( 0, 1, d_ )] * J[RM( 1, 0, d_ )];

  // Check determinant
  if ( fabs( det ) < DOLFIN_EPS )
  {
    error( "Map from reference cell is singular." );
  }

  // Compute inverse of Jacobian
  K[RM( 0, 0, d_ )] = +J[RM( 1, 1, d_ )] / det;
  K[RM( 0, 1, d_ )] = -J[RM( 0, 1, d_ )] / det;
  K[RM( 0, 2, d_ )] = +0.0;
  K[RM( 1, 0, d_ )] = -J[RM( 1, 0, d_ )] / det;
  K[RM( 1, 1, d_ )] = +J[RM( 0, 0, d_ )] / det;
  K[RM( 1, 2, d_ )] = +0.0;
  K[RM( 2, 0, d_ )] = +0.0;
  K[RM( 2, 1, d_ )] = +0.0;
  K[RM( 2, 2, d_ )] = +0.0;

  // Take absolute value of determinant
  det = fabs( det );
}
//-----------------------------------------------------------------------------
void EquiAffineMapping::updateTetrahedron( Cell const & cell )
{
  dolfin_assert( cell.dim() == 3 );

  // Get coordinates
  Array< size_t > const & vertices = cell.entities( 0 );
  MeshGeometry const &    geom     = cell.mesh().geometry();
  std::fill( &p0[0], &p0[d_], 0.0 );
  std::fill( &p1[0], &p1[d_], 0.0 );
  std::fill( &p2[0], &p2[d_], 0.0 );
  std::fill( &p3[0], &p3[d_], 0.0 );
  std::copy(
    &geom.x( vertices[0] )[0], &geom.x( vertices[0] )[0] + gdim_, &p0[0] );
  std::copy(
    &geom.x( vertices[1] )[0], &geom.x( vertices[1] )[0] + gdim_, &p1[0] );
  std::copy(
    &geom.x( vertices[2] )[0], &geom.x( vertices[2] )[0] + gdim_, &p2[0] );
  std::copy(
    &geom.x( vertices[3] )[0], &geom.x( vertices[3] )[0] + gdim_, &p3[0] );

  // Compute Jacobian of map
  real a = 1.0 / sqrt( 3.0 );
  real b = 1.0 / sqrt( 2.0 ) * a;
  real c = sqrt( 3.0 ) / sqrt( 2.0 );

  J[RM( 0, 0, d_ )] = -a * p0[0] + 2 * a * p1[0] - a * p2[0];
  J[RM( 0, 1, d_ )] = -p0[0] + p2[0];
  J[RM( 0, 2, d_ )] = -b * p0[0] - b * p1[0] - b * p2[0] + c * p3[0];
  J[RM( 1, 0, d_ )] = -a * p0[1] + 2 * a * p1[1] - a * p2[1];
  J[RM( 1, 1, d_ )] = -p0[1] + p2[1];
  J[RM( 1, 2, d_ )] = -b * p0[1] - b * p1[1] - b * p2[1] + c * p3[1];
  J[RM( 2, 0, d_ )] = -a * p0[2] + 2 * a * p1[2] - a * p2[2];
  J[RM( 2, 1, d_ )] = -p0[2] + p2[2];
  J[RM( 2, 2, d_ )] = -b * p0[2] - b * p1[2] - b * p2[2] + c * p3[2];

  // Compute sub-determinants
  real d00 = J[RM( 1, 1, d_ )] * J[RM( 2, 2, d_ )]
             - J[RM( 1, 2, d_ )] * J[RM( 2, 1, d_ )];
  real d01 = J[RM( 1, 2, d_ )] * J[RM( 2, 0, d_ )]
             - J[RM( 1, 0, d_ )] * J[RM( 2, 2, d_ )];
  real d02 = J[RM( 1, 0, d_ )] * J[RM( 2, 1, d_ )]
             - J[RM( 1, 1, d_ )] * J[RM( 2, 0, d_ )];

  real d10 = J[RM( 0, 2, d_ )] * J[RM( 2, 1, d_ )]
             - J[RM( 0, 1, d_ )] * J[RM( 2, 2, d_ )];
  real d11 = J[RM( 0, 0, d_ )] * J[RM( 2, 2, d_ )]
             - J[RM( 0, 2, d_ )] * J[RM( 2, 0, d_ )];
  real d12 = J[RM( 0, 1, d_ )] * J[RM( 2, 0, d_ )]
             - J[RM( 0, 0, d_ )] * J[RM( 2, 1, d_ )];

  real d20 = J[RM( 0, 1, d_ )] * J[RM( 1, 2, d_ )]
             - J[RM( 0, 2, d_ )] * J[RM( 1, 1, d_ )];
  real d21 = J[RM( 0, 2, d_ )] * J[RM( 1, 0, d_ )]
             - J[RM( 0, 0, d_ )] * J[RM( 1, 2, d_ )];
  real d22 = J[RM( 0, 0, d_ )] * J[RM( 1, 1, d_ )]
             - J[RM( 0, 1, d_ )] * J[RM( 1, 0, d_ )];

  // Compute determinant
  det =
    J[RM( 0, 0, d_ )] * d00 + J[RM( 1, 0, d_ )] * d10 + J[RM( 2, 0, d_ )] * d20;

  // Check determinant
  if ( fabs( det ) < DOLFIN_EPS )
  {
    error( "Map from reference cell is singular." );
  }

  // Compute inverse of Jacobian
  K[RM( 0, 0, d_ )] = d00 / det;
  K[RM( 0, 1, d_ )] = d10 / det;
  K[RM( 0, 2, d_ )] = d20 / det;
  K[RM( 1, 0, d_ )] = d01 / det;
  K[RM( 1, 1, d_ )] = d11 / det;
  K[RM( 1, 2, d_ )] = d21 / det;
  K[RM( 2, 0, d_ )] = d02 / det;
  K[RM( 2, 1, d_ )] = d12 / det;
  K[RM( 2, 2, d_ )] = d22 / det;

  // Take absolute value of determinant
  det = fabs( det );
}
//-----------------------------------------------------------------------------

}
