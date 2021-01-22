// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU GPL Version 2.

#include <dolfin/mesh/AffineMapping.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshGeometry.h>
#include <dolfin/mesh/celltypes/CellType.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Edge.h>
#include <dolfin/mesh/entities/Face.h>
#include <dolfin/mesh/entities/Vertex.h>

#include <algorithm>
#include <cmath>

namespace dolfin
{

#define RM( row, col ) ( ( row ) + ( ( d_ ) * ( col ) ) )

//-----------------------------------------------------------------------------

AffineMapping::AffineMapping( Mesh const & mesh )
  : Mapping()
  , det( 0.0 )
  , J( new real[d_ * d_] )
  , K( new real[d_ * d_] )
  , gdim_( mesh.geometry_dimension() )
{
  dolfin_assert( gdim_ <= d_ );
  std::fill( &J[0], &J[d_ * d_], 0.0 );
  std::fill( &K[0], &K[d_ * d_], 0.0 );
  for ( size_t ip = 0; ip < n_; ++ip )
  {
    std::fill( &p[ip][0], &p[ip][d_], 0.0 );
  }
}

//-----------------------------------------------------------------------------

AffineMapping::~AffineMapping()
{
  delete[] J;
  delete[] K;
}

//-----------------------------------------------------------------------------

void AffineMapping::update( Cell const & cell )
{
  dolfin_assert( n_ >= cell.num_entities( 0 ) );
  switch ( cell.type() )
  {
    case CellType::tetrahedron:
      updateTetrahedron( cell );
      break;
    case CellType::triangle:
      updateTriangle( cell );
      break;
    case CellType::interval:
      updateInterval( cell );
      break;
    default:
      error( "Unknown cell type for AffineMapping." );
      break;
  }
  switch ( gdim_ )
  {
    case 3:
      updateR3();
      break;
    case 2:
      updateR2();
      break;
    case 1:
      updateR1();
      break;
    default:
      error( "Unknown geometrical dimension type for AffineMapping." );
      break;
  }
}

//-----------------------------------------------------------------------------

void AffineMapping::map_from_reference_cell( real const * xref, real * x ) const
{
  //  x = J . xref
  for ( size_t i = 0; i < gdim_; ++i )
  {
    x[i] = p[0][i];
    for ( size_t j = 0; j < gdim_; ++j )
    {
      x[i] += J[RM( i, j )] * xref[j];
    }
  }
}

//-----------------------------------------------------------------------------

void AffineMapping::map_to_reference_cell( real const * x, real * xref ) const
{
  //  xref = K . (x - p0)
  for ( size_t i = 0; i < gdim_; i++ )
  {
    xref[i] = 0.0;
    for ( size_t j = 0; j < gdim_; j++ )
    {
      xref[i] += K[RM( i, j )] * ( x[j] - p[0][j] );
    }
  }
}

//-----------------------------------------------------------------------------

void AffineMapping::updateInterval( Cell const & cell )
{
  dolfin_assert( cell.dim() == 1 );
  std::fill( &J[0], &J[d_ * d_], 0.0 );
  std::fill( &K[0], &K[d_ * d_], 0.0 );
  std::fill( &p[0][0], &p[0][d_], 0.0 );
  std::fill( &p[1][0], &p[1][d_], 0.0 );

  // Get coordinates
  std::vector< size_t > const & vertices = cell.entities( 0 );
  MeshGeometry const &          geom     = cell.mesh().geometry();
  std::copy(
    &geom.x( vertices[0] )[0], &geom.x( vertices[0] )[0] + gdim_, &p[0][0] );
  std::copy(
    &geom.x( vertices[1] )[0], &geom.x( vertices[1] )[0] + gdim_, &p[1][0] );
}

//-----------------------------------------------------------------------------

void AffineMapping::updateTriangle( Cell const & cell )
{
  dolfin_assert( cell.dim() == 2 );
  std::fill( &J[0], &J[d_ * d_], 0.0 );
  std::fill( &K[0], &K[d_ * d_], 0.0 );
  std::fill( &p[0][0], &p[0][d_], 0.0 );
  std::fill( &p[1][0], &p[1][d_], 0.0 );
  std::fill( &p[2][0], &p[2][d_], 0.0 );

  // Get coordinates
  std::vector< size_t > const & vertices = cell.entities( 0 );
  MeshGeometry const &          geom     = cell.mesh().geometry();
  std::copy(
    &geom.x( vertices[0] )[0], &geom.x( vertices[0] )[0] + gdim_, &p[0][0] );
  std::copy(
    &geom.x( vertices[1] )[0], &geom.x( vertices[1] )[0] + gdim_, &p[1][0] );
  std::copy(
    &geom.x( vertices[2] )[0], &geom.x( vertices[2] )[0] + gdim_, &p[2][0] );
}

//-----------------------------------------------------------------------------

void AffineMapping::updateTetrahedron( Cell const & cell )
{
  dolfin_assert( cell.dim() == 3 );
  std::fill( &J[0], &J[d_ * d_], 0.0 );
  std::fill( &K[0], &K[d_ * d_], 0.0 );
  std::fill( &p[0][0], &p[0][d_], 0.0 );
  std::fill( &p[1][0], &p[1][d_], 0.0 );
  std::fill( &p[2][0], &p[2][d_], 0.0 );
  std::fill( &p[3][0], &p[3][d_], 0.0 );

  // Get coordinates
  std::vector< size_t > const & vertices = cell.entities( 0 );
  MeshGeometry const &          geom     = cell.mesh().geometry();
  std::copy(
    &geom.x( vertices[0] )[0], &geom.x( vertices[0] )[0] + gdim_, &p[0][0] );
  std::copy(
    &geom.x( vertices[1] )[0], &geom.x( vertices[1] )[0] + gdim_, &p[1][0] );
  std::copy(
    &geom.x( vertices[2] )[0], &geom.x( vertices[2] )[0] + gdim_, &p[2][0] );
  std::copy(
    &geom.x( vertices[3] )[0], &geom.x( vertices[3] )[0] + gdim_, &p[3][0] );
}

//-----------------------------------------------------------------------------

void AffineMapping::updateR1()
{
  // Compute Jacobian of map
  J[RM( 0, 0 )] = p[1][0] - p[0][0];

  // Compute determinant
  det = J[RM( 0, 0 )];

  // Check determinant
  if ( std::fabs( det ) < DOLFIN_EPS )
  {
    error( "Map from reference cell is singular." );
  }

  // Compute inverse of Jacobian
  K[RM( 0, 0 )] = +1.0 / det;

  // Take absolute value of determinant
  det = std::fabs( det );
}

//-----------------------------------------------------------------------------

void AffineMapping::updateR2()
{
  // Compute Jacobian of map
  J[RM( 0, 0 )] = p[1][0] - p[0][0];
  J[RM( 0, 1 )] = p[2][0] - p[0][0];
  J[RM( 1, 0 )] = p[1][1] - p[0][1];
  J[RM( 1, 1 )] = p[2][1] - p[0][1];

  // Compute determinant
  det = J[RM( 0, 0 )] * J[RM( 1, 1 )] - J[RM( 0, 1 )] * J[RM( 1, 0 )];

  // Check determinant
  if ( std::fabs( det ) < DOLFIN_EPS )
  {
    error( "Map from reference cell is singular." );
  }

  // Compute inverse of Jacobian
  K[RM( 0, 0 )] = +J[RM( 1, 1 )] / det;
  K[RM( 0, 1 )] = -J[RM( 0, 1 )] / det;
  K[RM( 1, 0 )] = -J[RM( 1, 0 )] / det;
  K[RM( 1, 1 )] = +J[RM( 0, 0 )] / det;

  // Take absolute value of determinant
  det = std::fabs( det );
}

//-----------------------------------------------------------------------------

void AffineMapping::updateR3()
{
  // Compute Jacobian of map
  J[RM( 0, 0 )] = p[1][0] - p[0][0];
  J[RM( 0, 1 )] = p[2][0] - p[0][0];
  J[RM( 0, 2 )] = p[3][0] - p[0][0];
  J[RM( 1, 0 )] = p[1][1] - p[0][1];
  J[RM( 1, 1 )] = p[2][1] - p[0][1];
  J[RM( 1, 2 )] = p[3][1] - p[0][1];
  J[RM( 2, 0 )] = p[1][2] - p[0][2];
  J[RM( 2, 1 )] = p[2][2] - p[0][2];
  J[RM( 2, 2 )] = p[3][2] - p[0][2];

  // Compute sub-determinants
  real d00 = J[RM( 1, 1 )] * J[RM( 2, 2 )] - J[RM( 1, 2 )] * J[RM( 2, 1 )];
  real d01 = J[RM( 1, 2 )] * J[RM( 2, 0 )] - J[RM( 1, 0 )] * J[RM( 2, 2 )];
  real d02 = J[RM( 1, 0 )] * J[RM( 2, 1 )] - J[RM( 1, 1 )] * J[RM( 2, 0 )];
  real d10 = J[RM( 0, 2 )] * J[RM( 2, 1 )] - J[RM( 0, 1 )] * J[RM( 2, 2 )];
  real d11 = J[RM( 0, 0 )] * J[RM( 2, 2 )] - J[RM( 0, 2 )] * J[RM( 2, 0 )];
  real d12 = J[RM( 0, 1 )] * J[RM( 2, 0 )] - J[RM( 0, 0 )] * J[RM( 2, 1 )];
  real d20 = J[RM( 0, 1 )] * J[RM( 1, 2 )] - J[RM( 0, 2 )] * J[RM( 1, 1 )];
  real d21 = J[RM( 0, 2 )] * J[RM( 1, 0 )] - J[RM( 0, 0 )] * J[RM( 1, 2 )];
  real d22 = J[RM( 0, 0 )] * J[RM( 1, 1 )] - J[RM( 0, 1 )] * J[RM( 1, 0 )];

  // Compute determinant
  det = J[RM( 0, 0 )] * d00 + J[RM( 1, 0 )] * d10 + J[RM( 2, 0 )] * d20;

  // Check determinant
  if ( std::fabs( det ) < DOLFIN_EPS )
  {
    error( "Map from reference cell is singular." );
  }

  // Compute inverse of Jacobian
  K[RM( 0, 0 )] = d00 / det;
  K[RM( 0, 1 )] = d10 / det;
  K[RM( 0, 2 )] = d20 / det;
  K[RM( 1, 0 )] = d01 / det;
  K[RM( 1, 1 )] = d11 / det;
  K[RM( 1, 2 )] = d21 / det;
  K[RM( 2, 0 )] = d02 / det;
  K[RM( 2, 1 )] = d12 / det;
  K[RM( 2, 2 )] = d22 / det;

  // Take absolute value of determinant
  det = std::fabs( det );
}

} // namespace dolfin
