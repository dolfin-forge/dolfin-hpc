// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/IntersectionDetector.h>

#include <dolfin/log/log.h>
#include <dolfin/mesh/entities/Cell.h>
#include <dolfin/mesh/entities/Facet.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/entities/Vertex.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------
IntersectionDetector::IntersectionDetector( Mesh & mesh )
  : gts( mesh )
{
}
//-----------------------------------------------------------------------------
IntersectionDetector::~IntersectionDetector() = default;
//-----------------------------------------------------------------------------
void IntersectionDetector::overlap( std::vector< Point > const & points,
                                    std::vector< size_t > &      cells ) const
{
  // Intersect each segment with mesh
  std::vector< size_t > cc;
  for ( size_t i = 0; i < points.size() - 1; i++ )
  {
    gts.overlap( points[i], points[i + 1], cc );
  }

  // sort cells
  std::sort( cc.begin(), cc.end() );

  // remove repeated cells
  cells.clear();
  cells.push_back( cc[0] );
  size_t k = cc[0];
  for ( size_t i = 1; i < cc.size(); i++ )
  {
    if ( cc[i] > k )
    {
      cells.push_back( cc[i] );
      k = cc[i];
    }
  }
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
