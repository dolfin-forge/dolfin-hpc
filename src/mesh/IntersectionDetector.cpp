// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Johan Jansson 2006.
// Modified by Ola Skavhaug 2006.
// Modified by Dag Lindbo 2008.
//
// First added:  2006-06-21
// Last changed: 2008-02-18

#include <dolfin/mesh/IntersectionDetector.h>

#include <dolfin/common/Array.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>

#include <algorithm>

namespace dolfin
{

//-----------------------------------------------------------------------------
IntersectionDetector::IntersectionDetector(Mesh& mesh) :
    gts(mesh)
{
}
//-----------------------------------------------------------------------------
IntersectionDetector::~IntersectionDetector()
{
}
//-----------------------------------------------------------------------------
void IntersectionDetector::overlap(Array<Point> const& points,
                                   Array<uint>& cells) const
{
  // Intersect each segment with mesh
  Array<uint> cc;
  for (uint i = 0; i < points.size() - 1; i++)
  {
    gts.overlap(points[i], points[i + 1], cc);
  }

  // sort cells
  std::sort(cc.begin(), cc.end());

  // remove repeated cells
  cells.clear();
  cells.push_back(cc[0]);
  uint k = cc[0];
  for (uint i = 1; i < cc.size(); i++)
  {
    if (cc[i] > k)
    {
      cells.push_back(cc[i]);
      k = cc[i];
    }
  }
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */

