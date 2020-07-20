// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_INTERSECTION_DETECTOR_H
#define __DOLFIN_INTERSECTION_DETECTOR_H

#include <dolfin/mesh/GTSInterface.h>

namespace dolfin
{

class Mesh;
class Cell;

class IntersectionDetector
{

public:

  /// Constructor
  IntersectionDetector(Mesh& mesh);

  /// Destructor
  ~IntersectionDetector();

  /// Compute overlap with mesh
  inline void overlap(Cell& c, Array<uint>& cells) const
  {
    gts.overlap(c, cells);
  }

  /// Compute overlap with point
  inline void overlap(Point const& p, Array<uint>& cells) const
  {
    gts.overlap(p, cells);
  }

  /// Compute overlap with bounding box
  inline void overlap(Point const& p1, Point const& p2, Array<uint>& cells) const
  {
    gts.overlap(p1, p2, cells);
  }

  /// Compute which cells are intersected by a polygon (defined by points)
  void overlap(Array<Point> const& points, Array<uint>& overlap) const;

private:

  GTSInterface gts;
};

} /* namespace dolfin */

#endif /* __DOLFIN_INTERSECTION_DETECTOR_H */
