// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-06-21
// Last changed: 2008-02-18

#ifndef __DOLFIN_INTERSECTION_DETECTOR_H
#define __DOLFIN_INTERSECTION_DETECTOR_H

#include <dolfin/mesh/GTSInterface.h>

namespace dolfin
{

class Mesh;
class Cell;
class Point;

template<class T> class Array;

class IntersectionDetector
{

public:

  /// Constructor
  IntersectionDetector(Mesh& mesh);

  /// Destructor
  ~IntersectionDetector();

  /// Compute overlap with mesh
  void overlap(Cell& c, Array<uint>& overlap) const;

  /// Compute overlap with point
  void overlap(Point const& p, Array<uint>& overlap) const;

  /// Compute overlap with bounding box
  void overlap(Point const& p1, Point const& p2, Array<uint>& overlap) const;

  /// Compute which cells are intersected by a polygon (defined by points)
  void overlap(Array<Point> const& points, Array<uint>& overlap) const;

private:

  GTSInterface gts;
};

} /* namespace dolfin */

#endif /* __DOLFIN_INTERSECTION_DETECTOR_H */
