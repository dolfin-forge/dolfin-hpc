// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_GTS_INTERFACE_H
#define __DOLFIN_GTS_INTERFACE_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/mesh/Point.h>

// Forward declarations
struct _GtsBBox;
using GtsBBox = _GtsBBox;
struct _GNode;
using GNode = _GNode;

namespace dolfin
{

class Mesh;
class Cell;

class GTSInterface
{
public:

  ///
  GTSInterface(Mesh& mesh);

  ///
  ~GTSInterface();

  /// Compute cells overlapping c
  void overlap(Cell& c, Array<uint>& cells) const;

  /// Compute cells overlapping p
  void overlap(Point const& p, Array<uint>& cells) const;

  /// Compute cells overlapping bounding box constructed from p1 and p2
  void overlap(Point const& p1, Point const& p2, Array<uint>& cells) const;

private:

  /// Construct bounding box of cell
  GtsBBox* bboxCell(Cell&) const;

  /// Construct bounding box of a single point
  GtsBBox* bboxPoint(Point const& p) const;

  /// Construct bounding box of a pair of points
  GtsBBox* bboxPoint(Point const& p1, Point const& p2) const;

  /// Construct hierarchical space partition tree of mesh
  void buildCellTree(void);

  Mesh& mesh_;
  mutable GNode * tree_;
};

}
#endif
