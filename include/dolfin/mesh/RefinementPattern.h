// Copyright (C) 2015 Aurelien Larcher
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_REFINEMENT_PATTERN_H
#define __DOLFIN_REFINEMENT_PATTERN_H

#include <dolfin/common/types.h>

namespace dolfin
{

class Cell;
class Mesh;
class MeshEditor;

/**
 *  @class  RefinementPattern
 *
 *  @brief  This class provides an interface for cell refinement patterns.
 *
 */

class RefinementPattern
{

public:

  ///
  RefinementPattern() = default;

  ///--- INTERFACE ------------------------------------------------------------

  /// Return if the cell to which the pattern applies is valid
  virtual bool pattern_applies(Cell& cell) const = 0;

  /// Refine cell uniformly
  virtual void refine_cell(Cell& cell, MeshEditor& editor, uint& current_cell) const = 0;

  /// Number of cells created by refinement pattern
  virtual uint num_refined_cells() const = 0;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  virtual uint num_refined_vertices(uint dim) const = 0;

  ///--------------------------------------------------------------------------

  /// Return the number of vertices after refinement
  uint num_refined_vertices(Mesh const& mesh) const;

  /// Return the number of cells after refinement
  uint num_refined_cells(Mesh const& mesh) const;

protected:

  ///
  virtual ~RefinementPattern() = default;

};

} /* namespace dolfin */

#endif /*  __DOLFIN_REFINEMENT_PATTERN_H */
