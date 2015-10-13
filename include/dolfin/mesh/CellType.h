// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Aurelien Larcher, 2015.
//
// First added:  2006-06-05
// Last changed: 2008-06-20

#ifndef __DOLFIN_CELL_TYPE_H
#define __DOLFIN_CELL_TYPE_H

#include <dolfin/common/types.h>

#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/RefinementPattern.h>
#include <dolfin/ufl/UFLCell.h>

#include <string>

namespace dolfin
{

class Cell;
class Mesh;
class MeshEditor;
class MeshEntity;
class Point;

/**
 *  @class  CellType
 *
 *  @brief  This class provides a common interface for different cell types.
 *          Each cell type implements mesh functionality that is specific to
 *          a certain type of cell.
 *          CellType implements the RefinementPattern interface which provides
 *          the default refinement pattern for the cell.
 */

class CellType : public RefinementPattern
{
public:

  /// Enum for different cell types
  /// NOTE: Type index is set explicitly to ensure consistency and stable. but
  ///       the value of the integral type itself should not be relied on.
  enum Type { point         = 0,
              interval      = 1,
              triangle      = 2,
              tetrahedron   = 3,
              quadrilateral = 4,
              hexahedron    = 6 };

  /// Constructor
  CellType(CellType::Type cell_type, CellType::Type facet_type);

  /// Destructor
  virtual ~CellType();

  /// Clone pattern
  virtual CellType* clone() const = 0;

  /// Create cell type from type (factory function)
  static CellType* create(CellType::Type type);

  /// Create cell type from string (factory function)
  static CellType* create(std::string type);

  /// Convert from string to cell type
  static CellType::Type string2type(std::string type);

  /// Convert from cell type to string
  static std::string type2string(CellType::Type type);

  /// Convert from cell type to UFL cell type
  static ufl::Domain::Type type2ufldomain(CellType::Type type);

  /// Return type of cell
  inline CellType::Type cellType() const { return cell_type; }

  /// Return type of cell for facets
  inline CellType::Type facetType() const { return facet_type; }

  /// Return topological dimension of cell
  virtual uint dim() const = 0;

  /// Return number of entitites of given topological dimension
  virtual uint numEntities(uint dim) const = 0;

  /// Return number of vertices for entity of given topological dimension
  virtual uint numVertices(uint dim) const = 0;

  /// Return orientation of the cell
  virtual uint orientation(Cell const& cell) const = 0;

  /// Create entities e of given topological dimension from vertices v
  virtual void createEntities(uint** e, uint dim, uint const* v) const = 0;

  /// Order entities locally
  virtual void orderEntities(Cell& cell) const = 0;

  //--- REFINEMENT PATTERN ----------------------------------------------------

  /// Return the cell type to which the pattern applies
  bool pattern_applies(Cell& cell) const;

  /// Refine cell uniformly
  virtual void refine_cell(Cell& cell, MeshEditor& editor, uint& current_cell) const = 0;

  /// Number of cells created by refinement pattern
  virtual uint num_refined_cells() const = 0;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  virtual uint num_refined_vertices(uint dim) const = 0;

  /// Return if refinement pattern requires entities of given dimension
  virtual bool refinement_needs_entities(uint dim) const = 0;

  //---------------------------------------------------------------------------

  /// Compute (generalized) volume of mesh entity
  virtual real volume(MeshEntity const& entity) const = 0;

  /// Compute diameter of mesh entity
  virtual real diameter(MeshEntity const& entity) const = 0;

  /// Compute circumradius of mesh entity
  virtual real circumradius(MeshEntity const& entity) const = 0;

  /// Compute coordinates of midpoint
  virtual Point midpoint(MeshEntity const& entity) const = 0;

  /// Compute of given facet with respect to the cell
  virtual Point normal(Cell const& cell, uint facet) const = 0;

  /// Compute the area/length of given facet with respect to the cell
  virtual real facetArea(Cell const& cell, uint facet) const = 0;

  /// Check if point p intersects the entity
  virtual bool intersects(MeshEntity const& e, Point const& p) const = 0;

  /// Check if points line connecting p1 and p2 cuts the entity
  virtual bool intersects(MeshEntity const& e, Point const& p1, Point const& p2) const = 0;

  /// Check if cell c intersects the cell
  virtual bool intersects(MeshEntity& entity, Cell& c) const;

  /// Create a mesh consisting of the reference cell
  virtual Mesh create_reference_cell() const = 0;

  /// Return description of cell type
  virtual std::string description() const = 0;

  /// Display information
  virtual void disp() const = 0;

  /// Common cell type check
  /// ASSERTION: cell vertices in ascending order
  /// ASSERTION: edge vertices in ascending order
  virtual void check(Cell& cell) const = 0;

  /// UFL binding
  operator ufl::Cell const&() const { return ufl_; }

protected:

  CellType::Type cell_type;
  CellType::Type facet_type;

private:

  /// Implementation detail after C++11 <algorithm>
  static uint const * is_sorted_until(uint const * begin, uint const * end);

  /// Implementation detail after C++11 <algorithm>
  static bool is_sorted(uint const * begin, uint const * end);

  ufl::Cell ufl_;

};

//-----------------------------------------------------------------------------

}

#endif /* __DOLFIN_CELL_TYPE_H */
