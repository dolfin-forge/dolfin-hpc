// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-06-05
// Last changed: 2008-06-20

#ifndef __DOLFIN_CELL_TYPE_H
#define __DOLFIN_CELL_TYPE_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/RefinementPattern.h>
#include <dolfin/ufl/UFLCell.h>

#include <string>

namespace dolfin
{

class Cell;
class Facet;
class Mesh;
class MeshEditor;
class MeshEntity;
class MeshTopology;

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
  CellType(std::string const& name, CellType::Type cell_type, CellType::Type facet_type);

  /// Destructor
  virtual ~CellType();

  /// Clone pattern
  virtual CellType* clone() const = 0;

  /// Return type of cell
  inline CellType::Type cellType() const { return cell_type; }

  /// Return type of cell for facets
  inline CellType::Type facetType() const { return facet_type; }

  /// Return string of cell type
  virtual std::string const& str() const;

  /// Return topological dimension of cell
  virtual uint dim() const = 0;

  /// Return topological dimension of facet
  inline uint facet_dim() const { return (this->dim() ? this->dim() - 1 : 0); }

  /// Return dimension of Euclidean space
  inline uint space_dim() const { return std::max(this->dim(), 1u); }

  /// Return number of entitites of given topological dimension
  virtual uint num_entities(uint dim) const = 0;

  /// Return number of vertices for entity of given topological dimension
  virtual uint num_vertices(uint dim) const = 0;

  /// Return orientation of the cell
  virtual uint orientation(Cell const& cell) const = 0;

  /// Create entities e of given topological dimension from vertices v
  virtual void create_entities(uint** e, uint dim, uint const* v) const = 0;

  /// Order entities locally
  virtual void order_entities(MeshTopology& topology, uint i) const = 0;

  /// Order vertices such that the facet is right-oriented w.r.t. facet normal
  virtual void order_facet(uint vertices[], Facet& facet) const = 0;

  /// Return if mesh connectivities require ordering
  virtual bool connectivity_needs_ordering(uint d0, uint d1) const = 0;

  /// Initialize mesh connectivities required by ordering
  virtual void initialize_connectivities(Mesh& mesh) const = 0;

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

  /// Compute inradius of mesh entity
  virtual real inradius(MeshEntity const& entity) const = 0;

  /// Compute coordinates of midpoint
  virtual void midpoint(MeshEntity const& entity, real * p) const = 0;

  /// Compute of given facet with respect to the cell
  virtual void normal(Cell const& cell, uint facet, real * n) const = 0;

  /// Compute the area/length of given facet with respect to the cell
  virtual real facet_area(Cell const& cell, uint facet) const = 0;

  /// Check if point p intersects the entity
  virtual bool intersects(MeshEntity const& e, Point const& p) const = 0;

  /// Check if points line connecting p1 and p2 cuts the entity
  virtual bool intersects(MeshEntity const& e, Point const& p1, Point const& p2) const = 0;

  /// Check if cell c intersects the cell
  virtual bool intersects(MeshEntity& entity, Cell& c) const;

  //--- REFERENCE CELL --------------------------------------------------------

  /// Create a mesh consisting of the reference cell
  virtual void create_reference_cell(Mesh& mesh) const = 0;

  /// Return coordinates of vertices in the reference cell
  virtual real const * reference_vertex(uint i) const = 0;

  //---------------------------------------------------------------------------

  /// Return description of cell type
  virtual std::string description() const = 0;

  /// Display information
  virtual void disp() const = 0;

  /// Common cell type check
  /// ASSERTION: cell vertices in ascending order
  /// ASSERTION: edge vertices in ascending order
  virtual bool check(Cell& cell) const = 0;

  /// UFL binding
  operator ufl::Cell const&() const { return ufl_; }

  //---------------------------------------------------------------------------

  /// Create cell type from type (factory function)
  static CellType* create(CellType::Type type);

  /// Create cell type from string (factory function)
  static CellType* create(std::string const& type);

  /// Create cell type from UFL type (factory function)
  static CellType* create(ufl::Cell const& cell);

  /// Create cell types
  static Array<CellType*> create_all();

  /// Create cell types
  static Array<CellType*> create_simplex();

  /// Create simplicial cell types
  static CellType* create_simplex(uint dim);

  /// Create cell types
  static Array<CellType*> create_hypercube();

  /// Create hypercube cell types
  static CellType* create_hypercube(uint dim);

protected:

  /// Convert from string to cell type
  static CellType::Type type(std::string const& type);

  /// Convert from cell type to string
  static std::string str(CellType::Type type);

  /// Convert from cell type to UFL cell type
  static ufl::Domain::Type ufldomain(CellType::Type type);

  std::string const name_;
  CellType::Type const cell_type;
  CellType::Type const facet_type;

  /// Implementation detail after C++11 <algorithm>
  static uint const * is_sorted_until(uint const * begin, uint const * end);

  /// Implementation detail after C++11 <algorithm>
  static bool is_sorted(uint const * begin, uint const * end);

private:

  ufl::Cell ufl_;

};

//-----------------------------------------------------------------------------

}

#endif /* __DOLFIN_CELL_TYPE_H */
