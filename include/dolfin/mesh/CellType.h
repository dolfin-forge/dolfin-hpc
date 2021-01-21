// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_CELL_TYPE_H
#define __DOLFIN_CELL_TYPE_H

#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Point.h>
#include <dolfin/mesh/RefinementPattern.h>

#include <ufc.h>

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
  ~CellType() override;

  /// Clone pattern
  virtual auto clone() const -> CellType* = 0;

  /// Comparison operator
  auto operator==(CellType const& other) const -> bool
  {
    return (this->cell_type == other.cell_type);
  }
  auto operator!=(CellType const& other) const -> bool
  {
    return !(*this == other);
  }

  /// Return type of cell
  inline auto cellType() const -> CellType::Type { return cell_type; }

  /// Return type of cell for facets
  inline auto facetType() const -> CellType::Type { return facet_type; }

  /// Return string of cell type
  virtual auto str() const -> std::string const&;

  /// Return topological dimension of cell
  virtual auto dim() const -> uint = 0;

  /// Return topological dimension of facet
  inline auto facet_dim() const -> uint { return (this->dim() ? this->dim() - 1 : 0); }

  /// Return dimension of Euclidean space
  inline auto space_dim() const -> uint { return std::max(this->dim(), 1ul); }

  /// Return number of entities of given topological dimension
  virtual auto num_entities(uint dim) const -> uint = 0;

  /// Return number of entities of given topological dimensions
  virtual auto num_entities(uint d0, uint d1) const -> uint = 0;

  /// Return number of vertices for entity of given topological dimension
  virtual auto num_vertices(uint dim) const -> uint = 0;

  /// Return orientation of the cell
  virtual auto orientation(Cell const& cell) const -> uint = 0;

  /// Create entities e of given topological dimension from vertices v
  virtual void create_entities(uint** e, uint dim, uint const* v) const = 0;

  /// Order entities locally
  virtual void order_entities(MeshTopology& topology, uint i) const = 0;

  /// Order vertices such that the facet is right-oriented w.r.t. facet normal
  virtual void order_facet(uint vertices[], Facet& facet) const = 0;

  /// Return if mesh connectivities require ordering
  virtual auto connectivity_needs_ordering(uint d0, uint d1) const -> bool = 0;

  /// Initialize mesh connectivities required by ordering
  virtual void initialize_connectivities(Mesh& mesh) const = 0;

  //--- REFINEMENT PATTERN ----------------------------------------------------

  /// Return the cell type to which the pattern applies
  auto pattern_applies(Cell& cell) const -> bool override;

  /// Refine cell uniformly
  void refine_cell(Cell& cell, MeshEditor& editor, uint& current_cell) const override = 0;

  /// Number of cells created by refinement pattern
  auto num_refined_cells() const -> uint override = 0;

  /// Number of vertices created by refinement pattern restricted to each
  /// entity of given topological dimensions
  auto num_refined_vertices(uint dim) const -> uint override = 0;

  //---------------------------------------------------------------------------

  /// Compute (generalized) volume of mesh entity
  virtual auto volume(MeshEntity const& entity) const -> real = 0;

  /// Compute diameter of mesh entity
  virtual auto diameter(MeshEntity const& entity) const -> real = 0;

  /// Compute circumradius of mesh entity
  virtual auto circumradius(MeshEntity const& entity) const -> real = 0;

  /// Compute inradius of mesh entity
  virtual auto inradius(MeshEntity const& entity) const -> real = 0;

  /// Compute coordinates of midpoint
  virtual void midpoint(MeshEntity const& entity, real * p) const = 0;

  /// Compute of given facet with respect to the cell
  virtual void normal(Cell const& cell, uint facet, real * n) const = 0;

  /// Compute the area/length of given facet with respect to the cell
  virtual auto facet_area(Cell const& cell, uint facet) const -> real = 0;

  /// Check if point p intersects the entity
  virtual auto intersects(MeshEntity const& e, Point const& p) const -> bool = 0;

  /// Check if points line connecting p1 and p2 cuts the entity
  virtual auto intersects(MeshEntity const& e, Point const& p1, Point const& p2) const -> bool = 0;

  /// Check if cell c intersects the cell
  virtual auto intersects(MeshEntity& entity, Cell& c) const -> bool;

  //--- REFERENCE CELL --------------------------------------------------------

  /// Create a mesh consisting of the reference cell
  virtual void create_reference_cell(Mesh& mesh) const = 0;

  /// Return coordinates of vertices in the reference cell
  virtual auto reference_vertex(uint i) const -> real const * = 0;

  //---------------------------------------------------------------------------

  /// Return description of cell type
  virtual auto description() const -> std::string = 0;

  /// Display information
  virtual void disp() const = 0;

  /// Common cell type check
  /// ASSERTION: cell vertices in ascending order
  /// ASSERTION: edge vertices in ascending order
  virtual auto check(Cell& cell) const -> bool = 0;

  //---------------------------------------------------------------------------

  /// Create cell type from type (factory function)
  static auto create(CellType::Type type) -> CellType*;

  /// Create cell type from ufc::shape (factory function)
  static auto create(ufc::shape type) -> CellType*;

  /// Create cell type from string (factory function)
  static auto create(std::string const& type) -> CellType*;

  /// Create cell types
  static auto create_all() -> Array<CellType*>;

  /// Create cell types
  static auto create_simplex() -> Array<CellType*>;

  /// Create simplicial cell types
  static auto create_simplex(uint dim) -> CellType*;

  /// Create cell types
  static auto create_hypercube() -> Array<CellType*>;

  /// Create hypercube cell types
  static auto create_hypercube(uint dim) -> CellType*;

protected:

  /// Convert from string to cell type
  static auto type(std::string const& type) -> CellType::Type;

  /// Convert from cell type to string
  static auto str(CellType::Type type) -> std::string;

  std::string const name_;
  CellType::Type const cell_type;
  CellType::Type const facet_type;

  /// Implementation detail after C++11 <algorithm>
  static auto is_sorted_until(uint const * begin, uint const * end) -> uint const *;

  /// Implementation detail after C++11 <algorithm>
  static auto is_sorted(uint const * begin, uint const * end) -> bool;

};

// Helper function

template<class E> inline auto dimension(CellType const& c) -> uint;

//--- TEMPLATE SPECIALIZATION -------------------------------------------------

class Vertex;
template<>
inline auto dimension<Vertex>(CellType const&) -> uint { return 0; }

class Edge;
template<>
inline auto dimension<Edge>(CellType const&) -> uint { return 1; }

class Face;
template<>
inline auto dimension<Face>(CellType const&) -> uint { return 2; }

class Facet;
template<>
inline auto dimension<Facet>(CellType const& c) -> uint { return c.facet_dim(); }

class Cell;
template<>
inline auto dimension<Cell>(CellType const& c) -> uint { return c.dim(); }

//-----------------------------------------------------------------------------

}

#endif /* __DOLFIN_CELL_TYPE_H */
