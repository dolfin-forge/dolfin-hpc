// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_CELL_H_
#define __UFL_CELL_H_

#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLDomain.h>
#include <dolfin/ufl/UFLSpace.h>

// Geometrical quantities
#include <dolfin/ufl/UFLCellSurfaceArea.h>
#include <dolfin/ufl/UFLCellVolume.h>
#include <dolfin/ufl/UFLCircumradius.h>
#include <dolfin/ufl/UFLFacetArea.h>
#include <dolfin/ufl/UFLFacetNormal.h>
#include <dolfin/ufl/UFLSpatialCoordinate.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  Cell
 *
 *  @brief  Provides an interface complying with UFL Cell.
 */

class Cell : public Class
{

public:

  ///
  Cell(Domain::Type const& domain);

  ///
  Cell(Domain::Type const& domain, Space const& space);

  ///
  ~Cell();

  //--- INTERFACE -------------------------------------------------------------

  /// UFL geometry value: The global spatial coordinates
  SpatialCoordinate const& x() const;

  /// UFL geometry value: The facet normal on the cell boundary
  FacetNormal const& n() const;

  /// UFL geometry value: The volume of the cell
  CellVolume const& volume() const;

  /// UFL geometry value: The circumradius of the cell
  Circumradius const& circumradius() const;

  /// UFL geometry value: The area of a facet of the cell
  FacetArea const& facet_area() const;

  /// UFL geometry value: The total surface area of the cell
  CellSurfaceArea const& surface_area() const;

  //---------------------------------------------------------------------------

  /// Return whether this cell is undefined
  bool const& is_undefined() const;

  /// Return the domain of the cell
  Domain::Type const& domain() const;

  /// Return the domain of the facet of this cell
  Domain::Type const& facet_domain() const;

  /// Return the number of facets this cell has
  uint const num_facets() const;

  /// Return the dimension of the space this cell is embedded in
  uint const& geometric_dimension() const;

  /// Return the dimension of the topology of this cell
  uint const& topological_dimension() const;

  /// The dimension of the cell is only valid is the geometric and topological
  /// dimensions are the same which does not seem to be useful.
  ///uint const d() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const repr() const;

  /// __str__
  std::string const str() const;

private:

  Domain::Type const domain_;
  Space const space_;

  bool const invalid_;

  uint const geometric_dimension_;
  uint const topological_dimension_;

  /// Geometrical quantities
  CellSurfaceArea const cell_surface_area_;
  CellVolume const cell_volume_;
  Circumradius const circumradius_;
  FacetArea const facet_area_;
  FacetNormal const facet_normal_;
  SpatialCoordinate const x_;

  repr_t const repr_;
  std::string const str_;

};

} /* namespace ufl */
#endif /* __UFL_CELL_H */
