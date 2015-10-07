// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __DOLFIN_UFL_CELL_H
#define __DOLFIN_UFL_CELL_H

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

#include <dolfin/common/types.h>

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
  explicit Cell(Domain::Type const& type);

  ///
  explicit Cell(Domain const& domain);

  ///
  explicit Cell(Domain const& domain, Space const& space);

  ///
  explicit Cell(repr_t const& repr);

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

  /// UFL: Return whether this cell is undefined
  bool is_undefined() const;

  /// UFL: Return the domain of the cell
  Domain const& domain() const;

  /// UFL: Return the domain of the facet of this cell
  Domain const& facet_domain() const;

  /// UFL: Return the number of facets this cell has
  dolfin::uint num_facets() const;

  /// UFL: Return the dimension of the space this cell is embedded in
  dolfin::uint geometric_dimension() const;

  /// UFL: Return the dimension of the topology of this cell
  dolfin::uint topological_dimension() const;

  /// The dimension of the cell is only valid is the geometric and topological
  /// dimensions are the same which does not seem to be useful.
  ///dolfin::uint const d() const;

  //--- INTERFACE inherited from UFLClass -------------------------------------

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

  ///
  void display() const;

private:

  Domain const domain_;
  Space const space_;
  Domain const facet_domain_;

  // Do not put the initialization of these variables at the end as the
  // GeometricalQuantities require the representation string to be valid.
  repr_t const repr_;
  std::string const str_;

  bool const invalid_;

  dolfin::uint const geometric_dimension_;
  dolfin::uint const topological_dimension_;

  /// GeometricalQuantities
  CellSurfaceArea const cell_surface_area_;
  CellVolume const cell_volume_;
  Circumradius const circumradius_;
  FacetArea const facet_area_;
  FacetNormal const facet_normal_;
  SpatialCoordinate const x_;
};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_CELL_H */
