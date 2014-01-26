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

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLCell
 *
 *  @brief  Provides an interface complying with UFL Cell.
 */

class UFLCell : public UFLClass
{

public:

  ///
  UFLCell(UFLDomain::Type const& domain);

  ///
  UFLCell(UFLDomain::Type const& domain, UFLSpace const& space);

  ///
  ~UFLCell();

  //--- INTERFACE -------------------------------------------------------------

  /// UFL geometry value: The global spatial coordinates
  UFLSpatialCoordinate const& x() const;

  /// UFL geometry value: The facet normal on the cell boundary
  UFLFacetNormal const& n() const;

  /// UFL geometry value: The volume of the cell
  UFLCellVolume const& volume() const;

  /// UFL geometry value: The circumradius of the cell
  UFLCircumradius const& circumradius() const;

  /// UFL geometry value: The area of a facet of the cell
  UFLFacetArea const& facet_area() const;

  /// UFL geometry value: The total surface area of the cell
  UFLCellSurfaceArea const& surface_area() const;

  //---------------------------------------------------------------------------

  /// Return whether this cell is undefined
  bool const& is_undefined() const;

  /// Return the domain of the cell
  UFLDomain::Type const& domain() const;

  /// Return the domain of the facet of this cell
  UFLDomain::Type const& facet_domain() const;

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
  std::string const repr() const;

  /// __str__
  std::string const str() const;

private:

  UFLDomain::Type const domain_;
  UFLSpace const space_;

  bool const invalid_;

  uint const geometric_dimension_;
  uint const topological_dimension_;

  /// Geometrical quantities
  UFLCellSurfaceArea const cell_surface_area_;
  UFLCellVolume const cell_volume_;
  UFLCircumradius const circumradius_;
  UFLFacetArea const facet_area_;
  UFLFacetNormal const facet_normal_;
  UFLSpatialCoordinate const x_;

  std::string const repr_;
  std::string const str_;

};

} /* namespace dolfin */
#endif /* __UFL_CELL_H */
