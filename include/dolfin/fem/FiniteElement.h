// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FINITE_ELEMENT_H
#define __DOLFIN_FINITE_ELEMENT_H

#include <dolfin/common/Array.h>
#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/ufc/ufc.h>

#include <cstring>
#include <string>

namespace dolfin
{

class CellType;
class Form;
class Mesh;

/**
 *  DOCUMENTATION:
 *
 *  @class  FiniteElement
 *
 *  @brief  This class provides an interface to the UFC definition of a finite
 *          element as well as an extension to manage mixed elements.
 *
 *  @author Aurélien Larcher <larcher@kth.se>
 */

class FiniteElement : public ufc::finite_element
{

public:
  /// Create finite element from given coefficient space of form
  FiniteElement( CellType const & cell, Form & form, size_t const i );

  /// Create finite element from UFC object
  /// Ownership of the UFC object is transfered to the instance if the boolean
  /// is set to true, otherwise a clone of the finite element is created.
  /// In any case the instance the member attribute will be destroyed.
  FiniteElement( ufc::finite_element const & element, bool const owner );

  /// Create the i-th subelement of the given element
  FiniteElement( ufc::finite_element const & element, size_t const i );

  /// Create subelement of the given element for given subsystem
  FiniteElement( ufc::finite_element const & element,
                 Array< size_t > const &     sub_system );

  /// Copy constructor
  explicit FiniteElement( FiniteElement const & other );

  ///
  ~FiniteElement() override;

  /// Check if the element definitions are identical
  auto operator==( FiniteElement const & other ) const -> bool;
  auto operator!=( FiniteElement const & other ) const -> bool;

  //--- UFC INTERFACE ---------------------------------------------------------

  /// Return a string identifying the finite element
  auto signature() const -> const char * override;

  /// Return the cell shape
  auto cell_shape() const -> ufc::shape override;

  /// Return the topological dimension of the cell shape
  auto topological_dimension() const -> size_t override;

  /// Return the geometric dimension of the cell shape
  auto geometric_dimension() const -> size_t override;

  /// Return the dimension of the finite element function space
  auto space_dimension() const -> size_t override;

  /// Return the rank of the value space
  auto value_rank() const -> size_t override;

  /// Return the dimension of the value space for axis i
  auto value_dimension( size_t i ) const -> size_t override;

  /// Return the value size
  auto value_size() const -> size_t override;

  /// Return the rank of the reference value space
  auto reference_value_rank() const -> size_t override;

  /// Return the dimension of the reference value space for axis i
  auto reference_value_dimension( size_t i ) const -> size_t override;

  /// Return the number of components of the reference value space
  auto reference_value_size() const -> size_t override;

  /// Return the maximum polynomial degree of the finite element function space
  auto degree() const -> size_t override;

  /// Return the family of the finite element function space
  auto family() const -> const char * override;

  /// Evaluate all basis functions at given point X in reference cell
  void evaluate_reference_basis( double *       reference_values,
                                 size_t         num_points,
                                 double const * X ) const override;

  /// Evaluate specific order derivatives of all basis functions at given point
  /// X in reference cell
  void evaluate_reference_basis_derivatives( double *       reference_values,
                                             size_t         order,
                                             size_t         num_points,
                                             double const * X ) const override;
  /// Transform order n derivatives (can be 0) of all basis functions
  /// previously evaluated in points X in reference cell with given
  /// Jacobian J and its inverse K for each point
  void transform_reference_basis_derivatives(
    double *       values,
    size_t         order,
    size_t         num_points,
    double const * reference_values,
    double const * X,
    double const * J,
    double const * detJ,
    double const * K,
    int            cell_orientation ) const override;

  /// FIXME has to be reintroduced in ffc/ufc
  // /// Compute mapped coordinates for evaluate_basis()
  // void evaluate_basis_map_coordinates(double & X,
  //                                     double & Y,
  //                                     double & Z,
  //                                     const double* coordinates,
  //                                     const ufc::cell& c) const override;

  /// FIXME has to be reintroduced in ffc/ufc
  // /// Compute mapped coordinates for evaluate_basis()
  // void evaluate_basis_from_coordinates(const double X,
  //                                      const double Y,
  //                                      const double Z,
  //                                      double** values) const override;

  /// Evaluate basis function i at given point in cell
  void evaluate_basis(
    size_t                          i,
    double *                        values,
    double const *                  x,
    double const *                  coordinate_dofs,
    int                             cell_orientation,
    ufc::coordinate_mapping const * cm = nullptr ) const override;

  /// Evaluate all basis functions at given point in cell
  void evaluate_basis_all(
    double *                        values,
    double const *                  x,
    double const *                  coordinate_dofs,
    int                             cell_orientation,
    ufc::coordinate_mapping const * cm = nullptr ) const override;

  /// Evaluate order n derivatives of basis function i at given point in cell
  void evaluate_basis_derivatives(
    size_t                          i,
    size_t                          n,
    double *                        values,
    double const *                  x,
    double const *                  coordinate_dofs,
    int                             cell_orientation,
    ufc::coordinate_mapping const * cm = nullptr ) const override;

  /// Evaluate order n derivatives of all basis functions at given point in cell
  void evaluate_basis_derivatives_all(
    size_t                          n,
    double *                        values,
    double const *                  x,
    double const *                  coordinate_dofs,
    int                             cell_orientation,
    ufc::coordinate_mapping const * cm = nullptr ) const override;

  /// Evaluate linear functional for dof i on the function f
  auto evaluate_dof( size_t                          i,
                     ufc::function const &           f,
                     double const *                  coordinate_dofs,
                     int                             cell_orientation,
                     ufc::cell const &               c,
                     ufc::coordinate_mapping const * cm = nullptr ) const
    -> double override;

  /// Evaluate linear functionals for all dofs on the function f
  void evaluate_dofs(
    double *                        values,
    ufc::function const &           f,
    double const *                  coordinate_dofs,
    int                             cell_orientation,
    ufc::cell const &               c,
    ufc::coordinate_mapping const * cm = nullptr ) const override;

  /// Interpolate vertex values from dof values
  void interpolate_vertex_values(
    double *                        vertex_values,
    double const *                  dof_values,
    double const *                  coordinate_dofs,
    int                             cell_orientation,
    ufc::coordinate_mapping const * cm = nullptr ) const override;

  /// Map coordinate xhat from reference cell to coordinate x in cell
  // void map_from_reference_cell( double *          x,
  //                               const double *    xhat,
  //                               const ufc::cell & c ) const override;

  // /// Map from coordinate x in cell to coordinate xhat in reference cell
  // void map_to_reference_cell( double *          xhat,
  //                             const double *    x,
  //                             const ufc::cell & c ) const override;

  /// Tabulate the coordinates of all dofs on a cell
  void tabulate_dof_coordinates(
    double *                        dof_coordinates,
    double const *                  coordinate_dofs,
    ufc::coordinate_mapping const * cm = nullptr ) const override;

  /// Tabulate the coordinates of all dofs on a reference cell
  void tabulate_reference_dof_coordinates(
    double * reference_dof_coordinates ) const override;

  /// Return the number of sub elements (for a mixed element)
  auto num_sub_elements() const -> size_t override;

  /// Create a new finite element for sub element i (for a mixed element)
  auto create_sub_element( size_t i ) const -> ufc::finite_element * override;

  /// Create a new class instance
  auto create() const -> ufc::finite_element * override;

  //--- EXTENSION OF UFC INTERFACE --------------------------------------------

  /// Recursively extract sub finite element
  static auto create_sub_element( ufc::finite_element const & finite_element,
                                  Array< size_t > const &     sub_system )
    -> ufc::finite_element *;

  /// Create sub finite element of given finite element
  auto create_sub_element( Array< size_t > const & sub_system ) const
    -> ufc::finite_element *;

  /// Get value dimensions for sub spaces just one level down for axis i
  auto sub_value_dimensions( size_t i ) const -> Array< size_t > const &;

  /// Get value dimensions for sub spaces just one level down for axis i
  auto sub_value_offsets( size_t i ) const -> Array< size_t > const &;

  /// Get list of scalar finite elements ordered by entries
  auto flatten() const -> Array< ufc::finite_element const * > const &;

  /// Create flatten representation finite element (append sub elements)
  static void flatten( ufc::finite_element const *            element,
                       Array< ufc::finite_element const * > & stack,
                       size_t                                 maxlevel );

  /// Create flatten representation finite element (append sub elements)
  static void flatten( ufc::finite_element const *            element,
                       Array< ufc::finite_element const * > & stack );

  /// Check if the element can be seen as a vector element
  auto is_vectorizable() const -> bool;

  //---

  void disp() const;

private:
  void Initialize();

  //--- ATTRIBUTES ------------------------------------------------------------
  ufc::finite_element const * ufc_finite_element_;

  //
  Array< size_t > * sub_value_dims_;

  //
  Array< size_t > * sub_value_offs_;

  //
  mutable Array< ufc::finite_element const * > flattened_;
};

//-----------------------------------------------------------------------------

inline auto FiniteElement::signature() const -> const char *
{
  return ufc_finite_element_->signature();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::cell_shape() const -> ufc::shape
{
  return ufc_finite_element_->cell_shape();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::topological_dimension() const -> size_t
{
  return ufc_finite_element_->topological_dimension();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::geometric_dimension() const -> size_t
{
  return ufc_finite_element_->geometric_dimension();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::space_dimension() const -> size_t
{
  return ufc_finite_element_->space_dimension();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::value_rank() const -> size_t
{
  return ufc_finite_element_->value_rank();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::value_dimension( size_t i ) const -> size_t
{
  return ufc_finite_element_->value_dimension( i );
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::value_size() const -> size_t
{
  return ufc_finite_element_->value_size();
}

inline auto FiniteElement::reference_value_rank() const -> size_t
{
  return ufc_finite_element_->reference_value_rank();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::reference_value_dimension( size_t i ) const -> size_t
{
  return ufc_finite_element_->reference_value_dimension( i );
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::reference_value_size() const -> size_t
{
  return ufc_finite_element_->reference_value_size();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::degree() const -> size_t
{
  return ufc_finite_element_->degree();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::family() const -> const char *
{
  return ufc_finite_element_->family();
}

//-----------------------------------------------------------------------------

inline void FiniteElement::evaluate_reference_basis( double * reference_values,
                                                     size_t   num_points,
                                                     double const * X ) const
{
  ufc_finite_element_->evaluate_reference_basis(
    reference_values, num_points, X );
}

//-----------------------------------------------------------------------------

inline void FiniteElement::evaluate_reference_basis_derivatives(
  double *       reference_values,
  size_t         order,
  size_t         num_points,
  double const * X ) const
{
  ufc_finite_element_->evaluate_reference_basis_derivatives(
    reference_values, order, num_points, X );
}

//-----------------------------------------------------------------------------

inline void FiniteElement::transform_reference_basis_derivatives(
  double *       values,
  size_t         order,
  size_t         num_points,
  double const * reference_values,
  double const * X,
  double const * J,
  double const * detJ,
  double const * K,
  int            cell_orientation ) const
{
  ufc_finite_element_->transform_reference_basis_derivatives(
    values,
    order,
    num_points,
    reference_values,
    X,
    J,
    detJ,
    K,
    cell_orientation );
}

//-----------------------------------------------------------------------------

// /// Compute mapped coordinates for evaluate_basis()
// inline void FiniteElement::evaluate_basis_map_coordinates(double & X,
//                                                           double & Y,
//                                                           double & Z,
//                                                           const double*
//                                                           coordinates, const
//                                                           ufc::cell& c) const
// {
//   ufc_finite_element_->evaluate_basis_map_coordinates(X, Y, Z, coordinates,
//   c);
// }

// /// Compute mapped coordinates for evaluate_basis()
// inline void FiniteElement::evaluate_basis_from_coordinates(const double X,
//                                                            const double Y,
//                                                            const double Z,
//                                                            double** values)
//                                                            const
// {
//   ufc_finite_element_->evaluate_basis_from_coordinates(X, Y, Z, values);
// }

//-----------------------------------------------------------------------------

inline void
  FiniteElement::evaluate_basis( size_t         i,
                                 double *       values,
                                 double const * x,
                                 double const * coordinate_dofs,
                                 int            cell_orientation,
                                 ufc::coordinate_mapping const * cm ) const
{
  ufc_finite_element_->evaluate_basis(
    i, values, x, coordinate_dofs, cell_orientation, cm );
}

//-----------------------------------------------------------------------------

inline void
  FiniteElement::evaluate_basis_all( double *       values,
                                     double const * x,
                                     double const * coordinate_dofs,
                                     int            cell_orientation,
                                     ufc::coordinate_mapping const * cm ) const
{
  ufc_finite_element_->evaluate_basis_all(
    values, x, coordinate_dofs, cell_orientation, cm );
}

//-----------------------------------------------------------------------------

inline void FiniteElement::evaluate_basis_derivatives(
  size_t                          i,
  size_t                          n,
  double *                        values,
  double const *                  x,
  double const *                  coordinate_dofs,
  int                             cell_orientation,
  ufc::coordinate_mapping const * cm ) const
{
  ufc_finite_element_->evaluate_basis_derivatives(
    i, n, values, x, coordinate_dofs, cell_orientation, cm );
}

//-----------------------------------------------------------------------------

inline void FiniteElement::evaluate_basis_derivatives_all(
  size_t                          n,
  double *                        values,
  double const *                  x,
  double const *                  coordinate_dofs,
  int                             cell_orientation,
  ufc::coordinate_mapping const * cm ) const
{
  ufc_finite_element_->evaluate_basis_derivatives_all(
    n, values, x, coordinate_dofs, cell_orientation, cm );
}

//-----------------------------------------------------------------------------

inline auto
  FiniteElement::evaluate_dof( size_t                          i,
                               ufc::function const &           f,
                               double const *                  coordinate_dofs,
                               int                             cell_orientation,
                               ufc::cell const &               c,
                               ufc::coordinate_mapping const * cm ) const
  -> double
{
  return ufc_finite_element_->evaluate_dof(
    i, f, coordinate_dofs, cell_orientation, c, cm );
}

//-----------------------------------------------------------------------------

inline void
  FiniteElement::evaluate_dofs( double *              values,
                                ufc::function const & f,
                                double const *        coordinate_dofs,
                                int                   cell_orientation,
                                ufc::cell const &     c,
                                ufc::coordinate_mapping const * cm ) const
{
  ufc_finite_element_->evaluate_dofs(
    values, f, coordinate_dofs, cell_orientation, c, cm );
}

//-----------------------------------------------------------------------------

inline void FiniteElement::interpolate_vertex_values(
  double *                        vertex_values,
  double const *                  dof_values,
  double const *                  coordinate_dofs,
  int                             cell_orientation,
  ufc::coordinate_mapping const * cm ) const
{
  ufc_finite_element_->interpolate_vertex_values(
    vertex_values, dof_values, coordinate_dofs, cell_orientation, cm );
}

//-----------------------------------------------------------------------------

// inline void FiniteElement::map_from_reference_cell( double *          x,
//                                                     const double *    xhat,
//                                                     const ufc::cell & c )
//                                                     const
// {
//   ufc_finite_element_->map_from_reference_cell( x, xhat, c );
// }

//-----------------------------------------------------------------------------

// inline void FiniteElement::map_to_reference_cell( double *          xhat,
//                                                   const double *    x,
//                                                   const ufc::cell & c ) const
// {
//   ufc_finite_element_->map_to_reference_cell( xhat, x, c );
// }

//-----------------------------------------------------------------------------

inline void FiniteElement::tabulate_dof_coordinates(
  double *                        dof_coordinates,
  double const *                  coordinate_dofs,
  ufc::coordinate_mapping const * cm ) const
{
  ufc_finite_element_->tabulate_dof_coordinates(
    dof_coordinates, coordinate_dofs, cm );
}

//-----------------------------------------------------------------------------

inline void FiniteElement::tabulate_reference_dof_coordinates(
  double * reference_dof_coordinates ) const
{
  ufc_finite_element_->tabulate_reference_dof_coordinates(
    reference_dof_coordinates );
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::num_sub_elements() const -> size_t
{
  return ufc_finite_element_->num_sub_elements();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::create_sub_element( size_t i ) const
  -> ufc::finite_element *
{
  return ufc_finite_element_->create_sub_element( i );
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::create() const -> ufc::finite_element *
{
  return ufc_finite_element_->create();
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::operator==( FiniteElement const & other ) const
  -> bool
{
  return ( std::strcmp( this->signature(), other.signature() ) == 0 );
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::operator!=( FiniteElement const & other ) const
  -> bool
{
  return !( *this == other );
}

//-----------------------------------------------------------------------------

inline auto
  FiniteElement::create_sub_element( Array< size_t > const & sub_system ) const
  -> ufc::finite_element *
{
  return FiniteElement::create_sub_element( *ufc_finite_element_, sub_system );
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::sub_value_dimensions( size_t i ) const
  -> Array< size_t > const &
{
  return sub_value_dims_[i];
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::sub_value_offsets( size_t i ) const
  -> Array< size_t > const &
{
  return sub_value_offs_[i];
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::flatten() const
  -> Array< ufc::finite_element const * > const &
{
  if ( flattened_.empty() )
  {
    flatten( ufc_finite_element_, flattened_ );
  }
  return flattened_;
}
}

#endif
