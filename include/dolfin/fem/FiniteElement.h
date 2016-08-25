// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-09-12 (merged from branch larcher)
// Last changed: 2014-09-13

#ifndef __DOLFIN_FINITE_ELEMENT_H
#define __DOLFIN_FINITE_ELEMENT_H

#include <dolfin/config/dolfin_config.h>
#include <dolfin/common/types.h>
#include <dolfin/common/Array.h>

#include <ufc.h>

#include <cstring>

namespace ufl
{

class FiniteElementSpace;

}

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
  FiniteElement(CellType& cell, Form& form, uint const i);

  /// Create finite element from UFC object
  /// Ownership of the UFC object is transfered to the instance if the boolean
  /// is set to true, otherwise a clone of the finite element is created.
  /// In any case the instance the member attribute will be destroyed.
  FiniteElement(ufc::finite_element const& element, bool const owner);

  /// Create the i-th subelement of the given element
  FiniteElement(ufc::finite_element const& element, uint const i);

  /// Create subelement of the given element for given subsystem
  FiniteElement(ufc::finite_element const& element, Array<uint> const& sub_system);

  /// Create finite element from UFL object from pregenerated element (testing)
  explicit FiniteElement(ufl::FiniteElementSpace const& finite_element);

  /// Copy constructor
  explicit FiniteElement(FiniteElement const& other);

  ///
  ~FiniteElement();

  /// Check if the element definitions are identical
  bool operator ==(FiniteElement const& other) const;
  bool operator !=(FiniteElement const& other) const;

  //--- INTERFACE -----------------------------------------------------------
  /// @implements UFC v2.1.3-hpc

  /// Return a string identifying the finite element
  /// UFC @since 1.1
  const char* signature() const;

  /// Return the cell shape
  /// UFC @since 1.1
  ufc::shape cell_shape() const;

  /// Return the topological dimension of the cell shape
  /// UFC @since 2.1.1
  uint topological_dimension() const;

  /// Return the geometric dimension of the cell shape
  /// UFC @since 2.1.1
  uint geometric_dimension() const;

  /// Return the dimension of the finite element function space
  /// UFC @since 1.1
  uint space_dimension() const;

  /// Return the rank of the value space
  /// UFC @since 1.1
  uint value_rank() const;

  /// Return the dimension of the value space for axis i
  /// UFC @since 1.1
  uint value_dimension(uint i) const;

  /// Evaluate basis function i at given point in cell
  /// UFC @since 1.1
  void evaluate_basis(uint i, double* values, const double* coordinates,
                      const ufc::cell& c) const;

  /// Evaluate all basis functions at given point in cell
  /// UFC @since 1.1 but not implemented
  void evaluate_basis_all(double* values, const double* coordinates,
                          const ufc::cell& c) const;

  /// Evaluate order n derivatives of basis function i at given point in cell
  /// UFC @since 1.1
  void evaluate_basis_derivatives(uint i, uint n, double* values,
                                  const double* coordinates,
                                  const ufc::cell& c) const;

  /// Evaluate order n derivatives of all basis functions at given point in cell
  /// UFC @since 1.1 but not implemented
  void evaluate_basis_derivatives_all(uint n, double* values,
                                      const double* coordinates,
                                      const ufc::cell& c) const;

  /// Evaluate linear functional for dof i on the function f
  /// UFC @since 1.1
  double evaluate_dof(uint i, const ufc::function& f, const ufc::cell& c) const;

  /// Evaluate linear functionals for all dofs on the function f
  /// UFC @since 1.1 but not implemented
  void evaluate_dofs(double* values, const ufc::function& f,
                     const ufc::cell& c) const;

  /// Interpolate vertex values from dof values
  /// UFC @since 1.1
  void interpolate_vertex_values(double* vertex_values,
                                 const double* dof_values,
                                 const ufc::cell& c) const;


  /// Map coordinate xhat from reference cell to coordinate x in cell
  /// UFC @since 2.1.1
  void map_from_reference_cell(double* x, const double* xhat,
                               const ufc::cell& c) const;

  /// Map from coordinate x in cell to coordinate xhat in reference cell
  /// UFC @since 2.1.1
  void map_to_reference_cell(double* xhat, const double* x,
                             const ufc::cell& c) const;

  /// Return the number of sub elements (for a mixed element)
  /// UFC @since 1.1 + FIAT + UFL
  uint num_sub_elements() const;

  /// Create a new finite element for sub element i (for a mixed element)
  /// UFC @since 1.1
  ufc::finite_element* create_sub_element(uint i) const;

  /// Create a new class instance
  /// UFC @since 2.1.1
  ufc::finite_element* create() const;

  //--- EXTENSION OF UFC INTERFACE --------------------------------------------

  /// Return the value size
  uint value_size() const;

  /// Recursively extract sub finite element
  static ufc::finite_element* create_sub_element(
      ufc::finite_element const& finite_element, Array<uint> const& sub_system);

  /// Create sub finite element of given finite element
  ufc::finite_element* create_sub_element(Array<uint> const& sub_system) const;

  /// Get value dimensions for sub spaces just one level down for axis i
  Array<uint> const& sub_value_dimensions(uint i) const;

  /// Get value dimensions for sub spaces just one level down for axis i
  Array<uint> const& sub_value_offsets(uint i) const;

  /// Get list of scalar finite elements ordered by entries
  Array<ufc::finite_element const *> const& flatten() const;

  /// Create flatten representation finite element (append sub elements)
  static void flatten(ufc::finite_element const * element,
                      Array<ufc::finite_element const *>& stack,
                      uint maxlevel);

  /// Create flatten representation finite element (append sub elements)
  static void flatten(ufc::finite_element const * element,
                      Array<ufc::finite_element const *>& stack);

  /// Check if the element can be seen as a vector element
  bool is_vectorizable() const;

  //---

  void disp() const;

private:

  void Initialize();

  //--- ATTRIBUTES ------------------------------------------------------------
  ufc::finite_element const * ufc_finite_element_;

  //
  Array<uint> * sub_value_dims_;

  //
  Array<uint> * sub_value_offs_;

  //
  mutable Array<ufc::finite_element const *> flattened_;

};

//-----------------------------------------------------------------------------
inline const char* FiniteElement::signature() const
{
  return ufc_finite_element_->signature();
}

//-----------------------------------------------------------------------------
inline ufc::shape FiniteElement::cell_shape() const
{
  return ufc_finite_element_->cell_shape();
}

//-----------------------------------------------------------------------------
inline uint FiniteElement::topological_dimension() const
{
  return ufc_finite_element_->topological_dimension();
}

//-----------------------------------------------------------------------------
inline uint FiniteElement::geometric_dimension() const
{
  return ufc_finite_element_->geometric_dimension();
}

//-----------------------------------------------------------------------------
inline uint FiniteElement::space_dimension() const
{
  return ufc_finite_element_->space_dimension();
}

//-----------------------------------------------------------------------------
inline uint FiniteElement::value_rank() const
{
  return ufc_finite_element_->value_rank();
}

//-----------------------------------------------------------------------------
inline uint FiniteElement::value_dimension(uint i) const
{
  return ufc_finite_element_->value_dimension(i);
}

//-----------------------------------------------------------------------------
inline void FiniteElement::evaluate_basis(uint i, double* values,
                                          const double* coordinates,
                                          const ufc::cell& c) const
{
  ufc_finite_element_->evaluate_basis(i, values, coordinates, c);
}

//-----------------------------------------------------------------------------
inline void FiniteElement::evaluate_basis_all(double* values,
                                              const double* coordinates,
                                              const ufc::cell& c) const
{
  ufc_finite_element_->evaluate_basis_all(values, coordinates, c);
}

//-----------------------------------------------------------------------------
inline void FiniteElement::evaluate_basis_derivatives(uint i, uint n,
                                                      double* values,
                                                      const double* coordinates,
                                                      const ufc::cell& c) const
{
  ufc_finite_element_->evaluate_basis_derivatives(i, n, values, coordinates, c);
}

//-----------------------------------------------------------------------------
inline void FiniteElement::evaluate_basis_derivatives_all(
    uint n, double* values, const double* coordinates, const ufc::cell& c) const
{
  ufc_finite_element_->evaluate_basis_derivatives_all(n, values, coordinates,
                                                      c);
}

//-----------------------------------------------------------------------------
inline double FiniteElement::evaluate_dof(uint i, const ufc::function& f,
                                          const ufc::cell& c) const
{
  return ufc_finite_element_->evaluate_dof(i, f, c);
}

//-----------------------------------------------------------------------------
inline void FiniteElement::evaluate_dofs(double* values, const ufc::function& f,
                                         const ufc::cell& c) const
{
  ufc_finite_element_->evaluate_dofs(values, f, c);
}

//-----------------------------------------------------------------------------
inline void FiniteElement::interpolate_vertex_values(double* vertex_values,
                                                     const double* dof_values,
                                                     const ufc::cell& c) const
{
  ufc_finite_element_->interpolate_vertex_values(vertex_values, dof_values, c);
}

//-----------------------------------------------------------------------------
inline void FiniteElement::map_from_reference_cell(double* x,
                                                   const double* xhat,
                                                   const ufc::cell& c) const
{
  ufc_finite_element_->map_from_reference_cell(x, xhat, c);
}

//-----------------------------------------------------------------------------
inline void FiniteElement::map_to_reference_cell(double* xhat, const double* x,
                                                 const ufc::cell& c) const
{
  ufc_finite_element_->map_to_reference_cell(xhat, x, c);
}

//-----------------------------------------------------------------------------
inline uint FiniteElement::num_sub_elements() const
{
  return ufc_finite_element_->num_sub_elements();
}

//-----------------------------------------------------------------------------
inline ufc::finite_element* FiniteElement::create_sub_element(uint i) const
{
  return ufc_finite_element_->create_sub_element(i);
}

//-----------------------------------------------------------------------------
inline ufc::finite_element* FiniteElement::create() const
{
  return ufc_finite_element_->create();
}

}

#endif
