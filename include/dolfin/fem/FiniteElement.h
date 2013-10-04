// Copyright (C) 2013 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2013-09-12
// Last changed: 2013-09-12

#ifndef __FINITE_ELEMENT_H
#define __FINITE_ELEMENT_H

#include <dolfin/config/dolfin_config.h>
#include <dolfin/elements/ElementLibrary.h>
#include <dolfin/fem/Form.h>

#include <ufc.h>

#include <cstring>

namespace dolfin
{

class Mesh;

///

class FiniteElement
{

  /// Dirty trick before ensuring that there is no critical performance hit
  friend class DiscreteFunction;

public:

  FiniteElement(std::string const& signature);

  FiniteElement(Mesh& mesh, Form& form, uint i);

  FiniteElement(ufc::finite_element& finite_element,
                bool const finite_element_local);

  ~FiniteElement();

  //--- INTERFACE -----------------------------------------------------------
  // Implements UFC v1.1, extension to v2.1.1 forseeable.
  /// Return a string identifying the finite element
  const char* signature() const;

  /// Return the cell shape
  ufc::shape cell_shape() const;

  /// Return the topological dimension of the cell shape (UFC v2.1.1)
  unsigned int topological_dimension() const;

  /// Return the geometric dimension of the cell shape  (UFC v2.1.1)
  unsigned int geometric_dimension() const;

  /// Return the dimension of the finite element function space
  unsigned int space_dimension() const;

  /// Return the rank of the value space
  unsigned int value_rank() const;

  /// Return the dimension of the value space for axis i
  unsigned int value_dimension(unsigned int i) const;

  /// Evaluate basis function i at given point in cell
  void evaluate_basis(unsigned int i, double* values, const double* coordinates,
                      const ufc::cell& c) const;

//#ifndef UFC_BACKWARD_COMPATIBILITY
//  /// Evaluate all basis functions at given point in cell
//  virtual void evaluate_basis_all(double* values, const double* coordinates,
//                                  const cell& c) const = 0;
//#else
//  /// Evaluate all basis functions at given point in cell
//  virtual void evaluate_basis_all(double* values,
//      const double* coordinates,
//      const cell& c) const
//  { throw std::runtime_error("Not implemented (introduced in UFC v1.1).");}
//#endif
//
  /// Evaluate order n derivatives of basis function i at given point in cell
  void evaluate_basis_derivatives(unsigned int i, unsigned int n,
                                  double* values, const double* coordinates,
                                  const ufc::cell& c) const;

//#ifndef UFC_BACKWARD_COMPATIBILITY
//  /// Evaluate order n derivatives of all basis functions at given point in cell
//  virtual void evaluate_basis_derivatives_all(unsigned int n, double* values,
//                                              const double* coordinates,
//                                              const cell& c) const = 0;
//#else
//  /// Evaluate order n derivatives of all basis functions at given point in cell
//  virtual void evaluate_basis_derivatives_all(unsigned int n,
//      double* values,
//      const double* coordinates,
//      const cell& c) const
//  { throw std::runtime_error("Not implemented (introduced in UFC v1.1).");}
//#endif

  /// Evaluate linear functional for dof i on the function f
  double evaluate_dof(unsigned int i, const ufc::function& f,
                      const ufc::cell& c) const;

  /// Evaluate linear functionals for all dofs on the function f
  void evaluate_dofs(double* values, const ufc::function& f,
                     const ufc::cell& c) const;

  /// Interpolate vertex values from dof values
  void interpolate_vertex_values(double* vertex_values,
                                 const double* dof_values,
                                 const ufc::cell& c) const;

//#ifndef UFC_BACKWARD_COMPATIBILITY
//  // omitted for backward compatibility code -------------
//  /// Map coordinate xhat from reference cell to coordinate x in cell
//  virtual void map_from_reference_cell(double* x, const double* xhat,
//                                       const cell& c) const = 0;
//
//  /// Map from coordinate x in cell to coordinate xhat in reference cell
//  virtual void map_to_reference_cell(double* xhat, const double* x,
//                                     const cell& c) const = 0;
//
//  // end omit ---------------------------------------------
//#endif

  /// Return the number of sub elements (for a mixed element)
  unsigned int num_sub_elements() const;

  /// Create a new finite element for sub element i (for a mixed element)
  ufc::finite_element* create_sub_element(unsigned int i) const;

  // Recursively extract sub finite element
  static ufc::finite_element* create_sub_element(
      ufc::finite_element const& finite_element, Array<uint> const& sub_system);

//#ifndef UFC_BACKWARD_COMPATIBILITY
//  // omitted for backward compatibility code -------------
//  /// Create a new class instance
//  virtual finite_element* create() const = 0;
//  // end omit ---------------------------------------------
//#endif

  //--- EXTENSION OF UFC INTERFACE --------------------------------------------

  /// Create sub finite element of given finite element
  ufc::finite_element* create_sub_element(Array<uint> const& sub_system) const;

  /// Get value dimensions for sub spaces just one level down for axis i
  Array<uint> const& sub_value_dimensions(uint i) const;

  /// Get value dimensions for sub spaces just one level down for axis i
  Array<uint> const& sub_value_offsets(uint i) const;

#if ENABLE_UFL

  /// Returns the family of the finite element
  std::string const& family() const;

  /// Returns the type of the finite element
  std::string const& type() const;

  /// Returns the degree of the finite element
  uint const degree() const;

  void info() const;

#endif

private:

  void init();

  //--- ATTRIBUTES ------------------------------------------------------------
  mutable ufc::finite_element * ufc_finite_element_;
  bool const finite_element_local_;
  Array<uint> * sub_value_dims_;
  Array<uint> * sub_value_offs_;

  std::string type_;
  std::string family_;
  std::string strshape_;
  uint topo_dim_;
  uint geom_dim_;
  uint degree_;

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
inline unsigned int FiniteElement::topological_dimension() const
{
  return topo_dim_;
}

//-----------------------------------------------------------------------------
inline unsigned int FiniteElement::geometric_dimension() const
{
  return geom_dim_;
}

//-----------------------------------------------------------------------------
inline unsigned int FiniteElement::space_dimension() const
{
  return ufc_finite_element_->space_dimension();
}

//-----------------------------------------------------------------------------
inline unsigned int FiniteElement::value_rank() const
{
  return ufc_finite_element_->value_rank();
}

//-----------------------------------------------------------------------------
inline unsigned int FiniteElement::value_dimension(unsigned int i) const
{
  return ufc_finite_element_->value_dimension(i);
}

//-----------------------------------------------------------------------------
inline void FiniteElement::evaluate_basis(unsigned int i, double* values,
                                          const double* coordinates,
                                          const ufc::cell& c) const
{
  ufc_finite_element_->evaluate_basis(i, values, coordinates, c);
}

//-----------------------------------------------------------------------------
inline void FiniteElement::evaluate_basis_derivatives(unsigned int i,
                                                      unsigned int n,
                                                      double* values,
                                                      const double* coordinates,
                                                      const ufc::cell& c) const
{
  ufc_finite_element_->evaluate_basis_derivatives(i, n, values, coordinates, c);
}

//-----------------------------------------------------------------------------
inline double FiniteElement::evaluate_dof(unsigned int i,
                                          const ufc::function& f,
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
inline unsigned int FiniteElement::num_sub_elements() const
{
  return ufc_finite_element_->num_sub_elements();
}

//-----------------------------------------------------------------------------
inline ufc::finite_element* FiniteElement::create_sub_element(
    unsigned int i) const
{
  return ufc_finite_element_->create_sub_element(i);
}

#if ENABLE_UFL

//-----------------------------------------------------------------------------
inline std::string const& FiniteElement::family() const
{
  return family_;
}

//-----------------------------------------------------------------------------
inline std::string const& FiniteElement::type() const
{
  return type_;
}

//-----------------------------------------------------------------------------
inline uint const FiniteElement::degree() const
{
  if (degree_ == -1)
  {
    error("Degree is not supported for MixedElement");
  }
  return degree_;
}

#endif

}

#endif
