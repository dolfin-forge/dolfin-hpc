// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_FINITE_ELEMENT_BASE_H_
#define __UFL_FINITE_ELEMENT_BASE_H_

#include <dolfin/ufl/UFLClass.h>

#include <dolfin/ufl/UFLCell.h>
#include <dolfin/ufl/UFLElementList.h>
#include <dolfin/ufl/UFLFamily.h>
#include <dolfin/ufl/UFLQuadratureScheme.h>
#include <dolfin/ufl/UFLValueArray.h>

#include <dolfin/common/types.h>

#include <vector>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLFiniteElement
 *
 *  @brief  Provides an interface complying with UFL FiniteElementBase.
 *
 */

class FiniteElementBase : public Class
{

public:

  //
  typedef std::vector<FiniteElementBase const *> List;

  /// Stupid factory function
  static FiniteElementBase * create(Object::repr_t const repr);

  ///
  virtual ~FiniteElementBase();

  /// Return quadrature scheme of finite element
  QuadratureScheme const& quadrature_scheme() const;

  /// Return the domain onto which the element is restricted
  //Domain::Type const domain_restriction() const; Not implemented

  ///
  dolfin::uint value_size() const;

  //--- INTERFACE -------------------------------------------------------------

  /// Return finite element family type
  virtual Family const& family() const = 0;

  /// Return the metatype of the finite element (C++ only)
  virtual Family::Type metatype() const = 0;

  /// Return cell of finite element
  virtual Cell const& cell() const = 0;

  /// Return polynomial degree of finite element
  /// Present in FIAT interface
  virtual type<dolfin::uint> const& degree() const = 0;

  /// Return the shape of the value space
  /// Present in FIAT interface
  virtual ValueArray const& value_shape() const = 0;

  /// Return whether the basis functions of this element is spatially constant
  /// over each cell
  virtual bool const is_cellwise_constant() const = 0;

  /// Return the symmetry dict, which is a mapping c0 -> c1 meaning that
  /// component c0 is represented by component c1
  virtual std::map<dolfin::uint, dolfin::uint> const symmetry() const = 0;

  /// Extract direct subelement index and subelement relative component index
  /// for a given component index
  virtual std::pair<ValueArray, ValueArray> const extract_subelement_component(
      ValueArray const& i) const = 0;

  /// Recursively extract component index relative to a (simple) element and
  /// that element for given value component index
  virtual std::pair<dolfin::uint, FiniteElementBase const *> const extract_component(
      ValueArray const& i) const = 0;

  /// Return number of sub elements
  virtual dolfin::uint const num_sub_elements() const = 0;

  /// Return list of sub elements
  virtual List const& sub_elements() const = 0;

  /// Operator: equality
  /// __eq__
  /// TODO: Unimplemented

  /// Add two elements, creating an enriched element
  /// Operator: addition
  /// __add__
  /// TODO: Unimplemented

  /// Multiply two elements, creating a mixed element
  /// Operator: multiplication
  /// __mult__
  /// TODO: Unimplemented

  /// __repr__
  virtual repr_t const& repr() const = 0;

  /// __str__
  virtual std::string const& str() const = 0;

  ///
  virtual void display() const;

protected:

  /// Base constructor for finite element definition
  FiniteElementBase(std::string const& name, QuadratureScheme quad_scheme =
                        QuadratureScheme());

  /// Representation-based constructor
  FiniteElementBase(std::string const& name, repr_t repr);

  ///
  bool component_is_valid(ValueArray const& i) const;

  ///
  void check_component(ValueArray const& i) const;

  ///
  Cell const get_cell(List const& elements) const;

  ///
  dolfin::uint const get_degree_max(List const& elements) const;

  /// Get the value size from subelements
  dolfin::uint get_value_size(List const& elements) const;

private:

  QuadratureScheme const quad_scheme_;

};

} /* namespace ufl */
#endif /* __UFL_FINITE_ELEMENT_BASE_H_ */
