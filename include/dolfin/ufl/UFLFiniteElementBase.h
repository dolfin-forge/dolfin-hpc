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

#include <dolfin/common/types.h>

#include <vector>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLFiniteElement
 *
 *  @brief  Provides an interface complying with UFL FiniteElementBase.
 */

///TODO: Implement quadrature scheme
typedef std::string UFLQuadratureScheme;
typedef std::vector<uint> ValueShape;

class UFLFiniteElementBase : public UFLClass
{

public:

  /// Return finite element family
  UFLElementList::FamilyType const family() const;

  /// Return cell of finite element
  UFLCell const cell() const;

  /// Return polynomial degree of finite element
  /// Present in FIAT interface
  uint const degree() const;

  /// Return quadrature scheme of finite element
  UFLQuadratureScheme const quadrature_scheme() const;

  /// Return the shape of the value space
  /// Present in FIAT interface
  ValueShape const value_shape() const;

  /// Return the domain onto which the element is restricted
  //UFLDomain::Type const domain_restriction() const;

  //--- INTERFACE -------------------------------------------------------------

  /// Return whether the basis functions of this element is spatially constant
  /// over each cell
  virtual bool const is_cellwise_constant() const = 0;

  /// Return the symmetry dict, which is a mapping c0 -> c1 meaning that
  /// component c0 is represented by component c1
  virtual std::map<uint, uint> const symmetry() const = 0;

  /// Extract direct subelement index and subelement relative component index
  /// for a given component index
  virtual std::pair<uint, uint> const extract_subelement_component(
      uint i) const = 0;

  /// Recursively extract component index relative to a (simple) element and
  /// that element for given value component index
  virtual uint const extract_component(uint i) const = 0;

  /// Return number of sub elements
  virtual uint const num_sub_elements() const = 0;

  /// Return list of sub elements
  virtual std::vector<UFLFiniteElementBase const *> const& sub_elements() const = 0;

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
  virtual std::string const repr() const = 0;

  /// __str__
  virtual std::string const str() const = 0;

protected:

  ///
  UFLFiniteElementBase(UFLElementList::FamilyType family,
                       UFLCell const& cell,
                       uint const degree,
                       UFLQuadratureScheme quad_scheme = "None",
                       ValueShape value_shape = ValueShape());

  ///
  virtual ~UFLFiniteElementBase();

  ///
  bool component_is_valid(std::vector<uint> const i);

  ///
  UFLCell const get_cell(std::vector<UFLFiniteElementBase const *> const& elements);

  ///
  uint const get_degree_max(std::vector<UFLFiniteElementBase const *> const& elements);

private:

  UFLElementList::FamilyType const family_;
  UFLCell const cell_;
  uint const degree_;
  UFLQuadratureScheme const quad_scheme_;
  ValueShape const value_shape_;

};

} /* namespace dolfin */
#endif /* __UFL_FINITE_ELEMENT_BASE_H_ */
