// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_FINITE_ELEMENT_BASE_H_
#define __UFL_FINITE_ELEMENT_BASE_H_

#include <dolfin/ufl/UFLElementList.h>

#include <dolfin/common/types.h>
#include <dolfin/mesh/CellType.h>

namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLFiniteElement
 *
 *  @brief  Provides an interface complying with UFL FiniteElementBase.
 */

class UFLFiniteElementBase
{

public:

  /// Return finite element family
  virtual UFLElementList::FamilyType const family() const = 0;

  /// Return cell of finite element
  virtual CellType::Type const cell() const = 0;

  /// Return whether the basis functions of this element is spatially constant
  /// over each cell
  virtual bool const is_cellwise_constant() const = 0;

  /// Return polynomial degree of finite element
  /// Present in FIAT interface
  virtual uint const degree() const = 0;

  /// Return quadrature scheme of finite element
  virtual std::string const quadrature_scheme() const = 0;

  /// Return the shape of the value space
  /// Present in FIAT interface
  virtual uint const value_shape() const = 0;

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

  /// Return the domain onto which the element is restricted
  virtual domain_restriction() const = 0;

  /// Return number of sub elements
  virtual uint const num_sub_elements() const = 0;

  /// Return list of sub elements
  virtual Array<UFLFiniteElementBase const *> sub_elements() const = 0;

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

  /// Format as string for signature evaluation
  /// __repr__
  virtual std::string signature() = 0;

protected:

  UFLFiniteElementBase()
  {
  }

  virtual ~UFLFiniteElementBase()
  {
  }

private:

  UFLElementList::FamilyType family_;

  CellType::Type cell_;

};

} /* namespace dolfin */
#endif /* __UFL_FINITE_ELEMENT_BASE_H_ */
