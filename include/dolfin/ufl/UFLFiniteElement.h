// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_FINITE_ELEMENT_H_
#define __UFL_FINITE_ELEMENT_H_

#include <dolfin/ufl/UFLCell.h>
#include <dolfin/ufl/UFLClass.h>
#include <dolfin/ufl/UFLElementList.h>

#include <dolfin/common/types.h>


namespace dolfin
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLFiniteElement
 *
 *  @brief  Provides an interface complying with UFL FiniteElementBase.
 */

using UFLElementList::FamilyType;

class UFLFiniteElementBase : public UFLClass
{

public:

  /// Return finite element family
  virtual FamilyType const family() const = 0;

  /// Return cell of finite element
  virtual UFLCell const cell() const = 0;

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
  /// TODO: Unimplemented

  /// Add two elements, creating an enriched element
  /// Operator: addition
  /// TODO: Unimplemented

  /// Multiply two elements, creating a mixed element
  /// Operator: multiplication
  /// TODO: Unimplemented

  /// Format as string for signature evaluation


protected:

  UFLFiniteElementBase();

  virtual ~UFLFiniteElementBase();

};

} /* namespace dolfin */
#endif /* __UFL_FINITE_ELEMENT_H_ */
