// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_FINITE_ELEMENT_H_
#define __UFL_FINITE_ELEMENT_H_

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

class UFLFiniteElement
{

public:

  enum FamilyType
  {
    Lagrange, DG, BDM
  };

  ///
  virtual FamilyType const family() const = 0;

  ///
  virtual CellType const cell() const = 0;

  ///
  virtual bool const is_cellwise_constant() const = 0;

  /// +FIAT
  virtual uint const degree() const = 0;

  ///
  virtual std::string const quadrature_scheme() const = 0;

  /// +FIAT
  virtual uint const value_shape() const = 0;

  ///
  virtual std::map<uint, uint> const symmetry() const = 0;

  ///
  virtual std::pair<uint, uint> const extract_subelement_component(
      uint i) const = 0;

  ///
  virtual uint const extract_component(uint i) const = 0;

  ///
  virtual domain_restriction() const = 0;

  /// +FIAT
  virtual uint const num_sub_elements() const = 0;

  ///
  virtual Array<UFLFiniteElement const *> sub_elements() const = 0;

  /// Operator: equality

  /// Operator: addition

  /// Operator: multiplication

protected:

  UFLFiniteElement();

  virtual ~UFLFiniteElement();

};

} /* namespace dolfin */
#endif /* __UFL_FINITE_ELEMENT_H_ */
