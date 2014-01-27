// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2014-01-21
// Last changed: 2014-01-21

#ifndef __UFL_RESTRICTED_ELEMENT_H_
#define __UFL_RESTRICTED_ELEMENT_H_

#include <dolfin/ufl/UFLFiniteElementBase.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLRestrictedElement
 *
 *  @brief  Provides an interface complying with UFL RestrictedElement.
 */

class RestrictedElement : public FiniteElementBase
{

public:

  ///
  RestrictedElement(FiniteElementBase const& element,
                       Domain::Type const domain);

  ///
  ~RestrictedElement();

  //--- INTERFACE -------------------------------------------------------------

  /// Return whether the basis functions of this element is spatially constant
  /// over each cell
  bool const is_cellwise_constant() const;

  /// Return the symmetry dict, which is a mapping c0 -> c1 meaning that
  /// component c0 is represented by component c1
  std::map<uint, uint> const symmetry() const;

  /// Extract direct subelement index and subelement relative component index
  /// for a given component index
  std::pair<ValueArray, ValueArray> const extract_subelement_component(
      ValueArray const& i) const;

  /// Recursively extract component index relative to a (simple) element and
  /// that element for given value component index
  std::pair<uint, FiniteElementBase const * const> const extract_component(ValueArray const& i) const;

  /// Return number of sub elements
  uint const num_sub_elements() const;

  /// Return list of sub elements
  FiniteElementBaseList const& sub_elements() const;

  //---------------------------------------------------------------------------

  ///
  FiniteElementBase const& element();

  /// __repr__
  std::string const repr() const;

  /// __str__
  std::string const str() const;

protected:

  ValueArray const value_shape_;
  std::map<uint, uint> const symmetry_;
  FiniteElementBase const& element_;

  mutable std::string repr_;
  mutable std::string str_;

};

} /* namespace ufl */
#endif /* __UFL_RESTRICTED_ELEMENT_H_ */
