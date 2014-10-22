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
  RestrictedElement(FiniteElementBase const& element, Domain const& domain);

  ///
  ~RestrictedElement();

  //--- INTERFACE -------------------------------------------------------------

  /// Return finite element family type
  Family const& family() const;

  /// Return the metatype of the finite element (C++ only)
  Family::Type metatype() const;

  /// Return cell of finite element
  Cell const& cell() const;

  /// Return polynomial degree of finite element
  /// Present in FIAT interface
  type<dolfin::uint> const& degree() const;

  /// Return the shape of the value space
  /// Present in FIAT interface
  ValueArray const& value_shape() const;

  /// Return whether the basis functions of this element is spatially constant
  /// over each cell
  bool const is_cellwise_constant() const;

  /// Return the symmetry dict, which is a mapping c0 -> c1 meaning that
  /// component c0 is represented by component c1
  std::map<dolfin::uint, dolfin::uint> const symmetry() const;

  /// Extract direct subelement index and subelement relative component index
  /// for a given component index
  std::pair<ValueArray, ValueArray> const extract_subelement_component(
      ValueArray const& i) const;

  /// Recursively extract component index relative to a (simple) element and
  /// that element for given value component index
  std::pair<dolfin::uint, FiniteElementBase const *> const extract_component(
      ValueArray const& i) const;

  /// Return number of sub elements
  dolfin::uint const num_sub_elements() const;

  /// Return list of sub elements
  List const& sub_elements() const;

  //---------------------------------------------------------------------------

  ///
  FiniteElementBase const& element();

  /// __repr__
  repr_t const repr() const;

  /// __str__
  std::string const str() const;

protected:

  FiniteElementBase const& element_;
  Family const family_;
  Cell const cell_;
  type<dolfin::uint> const degree_;
  ValueArray const value_shape_;
  std::map<dolfin::uint, dolfin::uint> const symmetry_;

  mutable repr_t repr_;
  mutable std::string str_;

};

} /* namespace ufl */
#endif /* __UFL_RESTRICTED_ELEMENT_H_ */
