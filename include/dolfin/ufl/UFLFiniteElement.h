// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_FINITE_ELEMENT_H
#define __DOLFIN_UFL_FINITE_ELEMENT_H

#include <dolfin/ufl/UFLFiniteElementSpace.h>

#include <dolfin/mesh/CellType.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  FiniteElement
 *
 *  @brief  Provides an interface complying with UFL FiniteElement.
 */

class FiniteElement : public FiniteElementSpace
{

public:

  //--- UFL -------------------------------------------------------------------
  ///
  FiniteElement(Family::Type family, Cell const& cell,
                dolfin::uint const degree);

  ///
  explicit FiniteElement(repr_t const& repr);

  ///
  ~FiniteElement() override;

  //--- INTERFACE -------------------------------------------------------------

  /// Return finite element family type
  Family const& family() const override;

  /// Return the metatype of the finite element (C++ only)
  Family::Type metatype() const override;

  /// Return cell of finite element
  Cell const& cell() const override;

  /// Return polynomial degree of finite element
  /// Present in FIAT interface
  type<dolfin::uint> const& degree() const override;

  /// Return the shape of the value space
  /// Present in FIAT interface
  ValueArray const& value_shape() const override;

  /// Return whether the basis functions of this element is spatially constant
  /// over each cell
  bool is_cellwise_constant() const override;

  /// Return the symmetry dict, which is a mapping c0 -> c1 meaning that
  /// component c0 is represented by component c1
  dolfin::_ordered_map<dolfin::uint, dolfin::uint> const& symmetry() const override;

  /// Extract direct subelement index and subelement relative component index
  /// for a given component index
  std::pair<ValueArray, ValueArray> extract_subelement_component(
      ValueArray const& i) const override;

  /// Recursively extract component index relative to a (simple) element and
  /// that element for given value component index
  std::pair<dolfin::uint, FiniteElementSpace const *> extract_component(
      ValueArray const& i) const override;

  /// Return number of sub elements
  dolfin::uint num_sub_elements() const override;

  /// Return list of sub elements
  List const& sub_elements() const override;

  /// __repr__
  repr_t const& repr() const override;

  /// __str__
  std::string const& str() const override;

protected:

  Family const family_;
  Cell const cell_;
  type<dolfin::uint> const degree_;
  ValueArray const value_shape_;
  dolfin::_ordered_map<dolfin::uint, dolfin::uint> const symmetry_;
  FiniteElementSpace::List const sub_elements_;

  repr_t const repr_;
  mutable std::string str_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_FINITE_ELEMENT_H */
