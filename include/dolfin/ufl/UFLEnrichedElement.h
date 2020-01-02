// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFL_ENRICHED_ELEMENT_H
#define __DOLFIN_UFL_ENRICHED_ELEMENT_H

#include <dolfin/ufl/UFLFiniteElementSpace.h>
#include <dolfin/ufl/UFLCell.h>

namespace ufl
{

/**
 *  DOCUMENTATION:
 *
 *  @class  UFLEnrichedElement
 *
 *  @brief  Provides an interface complying with UFL FiniteElement.
 */

class EnrichedElement : public FiniteElementSpace
{

public:

  ///
  EnrichedElement(List const& elements);

  ///
  ~EnrichedElement();

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
  bool is_cellwise_constant() const;

  /// Return the symmetry dict, which is a mapping c0 -> c1 meaning that
  /// component c0 is represented by component c1
  std::map<dolfin::uint, dolfin::uint> const& symmetry() const;

  /// Extract direct subelement index and subelement relative component index
  /// for a given component index
  std::pair<ValueArray, ValueArray> extract_subelement_component(
      ValueArray const& i) const;

  /// Recursively extract component index relative to a (simple) element and
  /// that element for given value component index
  std::pair<dolfin::uint, FiniteElementSpace const *> extract_component(
      ValueArray const& i) const;

  /// Return number of sub elements
  dolfin::uint num_sub_elements() const;

  /// Return list of sub elements
  List const& sub_elements() const;

  /// __repr__
  repr_t const& repr() const;

  /// __str__
  std::string const& str() const;

protected:

  FiniteElementSpace::List const sub_elements_;
  Family const family_;
  Cell const cell_;
  type<dolfin::uint> const degree_;
  ValueArray const value_shape_;
  std::map<dolfin::uint, dolfin::uint> const symmetry_;

  mutable repr_t repr_;
  mutable std::string str_;

};

} /* namespace ufl */
#endif /* __DOLFIN_UFL_ENRICHED_ELEMENT_H */
