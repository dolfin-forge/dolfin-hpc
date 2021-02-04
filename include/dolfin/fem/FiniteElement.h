// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FINITE_ELEMENT_H
#define __DOLFIN_FINITE_ELEMENT_H

#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/ufc/ufc.h>

#include <cstring>
#include <string>

namespace dolfin
{

class CellType;
class Form;
class Mesh;

/**
 *  DOCUMENTATION:
 *
 *  @class  FiniteElement
 *
 *  @brief  This class provides an interface to the UFC definition of a finite
 *          element as well as an extension to manage mixed elements.
 *
 *  @author Aurélien Larcher <larcher@kth.se>
 */

class FiniteElement
{

public:
  /// Create finite element from given coefficient space of form
  FiniteElement( CellType const & cell, Form & form, size_t const i );

  /// Create finite element from UFC object
  /// Ownership of the UFC object is transfered to the instance if the boolean
  /// is set to true, otherwise a clone of the finite element is created.
  /// In any case the instance the member attribute will be destroyed.
  FiniteElement( ufc::finite_element const & element, bool const owner );

  /// Create the i-th subelement of the given element
  FiniteElement( ufc::finite_element const & element, size_t const i );

  /// Create subelement of the given element for given subsystem
  FiniteElement( ufc::finite_element const &   element,
                 std::vector< size_t > const & sub_system );

  /// Copy constructor
  explicit FiniteElement( FiniteElement const & other );

  ///
  ~FiniteElement();

  /// Check if the element definitions are identical
  auto operator==( FiniteElement const & other ) const -> bool;
  auto operator!=( FiniteElement const & other ) const -> bool;

  auto operator()() const -> ufc::finite_element const & { return *ufc_finite_element_; }

  //--- EXTENSION OF UFC INTERFACE --------------------------------------------

  /// Recursively extract sub finite element
  static auto create_sub_element( ufc::finite_element const &   finite_element,
                                  std::vector< size_t > const & sub_system )
    -> ufc::finite_element *;

  /// Create sub finite element of given finite element
  auto create_sub_element( std::vector< size_t > const & sub_system ) const
    -> ufc::finite_element *;

  /// Get value dimensions for sub spaces just one level down for axis i
  auto sub_value_dimensions( size_t i ) const -> std::vector< size_t > const &;

  /// Get value dimensions for sub spaces just one level down for axis i
  auto sub_value_offsets( size_t i ) const -> std::vector< size_t > const &;

  /// Get list of scalar finite elements ordered by entries
  auto flatten() const -> std::vector< ufc::finite_element const * > const &;

  /// Create flatten representation finite element (append sub elements)
  static void flatten( ufc::finite_element const *                  element,
                       std::vector< ufc::finite_element const * > & stack,
                       size_t                                       maxlevel );

  /// Create flatten representation finite element (append sub elements)
  static void flatten( ufc::finite_element const *                  element,
                       std::vector< ufc::finite_element const * > & stack );

  /// Check if the element can be seen as a vector element
  auto is_vectorizable() const -> bool;

  //---

  void disp() const;

private:
  void Initialize();

  //--- ATTRIBUTES ------------------------------------------------------------
  ufc::finite_element const * ufc_finite_element_;

  //
  std::vector< size_t > * sub_value_dims_;

  //
  std::vector< size_t > * sub_value_offs_;

  //
  mutable std::vector< ufc::finite_element const * > flattened_;
};

//-----------------------------------------------------------------------------

inline auto FiniteElement::operator!=( FiniteElement const & other ) const
  -> bool
{
  return !( *this == other );
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::create_sub_element(
  std::vector< size_t > const & sub_system ) const -> ufc::finite_element *
{
  return FiniteElement::create_sub_element( *ufc_finite_element_, sub_system );
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::sub_value_dimensions( size_t i ) const
  -> std::vector< size_t > const &
{
  return sub_value_dims_[i];
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::sub_value_offsets( size_t i ) const
  -> std::vector< size_t > const &
{
  return sub_value_offs_[i];
}

//-----------------------------------------------------------------------------

inline auto FiniteElement::flatten() const
  -> std::vector< ufc::finite_element const * > const &
{
  if ( flattened_.empty() )
  {
    flatten( ufc_finite_element_, flattened_ );
  }
  return flattened_;
}

} // namespace dolfin

#endif
