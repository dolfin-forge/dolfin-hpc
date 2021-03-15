// Copyright (C) 2013 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FINITE_ELEMENT_H
#define __DOLFIN_FINITE_ELEMENT_H

#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>

#include <ufc.h>

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
  /// Create finite element from UFC object
  explicit FiniteElement( ufc::finite_element const & ufc_element );

  /// Create the i-th subelement of the given element
  explicit FiniteElement( ufc::finite_element const & ufc_element,
                          size_t const i );

  /// Create subelement of the given element for given subsystem
  FiniteElement( ufc::finite_element const &   ufc_element,
                 std::vector< size_t > const & sub_system );

  /// Copy constructor
  FiniteElement( FiniteElement const & copy );

  /// Move constructor
  FiniteElement( FiniteElement && move ) = delete;

  /// Destructor
  ~FiniteElement();

  auto operator=( FiniteElement const & copy ) = delete;
  auto operator=( FiniteElement && move ) = delete;

  /// Check if the element definitions are identical
  auto operator==( FiniteElement const & other ) const -> bool;
  auto operator!=( FiniteElement const & other ) const -> bool;

  /// access ufc interface
  auto ufc() const -> ufc::finite_element const & { return *element; }

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
                       std::vector< ufc::finite_element const * > & stack );

  /// Check if the element can be seen as a vector element
  auto is_vectorizable() const -> bool;

  //---

  void disp() const;

private:
  ufc::finite_element const * element{ nullptr };

public:
  std::string const signature{};
  std::string const family{};

  ufc::shape const shape{ ufc::shape::vertex };

  size_t const degree{ DOLFIN_SIZE_T_UNDEF };
  size_t const tdim{ DOLFIN_SIZE_T_UNDEF };
  size_t const gdim{ DOLFIN_SIZE_T_UNDEF };
  size_t const space_dim{ DOLFIN_SIZE_T_UNDEF };
  size_t const num_sub_elements{ DOLFIN_SIZE_T_UNDEF };

  size_t const                value_rank{ DOLFIN_SIZE_T_UNDEF };
  size_t const                value_size{ DOLFIN_SIZE_T_UNDEF };
  std::vector< size_t > const value_dims{};

  size_t const                ref_value_rank{ DOLFIN_SIZE_T_UNDEF };
  size_t const                ref_value_size{ DOLFIN_SIZE_T_UNDEF };
  std::vector< size_t > const ref_value_dims{};

private:
  //--- ATTRIBUTES ------------------------------------------------------------

  std::vector< std::vector< size_t > >       sub_value_dims_{};
  std::vector< std::vector< size_t > >       sub_value_offs_{};
  std::vector< ufc::finite_element const * > flattened_{};
};

//-----------------------------------------------------------------------------

inline auto FiniteElement::operator==( FiniteElement const & other ) const
  -> bool
{
  return other.signature == this->signature;
}

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
  return FiniteElement::create_sub_element( *element, sub_system );
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
  return flattened_;
}

} // namespace dolfin

#endif
