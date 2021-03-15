// Copyright (C) 2014 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_MESH_SPACE_H
#define __DOLFIN_MESH_SPACE_H

#include <dolfin/log/log.h>

#include <typeinfo>

namespace dolfin
{

//-----------------------------------------------------------------------------

struct Space
{
  static size_t const MAX_DIMENSION = 3;

  /// Constructor
  Space() = default;

  /// Destructor
  virtual ~Space() = default;

  /// Equality
  auto operator==( Space const & other ) const -> bool
  {
    if ( typeid( *this ) != typeid( other ) )
      return false;
    return ( this->dim() == other.dim() );
  }

  /// Non-equality
  auto operator!=( Space const & other ) const -> bool
  {
    return !( *this == other );
  }

  /// Space dimension
  virtual auto dim() const -> size_t = 0;

  /// Clone pattern
  virtual auto clone() const -> Space * = 0;

  /// Display info
  virtual void disp() const = 0;
};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_SPACE_H */
