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
  static uint const MAX_DIMENSION = 3;

  /// Constructor
  Space() {}

  /// Destructor
  virtual ~Space() {}

  /// Equality
  bool operator==(Space const& other) const
  {
    if(typeid(*this) != typeid(other)) return false;
    return (this->dim() == other.dim());
  }

  /// Non-equality
  bool operator!=(Space const& other) const
  { return !(*this == other); }

  /// Space dimension
  virtual uint dim() const = 0;

  /// Clone pattern
  virtual Space* clone() const = 0;

  /// Display info
  virtual void disp() const = 0;

};

//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_MESH_SPACE_H */
