// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
// Adapted from licorne.

#ifndef __DOLFIN_FUNCTION_OPERATORS_H_
#define __DOLFIN_FUNCTION_OPERATORS_H_

#include "impl/ops.h"
#include <dolfin/mesh/MeshValues.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
template <class Entity, class Operator, class Value>
inline MeshValues<real, Entity>&
operator<<(MeshValues<real, Entity>& v, EntityOp<Operator, Entity, Value>& o)
{
  for (typename Entity::iterator it(v.mesh()); !it.end(); ++it)
  {
    o(static_cast<Entity const&>(*it), v[it->index()]);
  }
  return v;
}
//-----------------------------------------------------------------------------

}

#endif /* __DOLFIN_FUNCTION_OPERATORS_H_ */
