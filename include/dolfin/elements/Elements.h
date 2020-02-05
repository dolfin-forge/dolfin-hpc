// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_ELEMENTS_H
#define __DOLFIN_ELEMENTS_H

#include <dolfin/ufl/UFLFiniteElement.h>
#include <dolfin/ufl/UFLVectorElement.h>
#include <dolfin/ufl/UFLTensorElement.h>

#include <dolfin/ufl/UFLEnrichedElement.h>
#include <dolfin/ufl/UFLMixedElement.h>
#include <dolfin/ufl/UFLRestrictedElement.h>

namespace dolfin
{

/// Unfortunately UFL is a big abstraction mess.
/// Wrapping is needed to avoid breakage and allow proper testing.

typedef ufl::Family::Type ElementType;

struct Elements
{
  struct cg { static ElementType const type = ufl::Family::CG; };

  struct dg { static ElementType const type = ufl::Family::DG; };
};

} /* namespace dolfin */

#endif /* __DOLFIN_ELEMENTS_H */
