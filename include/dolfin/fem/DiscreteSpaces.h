// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __DOLFIN_DISCRETE_SPACES_H
#define __DOLFIN_DISCRETE_SPACES_H

#include <dolfin/elements/Elements.h>

namespace dolfin
{

//--- SCALAR ------------------------------------------------------------------

template<class T> struct scalar : ufl::FiniteElement
{
  scalar(CellType const& cell, uint const degree) :
      ufl::FiniteElement(T::type, cell, degree)
  {
  }

  static ElementType const type = T::type;

  void disp() const { message("%s", this->repr().c_str()); }
};

typedef scalar<cg> sCG;
typedef scalar<dg> sDG;

//--- VECTOR ------------------------------------------------------------------

template<class T> struct vector : ufl::VectorElement
{
  vector(CellType const& cell, uint const degree, uint const dim) :
      ufl::VectorElement(T::type, cell, degree, dim)
  {
  }

  static ElementType const type = T::type;

  void disp() const { message("%s", this->repr().c_str()); }
};

typedef vector<cg> vCG;
typedef vector<dg> vDG;

//--- ELEMENTS ----------------------------------------------------------------

typedef scalar<cg> CG;
typedef scalar<dg> DG;

} /* namespace dolfin */

#endif /* __DOLFIN_DISCRETE_SPACES_H */
