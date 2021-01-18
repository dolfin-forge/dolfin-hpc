// Copyright (C) 2017 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_DISCRETE_SPACES_H
#define __DOLFIN_DISCRETE_SPACES_H

#include <dolfin/elements/Elements.h>

namespace dolfin
{

//--- SCALAR ------------------------------------------------------------------

// template<class T> struct scalar : ufl::FiniteElement
// {
//   scalar(CellType const& cell, uint const degree) :
//       ufl::FiniteElement(T::type, cell, degree)
//   {
//   }

//   static ElementType const type = T::type;

//   void disp() const { message("%s", this->repr().c_str()); }
// };

// using sCG = scalar<Elements::cg>;
// using sDG = scalar<Elements::dg>;

// //--- VECTOR ------------------------------------------------------------------

// template<class T> struct vector : ufl::VectorElement
// {
//   vector(CellType const& cell, uint const degree, uint const dim) :
//       ufl::VectorElement(T::type, cell, degree, dim)
//   {
//   }

//   static ElementType const type = T::type;

//   void disp() const { message("%s", this->repr().c_str()); }
// };

// using vCG = vector<Elements::cg>;
// using vDG = vector<Elements::dg>;

// //--- ELEMENTS ----------------------------------------------------------------

// using CG = scalar<Elements::cg>;
// using DG = scalar<Elements::dg>;

} /* namespace dolfin */

#endif /* __DOLFIN_DISCRETE_SPACES_H */
