// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_UFC_EXPRESSION_H
#define __DOLFIN_UFC_EXPRESSION_H

#include <dolfin/function/Expression.h>

#include <ufc.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
struct UFCExpression : ufc::function
{
  UFCExpression(Expression const& E) :
    E_(E)
  {
  }

  ///
  inline void evaluate(real* values, const real* coordinates,
                       const ufc::cell&) const override
  {
    E_.eval(values, coordinates);
  }

private:

  Expression const& E_;
};
//-----------------------------------------------------------------------------

} /* namespace dolfin */

#endif /* __DOLFIN_UFC_EXPRESSION_H */
