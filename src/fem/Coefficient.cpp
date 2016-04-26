// Copyright (C) 2014 Aurelien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#include <dolfin/fem/Coefficient.h>

namespace dolfin
{

//-----------------------------------------------------------------------------
uint Coefficient::value_size() const
{
  uint size = 1;
  for (uint i = 0; i < this->rank(); ++i)
  {
    size *= this->dim(i);
  }
  return size;
}
//-----------------------------------------------------------------------------
void Coefficient::interpolate(real* coefficients, const ufc::cell& cell,
                              const ufc::finite_element& finite_element,
                              const Cell& dolfin_cell) const
{
  finite_element.evaluate_dofs(coefficients, *this, cell);
}
//-----------------------------------------------------------------------------
void Coefficient::interpolate(real* coefficients, const ufc::cell& cell,
                              const ufc::finite_element& finite_element,
                              const Cell& dolfin_cell, uint facet) const
{
  finite_element.evaluate_dofs(coefficients, *this, cell);
}
//-----------------------------------------------------------------------------

} /* namespace dolfin */
