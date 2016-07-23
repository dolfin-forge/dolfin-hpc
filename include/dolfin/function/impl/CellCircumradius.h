#ifndef __DOLFIN_FUNCTION_IMPL_CELL_CIRCUMRADIUS_H_
#define __DOLFIN_FUNCTION_IMPL_CELL_CIRCUMRADIUS_H_

#include <dolfin/function/ValueSpace.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

struct CellDiameter : public ValueSpace<>
{
  ///
  void eval(real* values, const real* x, const UFCCell& cell) const
  {
    values[0] = cell.cell->circumradius();
  }

  ///
  real operator()(Cell& cell)
  {
    return cell.circumradius();
  }
};

//-----------------------------------------------------------------------------

}

#endif /* __DOLFIN_FUNCTION_IMPL_CELL_CIRCUMRADIUS_H_ */
