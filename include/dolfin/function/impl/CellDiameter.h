#ifndef __DOLFIN_FUNCTION_IMPL_CELL_DIAMETER_H_
#define __DOLFIN_FUNCTION_IMPL_CELL_DIAMETER_H_

#include <dolfin/function/ValueSpace.h>
#include <dolfin/fem/UFCCell.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

struct CellDiameter : public ValueSpace<>
{
  ///
  void evaluate(real* values, const real* x, const UFCCell& cell) const
  {
    values[0] = cell.cell->diameter();
  }

  ///
  real operator()(Cell& cell)
  {
    return cell.diameter();
  }
};

//-----------------------------------------------------------------------------

}

#endif /* __DOLFIN_FUNCTION_IMPL_CELL_DIAMETER_H_ */
