#ifndef __DOLFIN_FUNCTION_IMPL_CELL_VOLUME_H_
#define __DOLFIN_FUNCTION_IMPL_CELL_VOLUME_H_

#include <dolfin/function/ValueSpace.h>

namespace dolfin
{

//-----------------------------------------------------------------------------

struct CellVolume : public ValueSpace<>
{
  ///
  void evaluate(real* values, const real* x, const UFCCell& cell) const
  {
    values[0] = cell.cell->volume();
  }

  ///
  real operator()(Cell& cell)
  {
    return cell.volume();
  }
};

//-----------------------------------------------------------------------------

}

#endif /* __DOLFIN_FUNCTION_IMPL_CELL_VOLUME_H_ */
