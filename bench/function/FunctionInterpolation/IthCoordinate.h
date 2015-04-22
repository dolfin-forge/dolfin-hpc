// Copyright (C) 2013-15 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __ITH_COORDINATE_H_
#define __ITH_COORDINATE_H_

#include <dolfin/function/Function.h>

namespace dolfin
{

class IthCoordinate : public Function
{

public:

  ///
  IthCoordinate(Mesh& mesh, uint i) :
      Function(mesh),
      i_(i)
  {
    if (i >= mesh.geometry().dim())
    {
      error("IthCoordinate: invalid coordinate %d given geometry in %dd.", i,
            mesh.geometry().dim());
    }
  }

  ///
  void eval(real * values, real const * x) const
  {
    values[0] = x[i_];
  }

  ///
  uint dim(uint i) const
  {
    return 1;
  }

  ///
  uint rank() const
  {
    return 0;
  }

private:

  uint const i_;
};

#endif /* __ITH_COORDINATE_H_ */

} /* namespace dolfin */
