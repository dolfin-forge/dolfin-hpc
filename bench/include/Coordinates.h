// Copyright (C) 2013-15 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __COORDINATES_H_
#define __COORDINATES_H_

#include <dolfin/function/Function.h>

namespace dolfin
{

class Coordinates : public Function
{

public:

  ///
  Coordinates(Mesh& mesh) :
      Function(mesh),
      gdim_(mesh.geometry().dim())
  {
  }

  ///
  void eval(real * values, real const * x) const
  {
    std::copy(x, x + gdim_, values);
  }

  ///
  uint dim(uint i) const
  {
    return (i == 0 ? gdim_ : 1);
  }

  ///
  uint rank() const
  {
    return 1;
  }

private:

  uint const gdim_;

};

#endif /* __COORDINATES_H_ */

} /* namespace dolfin */
