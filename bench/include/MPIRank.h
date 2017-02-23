// Copyright (C) 2013-15 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __MPIRANK_H_
#define __MPIRANK_H_

#include <dolfin/function/Function.h>

namespace dolfin
{

class MPIRank : public Function
{

public:

  ///
  MPIRank(Mesh& mesh) :
      Function(mesh)
  {
  }

  ///
  void eval(real * values, real const * x) const
  {
    values[0] = dolfin::MPI::rank();
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

};

#endif /* __MPIRANK_H_ */

} /* namespace dolfin */
