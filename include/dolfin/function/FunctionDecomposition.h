// Copyright (C) 2015 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:
// Last changed:

#ifndef __FUNCTION_DECOMPOSITION_H_
#define __FUNCTION_DECOMPOSITION_H_

#include <dolfin/common/Array.h>

namespace dolfin
{

class Function;

class FunctionDecomposition
{

public:

  ///
  FunctionDecomposition(Function const& F, Array<Function *>& Si);

  ///
  static Array<Function *> compute(Function const& F);

private:

  Function const& F_;
  Array<Function *>& Si_;

};

} /* namespace dolfin */

#endif /* __FUNCTION_DECOMPOSITION_H_ */
