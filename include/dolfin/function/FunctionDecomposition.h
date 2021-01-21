// Copyright (C) 2015 Aurélien Larcher.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_FUNCTION_DECOMPOSITION_H
#define __DOLFIN_FUNCTION_DECOMPOSITION_H

#include <dolfin/common/Array.h>

namespace dolfin
{

class Function;

class FunctionDecomposition
{

public:
  ///
  static auto compute(Function const& F) -> Array<Function *>;
};

} /* namespace dolfin */

#endif /* __DOLFIN_FUNCTION_DECOMPOSITION_H */
