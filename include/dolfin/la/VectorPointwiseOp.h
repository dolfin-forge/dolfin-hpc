// Copyright (C) 2019 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_VECTOR_POINTWISE_OP_H
#define __DOLFIN_VECTOR_POINTWISE_OP_H

namespace dolfin
{

  // List of predefined pointwise operations

  enum VectorPointwiseOp {pw_min, pw_max, pw_mult, pw_div, pw_maxabs};

}

#endif
