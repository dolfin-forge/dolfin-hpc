// Copyright (C) 2006 Garth N. Wells
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-05-29
// Last changed: 

#include <dolfin/config/dolfin_config.h>

#ifndef NO_UBLAS

#ifndef __UBLAS_DENSE_MATRIX_H
#define __UBLAS_DENSE_MATRIX_H

#include "uBlasMatrix.h"

namespace dolfin
{

  typedef uBlasMatrix<ublas_dense_matrix> uBlasDenseMatrix;

}

#endif

#endif
