// Copyright (C) 2007 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LINEAR_ALGEBRA_FACTORY_H
#define __DOLFIN_LINEAR_ALGEBRA_FACTORY_H

#include "GenericMatrix.h"
#include "GenericSparsityPattern.h"
#include "GenericVector.h"

namespace dolfin
{

  class LinearAlgebraFactory
  {
    public:

    /// Constructor
    LinearAlgebraFactory() = default;

    /// Destructor
    virtual ~LinearAlgebraFactory() = default;

    /// Create empty matrix
    virtual auto createMatrix() const -> dolfin::GenericMatrix* = 0;

    /// Create empty vector
    virtual auto createVector() const -> dolfin::GenericVector* = 0;

    /// Create empty sparsity pattern 
    virtual auto createPattern() const -> dolfin::GenericSparsityPattern * = 0;

  };

}

#endif
