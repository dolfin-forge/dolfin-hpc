// Copyright (C) 2007 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-11-30
// Last changed: 2007-12-06

#ifndef __PETSC_FACTORY_H
#define __PETSC_FACTORY_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include "PETScMatrix.h"
#include "PETScVector.h"
#include "SparsityPattern.h"
#include "LinearAlgebraFactory.h"

namespace dolfin
{

  class PETScFactory : public LinearAlgebraFactory
  {
  public:

    /// Destructor
    virtual ~PETScFactory() {}

    /// Create empty matrix
    PETScMatrix* createMatrix() const;

    /// Create empty vector
    PETScVector* createVector() const;

    /// Create empty sparsity pattern 
    SparsityPattern* createPattern() const;

    /// Return singleton instance
    static PETScFactory& instance() 
    { return factory; }

  private:

    /// Private Constructor
    PETScFactory() {}
    static PETScFactory factory;

  };

}

#endif

#endif
