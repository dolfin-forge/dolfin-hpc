// Copyright (C) 2007 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_FACTORY_H
#define __DOLFIN_PETSC_FACTORY_H

#include <dolfin/config/dolfin_config.h>

#ifdef HAVE_PETSC

#include <dolfin/la/LinearAlgebraFactory.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/PETScVector.h>
#include <dolfin/la/SparsityPattern.h>

namespace dolfin
{

class PETScFactory : public LinearAlgebraFactory
{
public:
  /// Destructor
  virtual ~PETScFactory()
  {
  }

  /// Create empty matrix
  PETScMatrix * createMatrix() const;

  /// Create empty vector
  PETScVector * createVector() const;

  /// Create empty sparsity pattern
  SparsityPattern * createPattern() const;

  /// Return singleton instance
  static PETScFactory & instance()
  {
    return factory;
  }

private:
  /// Private Constructor
  PETScFactory()
  {
  }

  static PETScFactory factory;
};

}

#endif

#endif
