// Copyright (C) 2007 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_PETSC_FACTORY_H
#define __DOLFIN_PETSC_FACTORY_H

#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/LinearAlgebraFactory.h>
#include <dolfin/la/SparsityPattern.h>

#ifdef HAVE_PETSC

#include <dolfin/la/petsc/PETScMatrix.h>
#include <dolfin/la/petsc/PETScVector.h>

namespace dolfin
{

class PETScFactory : public LinearAlgebraFactory
{
public:
  /// Destructor
  ~PETScFactory() override = default;

  /// Create empty matrix
  auto createMatrix() const -> PETScMatrix * override;

  /// Create empty vector
  auto createVector() const -> PETScVector * override;

  /// Create empty sparsity pattern
  auto createPattern() const -> SparsityPattern * override;

  /// Return singleton instance
  static auto instance() -> PETScFactory &
  {
    return factory;
  }

private:
  /// Private Constructor
  PETScFactory() = default;

  static PETScFactory factory;
};

}

#endif

#endif
