// Copyright (C) 2021 Julian Hornich
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_TRILINOS_FACTORY_H
#define __DOLFIN_TRILINOS_FACTORY_H

#include <dolfin/config/dolfin_config.h>
#include <dolfin/la/LinearAlgebraFactory.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/trilinos/TrilinosMatrix.h>
#include <dolfin/la/trilinos/TrilinosVector.h>

#ifdef HAVE_TRILINOS

namespace dolfin
{

namespace trilinos
{

class Factory : public LinearAlgebraFactory
{
public:
  /// Destructor
  ~Factory() override = default;

  /// Create empty matrix
  auto createMatrix() const -> trilinos::Matrix * override;

  /// Create empty vector
  auto createVector() const -> trilinos::Vector * override;

  /// Create empty sparsity pattern
  auto createPattern() const -> SparsityPattern * override;

  /// Return singleton instance
  static auto instance() -> Factory &
  {
    return factory;
  }

private:
  /// Private Constructor
  Factory() = default;

  static Factory factory;
};

} // namespace trilinos

} // namespace dolfin

#endif // HAVE_TRILINOS

#endif // __DOLFIN_TRILINOS_FACTORY_H
