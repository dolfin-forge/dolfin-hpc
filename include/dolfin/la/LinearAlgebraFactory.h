// Copyright (C) 2007 Ola Skavhaug.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LINEAR_ALGEBRA_FACTORY_H
#define __DOLFIN_LINEAR_ALGEBRA_FACTORY_H

namespace dolfin
{

class GenericMatrix;
class GenericVector;
class GenericSparsityPattern;

//-----------------------------------------------------------------------------

class LinearAlgebraFactory
{
public:
  /// Constructor
  LinearAlgebraFactory() = default;

  /// Destructor
  virtual ~LinearAlgebraFactory() = default;

  /// Create empty matrix
  virtual auto createMatrix() const -> GenericMatrix * = 0;

  /// Create empty vector
  virtual auto createVector() const -> GenericVector * = 0;

  /// Create empty sparsity pattern
  virtual auto createPattern() const -> GenericSparsityPattern * = 0;
};

//-----------------------------------------------------------------------------

} // namespace dolfin

#endif
