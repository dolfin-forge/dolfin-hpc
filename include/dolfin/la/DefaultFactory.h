// Copyright (C) 2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_DEFAULT_FACTORY_H
#define __DOLFIN_DEFAULT_FACTORY_H

#include <dolfin/la/LinearAlgebraFactory.h>

namespace dolfin
{

class DefaultFactory : public LinearAlgebraFactory
{

public:

  // Return instance of default backend
  static LinearAlgebraFactory& factory();

  /// Constructor
  DefaultFactory() = default;

  /// Destructor
  ~DefaultFactory() override = default;

  /// Create empty matrix
  GenericMatrix* createMatrix() const override;

  /// Create empty vector
  GenericVector* createVector() const override;

  /// Create empty sparsity pattern
  GenericSparsityPattern * createPattern() const override;

};

} /* namespace dolfin */

#endif /* __DOLFIN_DEFAULT_FACTORY_H */
