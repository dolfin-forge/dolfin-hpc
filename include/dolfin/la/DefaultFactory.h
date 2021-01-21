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
  static auto factory() -> LinearAlgebraFactory&;

  /// Constructor
  DefaultFactory() = default;

  /// Destructor
  ~DefaultFactory() override = default;

  /// Create empty matrix
  auto createMatrix() const -> GenericMatrix* override;

  /// Create empty vector
  auto createVector() const -> GenericVector* override;

  /// Create empty sparsity pattern
  auto createPattern() const -> GenericSparsityPattern * override;

};

} /* namespace dolfin */

#endif /* __DOLFIN_DEFAULT_FACTORY_H */
