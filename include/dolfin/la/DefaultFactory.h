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
  DefaultFactory()
  {
  }

  /// Destructor
  virtual ~DefaultFactory()
  {
  }

  /// Create empty matrix
  virtual GenericMatrix* createMatrix() const;

  /// Create empty vector
  virtual GenericVector* createVector() const;

  /// Create empty sparsity pattern
  virtual GenericSparsityPattern * createPattern() const;

};

} /* namespace dolfin */

#endif /* __DOLFIN_DEFAULT_FACTORY_H */
