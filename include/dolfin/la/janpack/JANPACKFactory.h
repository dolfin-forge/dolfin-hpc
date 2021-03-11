// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_JANPACK_FACTORY_H
#define __DOLFIN_JANPACK_FACTORY_H

#include <dolfin/la/LinearAlgebraFactory.h>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/la/janpack/JANPACKMat.h>
#include <dolfin/la/janpack/JANPACKVec.h>

#ifdef HAVE_JANPACK

namespace dolfin
{

  class JANPACKFactory : public LinearAlgebraFactory
  {
  public:

    /// Destructor
    virtual ~JANPACKFactory() {}

    // Create empty matrix
    JANPACKMat* createMatrix() const;

    /// Create empty vector
    JANPACKVec* createVector() const;

    /// Create empty sparsity pattern
    SparsityPattern* createPattern() const
    { return new SparsityPattern(); }

    /// Return singleton instance
    static JANPACKFactory& instance()
    { return factory; }

  private:

    /// Private Constructor
    JANPACKFactory() {}
    static JANPACKFactory factory;

  };

}

#endif

#endif
