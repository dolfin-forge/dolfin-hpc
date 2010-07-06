// Copyright (C) 2010 Niclas Jansson
// Licensed under the GNU LGPL Version 2.1.
//

#ifndef __JANPACK_FACTORY_H
#define __JANPACK_FACTORY_H

#include "JANPACKMat.h"
#include "JANPACKVec.h"
#include "SparsityPattern.h"
#include "LinearAlgebraFactory.h"

#ifdef HAS_JANPACK

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
