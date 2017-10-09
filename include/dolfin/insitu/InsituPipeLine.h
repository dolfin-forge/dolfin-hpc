// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-10-09
// Last changed: 2017-10-09

#ifndef __DOLFIN_INSITU_PIPELINE_INTERFACE_H
#define __DOLFIN_INSITU_PIPELINE_INTERFACE_H


#include <dolfin/common/types.h>

namespace dolfin
{

/// This class defines a visualization pipeline to be executed by an
/// insitu backend

  class InsituPipeLine
  {
  public:
    
    /// Constructor
    InsituPipeLine();
    
    /// Destructor
    virtual ~InsituPipeLine();
    
    /// Execute visualization pipeline
    virtual void exec() const = 0;
  };
} 

#endif
