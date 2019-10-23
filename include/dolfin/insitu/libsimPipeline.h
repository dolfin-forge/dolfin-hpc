// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __DOLFIN_LIBSIM_PIPELINE_INTERFACE_H
#define __DOLFIN_LIBSIM_PIPELINE_INTERFACE_H


#include <dolfin/common/types.h>
#include <dolfin/config/dolfin_config.h>
//#include <dolfin/insitu/InsituPipeLine.h>

#ifdef HAVE_LIBSIM
#include <VisItControlInterface_V2.h>
#include <VisItDataInterface_V2.h>
#endif


namespace dolfin
{

  /// This class defines a visualization pipeline to be executed by libsim
  
  class libsimPipeline 
  {

  public:
    

    libsimPipeline() 
    {
    }

    ~libsimPipeline() {}

    /// Execute visualization pipeline
    virtual void exec(real t, uint tstep) const = 0;

  };
  
} /* namespace dolfin */

#endif
