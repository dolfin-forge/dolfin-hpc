// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-05-24
// Last changed: 2017-05-24

#ifndef __DOLFIN_LIBSIM_INTERFACE_H
#define __DOLFIN_LIBSIM_INTERFACE_H

#include <dolfin/common/types.h>
#include <dolfin/main/MPI.h>
#include <dolfin/config/dolfin_config.h>



#ifdef HAVE_LIBSIM
#include <VisItControlInterface_V2.h>
#endif


namespace dolfin
{
  /// This class provides an interface to VisIt/libsim

  class libsimInterface
  {
  public:
        
    static void initBatch();
    
    static void initInteractive();

    static void shutdown();

    static void batchRender(std::string filename);

  private:

    static int setupEnv();

#ifdef HAVE_MPI
    static MPI::Communicator comm;
#endif

  };

}

#endif
