// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-08-23
// Last changed: 2017-08-23

#ifndef __DOLFIN_LIBSIM_INTERFACE_H
#define __DOLFIN_LIBSIM_INTERFACE_H

#include <dolfin/common/Label.h>
#include <dolfin/common/types.h>
#include <dolfin/function/Function.h>
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

    static void ctrlLoop();

    static void addData(Function& function ,std::string name);
    
    static void addData(LabelList<Function>& functions);

  private:

    static int setupEnv();

#ifdef HAVE_MPI
    static MPI::Communicator comm;
#endif

    /// Simulation state (running)
    static int runflag;

    static Mesh& mesh_;

    static LabelList<Function> function_list_;

  };

  inline void libsimInterface::addData(Function& function ,std::string name)
  {
    Label<Function> item(function, name);
    function_list_.push_back(item);
  }
  
  inline void libsimInterface::addData(LabelList<Function>& functions)
  {
    function_list_ = functions;
  }
  


}

#endif
