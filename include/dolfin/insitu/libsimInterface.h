// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-08-23
// Last changed: 2017-10-02

#ifndef __DOLFIN_LIBSIM_INTERFACE_H
#define __DOLFIN_LIBSIM_INTERFACE_H

#include <dolfin/common/Label.h>
#include <dolfin/common/types.h>
#include <dolfin/function/Function.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/config/dolfin_config.h>


#ifdef HAVE_LIBSIM
#include <VisItControlInterface_V2.h>
#include <VisItDataInterface_V2.h>
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


    // Simulation state (running)
    static int runflag;

    struct libsimData 
    {
      double& t_;		
      uint tstep_;
      Mesh& mesh_;      
      LabelList<Function> function_list_;      
    };

    // Simulation (insitu) data
    static libsimData InsituData_;

//--- Callback functions  -----------------------------------------------------

    // Function to return meta data
    inline visit_handle libsimGetMetaData(void *data) 
    {
      libsimData *d = (libsimData *) data;

      visit_handle md = VISIT_INVALID_HANDLE;
      visit_handle msh = VISIT_INVALID_HANDLE;
      visit_handle vmd = VISIT_INVALID_HANDLE;

      if (VisIt_SimulationMetaData_alloc(&md) == VISIT_OKAY) 
      {
	VisIt_SimulationMetaData_setCycleTime(md, d->tstep_, d->t_);

	// Mesh meta data
	if (VisIt_MeshMetaData_alloc(&msh) == VISIT_OKAY)
	{
	  VisIt_MeshMetaData_setName(msh, "Mesh");
	  VisIt_MeshMetaData_setMeshType(msh, VISIT_MESHTYPE_UNSTRUCTURED);
	  VisIt_MeshMetaData_setSpatialDimension(msh, d->mesh_.geometry().dim());
	  VisIt_MeshMetaData_setTopologicalDimension(msh,
						     d->mesh_.topology().dim());
	  
	  VisIt_SimulationMetaData_addMesh(md, msh);
	}

	// Function meta data
	for (LabelList<Function>::iterator it = d->function_list_.begin(); 
	     it != d->function_list_.end(); it++)
	{
	  if (VisIt_VariableMetaData_alloc(&vmd) == VISIT_OKAY) 
	  {

	    VisIt_VariableMetaData_setName(vmd, it->second.c_str());
	    VisIt_VariableMetaData_setMeshName(vmd, "Mesh");
	    VisIt_VariableMetaData_setCentering(vmd, VISIT_VARCENTERING_NODE);

	    Function *u = it->first;
	    if (u->value_size() == 1)
	    {
	      VisIt_VariableMetaData_setType(vmd, VISIT_VARTYPE_SCALAR);
	    }
	    else if (u->value_size() > 1)
	    {
	      VisIt_VariableMetaData_setType(vmd, VISIT_VARTYPE_VECTOR);
	    }
	    else 
	    {
	      error("Invalid function");
	    }
	    
	    VisIt_SimulationMetaData_addVariable(md, vmd);
	  
	  }
	}

      }
      return md;
    }

//-----------------------------------------------------------------------------

  };
  
  inline void libsimInterface::addData(Function& function ,std::string name)
  {
    Label<Function> item(function, name);
    InsituData_.function_list_.push_back(item);
  }
  
  inline void libsimInterface::addData(LabelList<Function>& functions)
  {
    InsituData_.function_list_ = functions;
  }

}

#endif
