// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-08-23
// Last changed: 2017-10-09

#ifndef __DOLFIN_LIBSIM_INTERFACE_H
#define __DOLFIN_LIBSIM_INTERFACE_H

#include <dolfin/common/Label.h>
#include <dolfin/common/types.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/main/PE.h>
#include <dolfin/main/MPI.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/insitu/libsimPipeline.h>
#include <dolfin/config/dolfin_config.h>

#include <algorithm>
#include <cstring>

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

    static void batchRender();

    static void ctrlLoop();

    static void addData(GenericFunction& function ,std::string name);
    
    static void addData(LabelList<GenericFunction>& functions);

    static void addData(Mesh& mesh);

    static void addPipeline(libsimPipeline& pipeline);

  private:

    static int setupEnv();

#ifdef HAVE_MPI
    static MPI::Communicator comm;
#endif


    // Simulation state (running)
    static int runflag;

    struct libsimData
    {
      double *t_;		
      uint tstep_;
      Mesh *mesh_;      
      LabelList<GenericFunction> function_list_;
      Array<libsimPipeline*> pipelines_;      
      libsimData(): t_(NULL),mesh_(NULL){}
    };

    // Simulation (insitu) data
    static libsimData InsituData_;

#ifdef HAVE_LIBSIM

//--- Callback functions  -----------------------------------------------------

    // Function to return meta data
    inline static visit_handle libsimGetMetaData(void *data) 
    {
      libsimData *d = (libsimData *) data;

      visit_handle md = VISIT_INVALID_HANDLE;
      visit_handle msh = VISIT_INVALID_HANDLE;
      visit_handle vmd = VISIT_INVALID_HANDLE;

      if (VisIt_SimulationMetaData_alloc(&md) == VISIT_OKAY) 
      {

	if (d->t_ != NULL)
	  VisIt_SimulationMetaData_setCycleTime(md, d->tstep_, *(d->t_));

	// Mesh meta data
	if (VisIt_MeshMetaData_alloc(&msh) == VISIT_OKAY)
	{
	  VisIt_MeshMetaData_setName(msh, "Mesh");
	  VisIt_MeshMetaData_setMeshType(msh, VISIT_MESHTYPE_UNSTRUCTURED);
	  VisIt_MeshMetaData_setSpatialDimension(msh, d->mesh_->geometry_dimension());
	  VisIt_MeshMetaData_setTopologicalDimension(msh,
						     d->mesh_->topology_dimension());
	  VisIt_MeshMetaData_setNumDomains(msh, PE::size());
	  VisIt_SimulationMetaData_addMesh(md, msh);
	}
	
	// Function meta data
	for (LabelList<GenericFunction>::iterator it = 
	       d->function_list_.begin(); it != d->function_list_.end(); it++)
	  {
	  if (VisIt_VariableMetaData_alloc(&vmd) == VISIT_OKAY) 
	  {

	    VisIt_VariableMetaData_setName(vmd, it->second.c_str());
	    VisIt_VariableMetaData_setMeshName(vmd, "Mesh");
	    VisIt_VariableMetaData_setCentering(vmd, VISIT_VARCENTERING_NODE);

	    GenericFunction *u = it->first;
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

    // Function to return domains (PE partitions)
    inline static visit_handle libsimGetDomain(const char *name, void *data)
    {
      visit_handle dl = VISIT_INVALID_HANDLE;
      if (VisIt_DomainList_alloc(&dl) == VISIT_OKAY)
      {
	visit_handle hdl;
	int pe_rank = PE::rank();
	int pe_size = PE::size();
	VisIt_VariableData_alloc(&hdl);
	VisIt_VariableData_setDataI(hdl, VISIT_OWNER_COPY, 1, 1, &pe_rank);
	VisIt_DomainList_setDomains(dl, pe_size, hdl);
      }
      
      return dl;
      
    }
    
    // Function to return mesh data
    inline static visit_handle libsimGetMesh(int domain, 
					     const char *name, void *data)
    {
      libsimData *d = (libsimData *) data;

      visit_handle msh = VISIT_INVALID_HANDLE;
      visit_handle coords = VISIT_INVALID_HANDLE;
      visit_handle conn = VISIT_INVALID_HANDLE;

      // Currently we're limited to one mesh
      if (strcmp(name, "Mesh") != 0) return msh;
      
      if (VisIt_UnstructuredMesh_alloc(&msh) == VISIT_OKAY &&
	  VisIt_VariableData_alloc(&coords) == VISIT_OKAY &&
	  VisIt_VariableData_alloc(&conn) == VISIT_OKAY)
      {

	VisIt_VariableData_setDataD(coords, VISIT_OWNER_SIM,
				    d->mesh_->topology_dimension(),
				    d->mesh_->num_vertices(),
				    d->mesh_->geometry().coordinates());
	VisIt_UnstructuredMesh_setCoords(msh, coords);
	
	uint *dolfin_conn = d->mesh_->cells();
	uint cell_type = 0;
	switch(d->mesh_->type().cellType())
	  {
	  case CellType::triangle:
	    cell_type = VISIT_CELL_TRI;
	    break;
	  case CellType::tetrahedron:
	    cell_type = VISIT_CELL_TET;
	  default:
	    error("Unsupported (insitu) mesh cell type");
	    break;
	  }

	int nconn = d->mesh_->num_cells() * 
	  (d->mesh_->type().num_entities(0) + 1);
	int *visit_conn = new int[nconn];
	
	for (int *cp = &visit_conn[0], i = 0; i < d->mesh_->num_cells(); i++)
	{
	  *(cp++) = cell_type;
	  for (int j = 0; j < d->mesh_->type().num_entities(0); j++) 
	  {
	    *(cp++) = (int) *(dolfin_conn++);
	  }
	}
	
	VisIt_VariableData_setDataI(conn, VISIT_OWNER_VISIT, 1, 
				    nconn, visit_conn);
	VisIt_UnstructuredMesh_setConnectivity(msh, d->mesh_->num_cells(), conn);

      }
      
      return msh;
    }

    
    // Function to return function data
    inline static visit_handle libsimGetFunction(int domain, 
						 const char *name, void *data)
    {
      libsimData *d = (libsimData *) data;

      visit_handle func = VISIT_INVALID_HANDLE;


      if (VisIt_VariableData_alloc(&func) == VISIT_OKAY) 
      {
	LabelList<GenericFunction>::iterator it = d->function_list_.begin();
	for( ; it != d->function_list_.end(); it++)
	{
	  if (strcmp(it->second.c_str(), name) == 0) break;
	}

	if (it == d->function_list_.end())
	{
	  VisIt_VariableData_free(func);
	  return VISIT_INVALID_HANDLE;
	}

	GenericFunction *u = it->first;
	uint const num_cell_vertices = d->mesh_->type().num_entities(0);
	uint const num_cell_dofs = num_cell_vertices * u->value_size();
	real *vertex_values = new real[num_cell_dofs * d->mesh_->num_vertices()];

	u->interpolate_vertex_values(vertex_values);

	// TODO Check if VisIt can handle u->value_size() > 3
	real *values = new real[u->value_size() * d->mesh_->num_vertices()];
	memset(values, u->value_size() * d->mesh_->num_vertices(), sizeof(real));
	real *vp = values;

	for (VertexIterator v(*(d->mesh_)); !v.end(); ++v)
	{
	  for (uint i = 0; i < u->value_size(); i++) 
	  {
	    *(vp++) = vertex_values[v->index() + i * d->mesh_->num_vertices()];
	  }
	}
	delete[] vertex_values;
	
	VisIt_VariableData_setDataD(func, VISIT_OWNER_VISIT, u->value_size(),
				    d->mesh_->num_vertices(), values);

      }
      return func;
    }


//-----------------------------------------------------------------------------

#endif

  };
  
  inline void libsimInterface::addData(GenericFunction& function,
				       std::string name)
  {
    Label<GenericFunction> item(function, name);
    InsituData_.function_list_.push_back(item);
  }
  
  inline void libsimInterface::addData(LabelList<GenericFunction>& functions)
  {
    InsituData_.function_list_ = functions;
  }

  inline void libsimInterface::addData(Mesh& mesh)
  {
    InsituData_.mesh_ = &mesh;
  }

  inline void libsimInterface::addPipeline(libsimPipeline& pipeline)
  {
    InsituData_.pipelines_.push_back(&pipeline);
  }

}

#endif
