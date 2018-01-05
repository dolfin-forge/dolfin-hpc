// Copyright (C) 2017 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-08-23
// Last changed: 2018-01-01

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

    enum Mode
    {
      batch, interactive
    };

    static void init(Mode mode, bool debug = false);
        
    static void shutdown();

    static void batchRender(real t = 0.0, uint tstep = 0.0);

    static void ctrlLoop(real t = 0.0, uint tstep = 0.0, int blocking = 0);

    static void addData(GenericFunction& function, 
			std::string function_name,
			std::string mesh_name);
    
    static void addData(Mesh& mesh, std::string name);

    static void addPipeline(libsimPipeline& pipeline);

  private:

    static int initBatch();
    
    static int initInteractive();

    static int setupEnv();

    static int setupCallbacks();

    // Simulation state (running)
    static int runflag;

    struct libsimData
    {
      double t_;		
      uint tstep_;
      bool batch_;
      LabelList<Mesh> mesh_list_;
      LabelList<GenericFunction> function_list_;
      _map<std::string, std::string> function_mesh_map_;
      Array<libsimPipeline*> pipelines_;
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

	VisIt_SimulationMetaData_setCycleTime(md, d->tstep_, d->t_);

	// Mesh meta data
	for (LabelList<Mesh>::iterator it = 
	       d->mesh_list_.begin(); it != d->mesh_list_.end(); it++)
	{
	  if (VisIt_MeshMetaData_alloc(&msh) == VISIT_OKAY)
	  {
	    Mesh *mesh = it->first;
	    VisIt_MeshMetaData_setName(msh, it->second.c_str());
	    VisIt_MeshMetaData_setMeshType(msh, VISIT_MESHTYPE_UNSTRUCTURED);
	    VisIt_MeshMetaData_setSpatialDimension(msh, 
						   mesh->geometry().dim());
	    VisIt_MeshMetaData_setTopologicalDimension(msh,
						       mesh->topology().dim());
	    VisIt_MeshMetaData_setNumDomains(msh, PE::size());
	    VisIt_SimulationMetaData_addMesh(md, msh);
	  }
	}

	// Function meta data
	for (LabelList<GenericFunction>::iterator it = 
	       d->function_list_.begin(); it != d->function_list_.end(); it++)
	{
	  if (VisIt_VariableMetaData_alloc(&vmd) == VISIT_OKAY) 
	  {

	    const char *name = it->second.c_str();
	    const char *mesh_name = d->function_mesh_map_[name].c_str();
	    VisIt_VariableMetaData_setName(vmd, name);
	    VisIt_VariableMetaData_setMeshName(vmd, mesh_name);
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

      if (VisIt_UnstructuredMesh_alloc(&msh) == VISIT_OKAY &&
	  VisIt_VariableData_alloc(&coords) == VISIT_OKAY &&
	  VisIt_VariableData_alloc(&conn) == VISIT_OKAY)
      {
	
	LabelList<Mesh>::iterator it = d->mesh_list_.begin();
	for( ; it != d->mesh_list_.end(); it++)
	{
	  if (strcmp(it->second.c_str(), name) == 0) break;
	}
	
	if (it == d->mesh_list_.end())
	{
	  VisIt_VariableData_free(msh);
	  VisIt_VariableData_free(coords);
	  VisIt_VariableData_free(conn);
	  return VISIT_INVALID_HANDLE;
	}
	
	Mesh *mesh = it->first;
	
	VisIt_VariableData_setDataD(coords, VISIT_OWNER_SIM,
				    mesh->topology().dim(),
				    mesh->num_vertices(),
				    mesh->geometry().coordinates());
	VisIt_UnstructuredMesh_setCoords(msh, coords);
	
	uint *dolfin_conn = mesh->cells();
	uint cell_type = 0;
	switch(mesh->type().cellType())
	  {
	  case CellType::triangle:
	    cell_type = VISIT_CELL_TRI;
	    break;
	  case CellType::tetrahedron:
	    cell_type = VISIT_CELL_TET;
	    break;
	  default:
	    error("Unsupported (insitu) mesh cell type");
	    break;
	  }

	int nconn = mesh->num_cells() * 
	  (mesh->type().num_entities(0) + 1);
	int *visit_conn = new int[nconn];
	
	for (int *cp = &visit_conn[0], i = 0; i < mesh->num_cells(); i++)
	{
	  *(cp++) = cell_type;
	  for (int j = 0; j < mesh->type().num_entities(0); j++) 
	  {
	    *(cp++) = (int) *(dolfin_conn++);
	  }
	}
	
	VisIt_VariableData_setDataI(conn, VISIT_OWNER_VISIT, 1, 
				    nconn, visit_conn);
	VisIt_UnstructuredMesh_setConnectivity(msh, mesh->num_cells(), conn);

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

	// Fetch assoicated mesh
	const char *mesh_name = d->function_mesh_map_[name].c_str();
	LabelList<Mesh>::iterator mesh_it = d->mesh_list_.begin();
	for( ; mesh_it != d->mesh_list_.end(); mesh_it++)
	{
	  if (strcmp(mesh_it->second.c_str(), mesh_name) == 0) break;
	}

	if (mesh_it == d->mesh_list_.end())
	{
	  VisIt_VariableData_free(func);
	  return VISIT_INVALID_HANDLE;
	}

	GenericFunction *u = it->first;
	Mesh *mesh = mesh_it->first;
	uint const num_cell_vertices = mesh->type().num_entities(0);
	uint const num_cell_dofs = num_cell_vertices * u->value_size();
	real *vertex_values = new real[num_cell_dofs * mesh->num_vertices()];

	u->interpolate_vertex_values(vertex_values);

	// TODO Check if VisIt can handle u->value_size() > 3
	real *values = new real[u->value_size() * mesh->num_vertices()];
	real *vp = values;

	for (VertexIterator v(*(mesh)); !v.end(); ++v)
	{
	  for (uint i = 0; i < u->value_size(); i++) 
	  {
	    *(vp++) = vertex_values[v->index() + i * mesh->num_vertices()];
	  }
	}
	delete[] vertex_values;
	
	VisIt_VariableData_setDataD(func, VISIT_OWNER_VISIT, u->value_size(),
				    mesh->num_vertices(), values);

      }
      return func;
    }


//-----------------------------------------------------------------------------

#endif

  };
  
  inline void libsimInterface::addData(GenericFunction& function,
				       std::string function_name, 
				       std::string mesh_name)
  {
    Label<GenericFunction> item(function, function_name);
    InsituData_.function_list_.push_back(item);
    
    // Check if the mesh is already registered...
    LabelList<Mesh>::iterator it = InsituData_.mesh_list_.begin();
    for( ; it != InsituData_.mesh_list_.end(); it++)
    {
      if (it->second == mesh_name) 
      {
	if (it->first->hash() == function.mesh().hash())
	{      
	  break;
	}
	else
	{
	  error("A different mesh with the same name is already registered");
	}
      }
    }

    // ...if not, register it with libsim
    if (it == InsituData_.mesh_list_.end())
    {
      addData(function.mesh(), mesh_name);
    }

    InsituData_.function_mesh_map_[function_name] = mesh_name;


  }
  
  inline void libsimInterface::addData(Mesh& mesh, std::string name)
  {
    Label<Mesh> item(mesh, name);
    InsituData_.mesh_list_.push_back(item);
  }

  inline void libsimInterface::addPipeline(libsimPipeline& pipeline)
  {
    InsituData_.pipelines_.push_back(&pipeline);
  }

}

#endif
