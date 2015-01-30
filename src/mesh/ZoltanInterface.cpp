// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2015-01-30
// Last changed: 2015-01-30

#include <dolfin/config/dolfin_config.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshRenumber.h>
#include <dolfin/mesh/ZoltanInterface.h>
#include <dolfin/parameter/parameters.h>

#include <dolfin/mesh/Vertex.h>
#include <dolfin/mesh/Cell.h>


#ifdef HAVE_ZOLTAN
#include <zoltan_cpp.h>
#endif


using namespace dolfin;

#ifdef HAVE_ZOLTAN
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionCommonZoltan(Mesh& mesh, 
					 MeshFunction<uint>& partitions,
					 MeshFunction<uint>* weight)
{
  
  zz_ = new Zoltan(MPI::DOLFIN_COMM);

  // General query functions
  zz_->Set_Num_Obj_Fn(partitionZoltanNumObj, &mesh);
  zz_->Set_Obj_List_Fn(partitionZoltanObjList, &mesh);  
  
  /* TODO */

  delete zz_;
}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionGeomZoltan(Mesh& mesh, 
					MeshFunction<uint>& partitions)
{

  zz_ = new Zoltan(MPI::DOLFIN_COMM);

  // General query functions
  zz_->Set_Num_Obj_Fn(partitionZoltanNumObj, &mesh);
  zz_->Set_Obj_List_Fn(partitionZoltanObjList, &mesh);  

  /* TODO */

  delete zz_;
}
//-----------------------------------------------------------------------------
int ZoltanInterface::partitionZoltanNumObj(void *data, int *ierr) 
{

  Mesh *mesh = (Mesh *) data;
  *ierr = ZOLTAN_OK;

  return mesh->numCells();

}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionZoltanObjList(void *data, int num_gid_entries, 
					   int num_lid_entries, 
					   ZOLTAN_ID_PTR global_ids, 
					   ZOLTAN_ID_PTR local_ids, 
					   int wgt_dim, float *obj_wgts,
					   int *ierr) 
{
  Mesh *mesh = (Mesh *) data;
  *ierr = ZOLTAN_OK;
  
  for (uint i = 0; i < mesh->numCells(); i++) 
  {
    global_ids[i] = mesh->distdata().get_global(i, 0);
    local_ids[i] = i;
  }

  return;

}
//-----------------------------------------------------------------------------
int ZoltanInterface::partitionZoltanNumGeom(void *data, int *ierr)
{
  Mesh *mesh = (Mesh *) data;
  *ierr = ZOLTAN_OK;

  return mesh->geometry().dim();

}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionZoltanGeomCoords(void *data, int num_gid_entries,
					      int num_lid_entries, int num_obj, 
					      ZOLTAN_ID_PTR global_ids, 
					      ZOLTAN_ID_PTR local_ids, 
					      int num_dim, double *geom_vec, 
					      int *ierr) 
{
  /* TODO */
}
//-----------------------------------------------------------------------------
#else
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionCommonZoltan(Mesh& mesh, 
					 MeshFunction<uint>& partitions,
					 MeshFunction<uint>* weight)
{
  error("DOLFIN needs to be built with Zoltan support");
}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionGeomZoltan(Mesh& mesh, 
					MeshFunction<uint>& partitions)
{
  error("DOLFIN needs to be built with Zoltan support");
}
//-----------------------------------------------------------------------------
#endif

