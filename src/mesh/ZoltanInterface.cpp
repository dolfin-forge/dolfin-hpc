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

#ifndef ENABLE_P1_OPTIMIZATIONS
  Zoltan *zz_ = new Zoltan(MPI::DOLFIN_COMM);

  // Setup query functions for graph based partitioning
  zz_->Set_Num_Obj_Fn(partitionZoltanNumCells, &mesh);
  zz_->Set_Obj_List_Fn(partitionZoltanCellList, &mesh);

  // Use Zoltan's Parallel Hypergraph and Graph partitioner
  zz_->Set_Param( "LB_METHOD", "GRAPH");

  partitions.init(mesh, mesh.topology().dim());
  partitions = MPI::processNumber();

  partitionZoltanInternal(mesh, partitions, zz_);

  delete zz_;

#else
  error("Zoltan's graph partitioner needs shared face connectivity");
#endif
}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionGeomZoltan(Mesh& mesh,
					  MeshFunction<uint>& partitions)
{

  Zoltan *zz_ = new Zoltan(MPI::DOLFIN_COMM);

  // Setup query functions for geometry based partitioning
  zz_->Set_Num_Obj_Fn(partitionZoltanNumVertices, &mesh);
  zz_->Set_Obj_List_Fn(partitionZoltanVertexList, &mesh);
  zz_->Set_Num_Geom_Fn(partitionZoltanNumGeom, &mesh);
  zz_->Set_Geom_Multi_Fn(partitionZoltanGeomCoords, &mesh);

  // Use a Hilbert Space-Filling Curve
  zz_->Set_Param( "LB_METHOD", "HSFC");

  partitions.init(mesh, 0);
  partitions = MPI::processNumber();

  partitionZoltanInternal(mesh, partitions, zz_);

  delete zz_;
}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionZoltanInternal(Mesh& mesh,
					      MeshFunction<uint>& partitions,
					      Zoltan *zz_)
{
  ZOLTAN_ID_PTR import_global_ids, import_local_ids;
  ZOLTAN_ID_PTR export_global_ids, export_local_ids;

  int changes, num_gid_entries, num_lid_entries;
  int num_import, num_export;

  int *import_procs, *import_to_part;
  int *export_procs, *export_to_part;


  if (zz_->LB_Partition(changes, num_gid_entries, num_lid_entries, num_import,
			import_global_ids, import_local_ids,
			import_procs, import_to_part, num_export,
			export_global_ids, export_local_ids,
			export_procs, export_to_part) != ZOLTAN_OK)
  {
    error("Zoltan failed to partition the mesh");
  }

  // Fill meshfunction from partitions
  for (uint i = 0; i < num_export; i++)
  {
    partitions.set(export_local_ids[i], export_procs[i]);
  }

  Zoltan::LB_Free_Part(&import_global_ids, &import_local_ids,
		       &import_procs, &import_to_part);

  Zoltan::LB_Free_Part(&export_global_ids, &export_local_ids,
		       &export_procs, &export_to_part);

}
//-----------------------------------------------------------------------------
int ZoltanInterface::partitionZoltanNumCells(void *data, int *ierr)
{

  Mesh *mesh = (Mesh *) data;
  *ierr = ZOLTAN_OK;

  return mesh->numCells();

}
//-----------------------------------------------------------------------------
int ZoltanInterface::partitionZoltanNumVertices(void *data, int *ierr)
{

  Mesh *mesh = (Mesh *) data;
  *ierr = ZOLTAN_OK;

  return mesh->numVertices();

}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionZoltanCellList(void *data, int num_gid_entries,
					      int num_lid_entries,
					      ZOLTAN_ID_PTR global_ids,
					      ZOLTAN_ID_PTR local_ids,
					      int wgt_dim, float *obj_wgts,
					      int *ierr)
{
  Mesh *mesh = (Mesh *) data;
  *ierr = ZOLTAN_OK;

  uint i = 0;
  for (CellIterator cell(*mesh); !cell.end(); ++cell)
  {
    global_ids[i] = mesh->distdata().get_global(*cell);
    local_ids[i++] = cell->index();
  }

  return;

}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionZoltanVertexList(void *data,
						int num_gid_entries,
						int num_lid_entries,
						ZOLTAN_ID_PTR global_ids,
						ZOLTAN_ID_PTR local_ids,
						int wgt_dim, float *obj_wgts,
						int *ierr)
{
  Mesh *mesh = (Mesh *) data;
  *ierr = ZOLTAN_OK;

  uint i = 0;
  for (VertexIterator vertex(*mesh); !vertex.end(); ++vertex)
  {
    global_ids[i] = mesh->distdata().get_global(*vertex);
    local_ids[i++] = vertex->index();
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

  Mesh *mesh = (Mesh *) data;

  if (num_obj != mesh->numVertices())
  {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  uint i = 0;
  for (VertexIterator vertex(*mesh); !vertex.end(); ++vertex)
  {
    geom_vec[i] = vertex->point().x();
    geom_vec[i + 1] = vertex->point().y();
    if (num_dim > 2)
      geom_vec[i + 2] = vertex->point().z();
  }

  *ierr = ZOLTAN_OK;
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
