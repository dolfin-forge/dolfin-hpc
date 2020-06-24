// Copyright (C) 2015 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/ZoltanInterface.h>

#include <dolfin/config/dolfin_config.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshRenumber.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/parameters.h>

#ifdef HAVE_ZOLTAN
#include <zoltan_cpp.h>
#endif


using namespace dolfin;

#ifdef HAVE_ZOLTAN
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionCommonZoltan(Mesh& mesh,
                                            MeshValues<uint, Cell>& partitions,
                                            MeshValues<uint, Cell>* weight)
{

#ifndef ENABLE_P1_OPTIMIZATIONS
  Zoltan *zz_ = new Zoltan(MPI::DOLFIN_COMM);

  // Setup query functions for graph based partitioning
  zz_->Set_Num_Obj_Fn(partitionZoltanNumCells, &mesh);
  zz_->Set_Num_Edges_Multi_Fn(partitionZoltanNumEdges, &mesh);
  zz_->Set_Obj_List_Fn(partitionZoltanCellList, &mesh);
  zz_->Set_Edge_List_Multi_Fn(partitionZoltanEdgeList, &mesh);


  // Use Zoltan's Parallel Hypergraph and Graph partitioner
  zz_->Set_Param( "LB_METHOD", "GRAPH");

  partitions = MPI::rank();

  partitionZoltanInternal(mesh, partitions, zz_);

  delete zz_;

#else
  error("Zoltan's graph partitioner needs shared facet connectivity");
#endif
}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionGeomZoltan(Mesh& mesh,
                                          MeshValues<uint, Vertex>& partitions)
{

  Zoltan *zz_ = new Zoltan(MPI::DOLFIN_COMM);

  // Setup query functions for geometry based partitioning
  zz_->Set_Num_Obj_Fn(partitionZoltanNumVertices, &mesh);
  zz_->Set_Obj_List_Fn(partitionZoltanVertexList, &mesh);
  zz_->Set_Num_Geom_Fn(partitionZoltanNumGeom, &mesh);
  zz_->Set_Geom_Multi_Fn(partitionZoltanGeomCoords, &mesh);

  // Use a Hilbert Space-Filling Curve
  zz_->Set_Param( "LB_METHOD", "HSFC");

  partitions = MPI::rank();

  partitionZoltanInternal(mesh, partitions, zz_);

  delete zz_;
}
//-----------------------------------------------------------------------------
template <class Entity>
void ZoltanInterface::partitionZoltanInternal(Mesh& mesh,
					      MeshValues<uint, Entity> & partitions,
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
    partitions(export_local_ids[i]) = export_procs[i];
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

  return mesh->num_cells();

}
//-----------------------------------------------------------------------------
int ZoltanInterface::partitionZoltanNumVertices(void *data, int *ierr)
{

  Mesh *mesh = (Mesh *) data;
  *ierr = ZOLTAN_OK;

  return mesh->size(0);

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
    global_ids[i] = cell->global_index();
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
    global_ids[i] = vertex->global_index();
    local_ids[i++] = vertex->index();
  }

  return;

}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionZoltanNumEdges(void *data,
					      int num_gid_entries,
					      int num_lid_entries,
					      int num_obj,
					      ZOLTAN_ID_PTR global_ids,
					      ZOLTAN_ID_PTR local_ids,
					      int *num_edges,
					      int *ierr)
{
  Mesh *mesh = (Mesh *) data;
  uint const tdim = mesh->topology_dimension();
  if (num_obj != mesh->num_cells() || num_gid_entries > 1)
  {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  // Count number of dual graph edges
  uint i = 0;
  for (CellIterator c(*mesh); !c.end(); i++, ++c)
  {
    num_edges[i] = 0;
    for (FacetIterator f(*c); !f.end(); ++f)
    {
      // Filter out non-shared boundary facets
      if (f->num_entities(tdim) == 1 && !f->is_shared())
      {
        continue;
      }
      num_edges[i]++;
    }
  }


  *ierr = ZOLTAN_OK;
}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionZoltanEdgeList(void *data, int num_gid_entries,
					      int num_lid_entries, int num_obj,
					      ZOLTAN_ID_PTR global_ids,
					      ZOLTAN_ID_PTR local_ids,
					      int *num_edges,
					      ZOLTAN_ID_PTR nbor_global_id,
					      int *nbor_procs, int wgt_dim,
					      float *ewgts, int *ierr)
{

  Mesh *mesh = (Mesh *) data;
  uint const facet_dim = mesh->type().facet_dim();
  uint const tdim = mesh->type().dim();
  DistributedData const& dist = mesh->distdata()[facet_dim];
  if (num_obj != mesh->num_cells() || num_gid_entries > 1)
  {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  // Map between global facet index and neigh. Cell (global index)
  _map<uint, uint> facet_cell_map;

  uint rank = MPI::rank();
  uint pe_size = MPI::size();

  Array<uint> *glb_facet = new Array<uint>[pe_size];
  for (SharedIterator f(dist); f.valid(); ++f)
  {
    uint const adj_rank = *(f.adj().begin());
    glb_facet[adj_rank].push_back(f.global_index());
  }

  uint max_recv = 0;
  for (int i = 0; i < pe_size; i++)
  {
    max_recv = std::max(max_recv, (uint) glb_facet[i].size());
  }
  uint *recv_buff = new uint[max_recv];
  uint *send_buff = new uint[max_recv];

  for (int j = 1; j < (int) pe_size; ++j)
  {
    int src = (rank - j + pe_size) % pe_size;
    int dest = (rank + j) % pe_size;

    // Exchange global facet index
    int recv_count = MPI::sendrecv( &glb_facet[dest][0], glb_facet[dest].size(),
                                    dest, recv_buff, max_recv, src, 1 );

    for (uint k = 0; k < recv_count; k++)
    {
      Facet f(*mesh, dist.get_local(recv_buff[k]));
      send_buff[k] = mesh->distdata()[tdim].get_global(f.entities(tdim)[0]);
    }

    // Send back corresponding global cell index
    recv_count = MPI::sendrecv( &send_buff[0], recv_count, src,
                                recv_buff, max_recv, dest, 2 );

    for (uint k = 0; k < recv_count; k++)
    {
      facet_cell_map[glb_facet[dest][k]] = recv_buff[k];
    }
  }

  delete[] recv_buff;
  delete[] send_buff;

  uint neigh_idx= 0;
  uint i = 0;

  // Construct dual graph
  for (CellIterator c(*mesh); !c.end(); i++, ++c)
  {
    for (FacetIterator f(*c); !f.end(); ++f)
    {
      // Filter out non-shared boundary facets
      if (f->num_entities(tdim) == 1 && !f->is_shared())
      {
        continue;
      }
      else if (f->is_shared())
      {
        nbor_global_id[i] = facet_cell_map[f->global_index()];
      }
      else
      {
        if (f->entities(tdim)[0] != c->index())
        {
          neigh_idx = 0;
        }
        else if (f->entities(tdim)[1] != c->index())
        {
          neigh_idx = 1;
        }

        Cell neig_cell(*mesh, f->entities(tdim)[0]);
        nbor_global_id[i] = neig_cell.global_index();
      }
    }
  }

  *ierr = ZOLTAN_OK;
}
//-----------------------------------------------------------------------------
int ZoltanInterface::partitionZoltanNumGeom(void *data, int *ierr)
{
  Mesh *mesh = (Mesh *) data;
  *ierr = ZOLTAN_OK;

  return mesh->geometry_dimension();

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

  if (num_obj != mesh->size(0))
  {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (VertexIterator vertex(*mesh); !vertex.end(); geom_vec+=num_dim, ++vertex)
  {
    std::copy(vertex->x(), vertex->x() + num_dim, geom_vec);
  }

  *ierr = ZOLTAN_OK;
}
//-----------------------------------------------------------------------------
#else
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionCommonZoltan(Mesh&,
                                            MeshValues<uint, Cell>&,
                                            MeshValues<uint, Cell>*)
{
  error("DOLFIN needs to be built with Zoltan support");
}
//-----------------------------------------------------------------------------
void ZoltanInterface::partitionGeomZoltan(Mesh&, MeshValues<uint, Vertex> &)
{
  error("DOLFIN needs to be built with Zoltan support");
}
//-----------------------------------------------------------------------------
#endif
