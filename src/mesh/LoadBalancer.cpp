// Copyright (C) 2008 Niclas Jansson.
// Licensed under the GNU LGPL Version 2.1.

#include <dolfin/mesh/LoadBalancer.h>

#include <dolfin/common/Array.h>
#include <dolfin/main/MPI.h>
#include <dolfin/main/PE.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/CellIterator.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/EdgeIterator.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/parameters.h>

#include <string>
#include <cstring>

using namespace dolfin;

_ordered_map<Mesh *, MeshValues<uint, Cell> *> LoadBalancer::s_;

#ifdef HAVE_MPI
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<uint, Cell>& weight)
{
  begin("Load balancing with cell weights");

  // Repartition mesh
  message("Partition");
  MeshValues<uint, Cell> partitions(mesh);
  mesh.partition(partitions, weight);

  // Calculate process reassignment
  message("Process reassignment");
  uint max_sendrecv;
  process_reassignment(partitions, &max_sendrecv);

  // Distribute mesh according to new partition function
  message("Distribute");
  mesh.distribute(partitions);

  end();
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                           Type type)
{
  balance(mesh, cell_marker, 0.0, 0.0, 0.0, type);
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                           real tf, real tb, real ts, Type type)
{
  // Construct weight function
  MeshValues<uint, Cell> weight(mesh);
  uint w_local, w_sum, w_max;
  real w_avg;
  message("Estimate load imbalance");
  weight_function(mesh, cell_marker, weight, &w_local, type);

  // Preliminary evalution of load imbalance
  MPI::all_reduce<MPI::max>( w_local, w_max );
  MPI::all_reduce<MPI::sum>( w_local, w_sum );

  w_avg = (real) w_sum / (real) MPI::size();

  const real threshold = 1.04;
  real imbalance = (real) w_max / (real) w_avg;
  if (threshold > imbalance)
  {
    message("Load imbalance %0.2f percent, below threshold.",
            (imbalance - 1.0) * 100);
    return;
  }
  else
  {
    message("Repartitioning %0.2f percent load imbalance.",
            (imbalance - 1.0) * 100);
  }

  // Repartition mesh
  MeshValues<uint, Cell> partitions(mesh);
  mesh.partition(partitions, weight);

  // Calculate process reassignment
  uint max_sendrecv;
  message("Process reassignment");
  process_reassignment(partitions, &max_sendrecv);

  // Determine if the computation gains from repartitioning/reassignment
  if (tf > 0.0 && tb > 0.0 && ts > 0.0)
  {
    if (!computational_gain(mesh, weight, partitions, max_sendrecv, tf, tb, ts))
    {
      return;
    }
  }

  // Distribute mesh according to new partition function
  if (dolfin_get<bool>("Load balancer redistribute"))
  {
    MeshData D(mesh); D.add(cell_marker);
    message("Redistribute mesh");
    mesh.distribute(partitions, D);
  }
  else
  {
    swap( LoadBalancer::partitions(mesh), partitions );
  }

  if (dolfin_get<bool>("Load balancer report"))
  {
    weight_function(mesh, cell_marker, weight, &w_local, type);

  MPI::all_reduce<MPI::max>( w_local, w_max );
  MPI::all_reduce<MPI::sum>( w_local, w_sum );

    w_avg = (real) w_sum / (real) MPI::size();
    real new_imbalance = (real) w_max / (real) w_avg;
    message("%0.2f percent load imbalance after repartitioning.",
            (new_imbalance - 1.0) * 100);
  }
}

//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<uint, Cell>& weight,
                           Array<VertexFunctionPPair>& vertex_functions)
{
  begin("Load balancing with cell weights and vertex functions");

  // Repartition mesh
  MeshValues<uint, Cell> partitions(mesh);
  mesh.partition(partitions, weight);

  // Calculate process reassignment
  uint max_sendrecv;
  process_reassignment(partitions, &max_sendrecv);

  // Distribute mesh according to new partition function
  MeshData D(mesh);
  for (Array<VertexFunctionPPair>::iterator it = vertex_functions.begin();
       it != vertex_functions.end(); ++it)
  {
    it->second = it->first; D.add(*(it->second));
  }
  mesh.distribute(partitions, D);

  end();
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                           Array<VertexFunctionPPair>& vertex_functions,
                           Type type)
{
  balance(mesh, cell_marker, vertex_functions, 0.0, 0.0, 0.0, type);
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                           Array<VertexFunctionPPair>& vertex_functions,
                           real tf, real tb, real ts, Type type)
{

  // Construct weight function
  MeshValues<uint, Cell> weight(mesh);
  uint w_local, w_sum, w_max;
  real w_avg;
  weight_function(mesh, cell_marker, weight, &w_local, type);

  // Preliminary evalution of load imbalance
  MPI::all_reduce<MPI::max>( w_local, w_max );
  MPI::all_reduce<MPI::sum>( w_local, w_sum );

  w_avg = (real) w_sum / (real) MPI::size();

  const real threshold = 1.04;
  real imbalance = (real) w_max / (real) w_avg;
  if (threshold > imbalance)
  {
    message("Load imbalance %0.2f percent, below threshold.",
            (imbalance - 1.0) * 100);
    return;
  }
  else
  {
    message("Repartitioning %0.2f percent load imbalance.",
            (imbalance - 1.0) * 100);
  }

  // Repartition mesh
  MeshValues<uint, Cell> partitions(mesh);
  mesh.partition(partitions, weight);

  // Calculate process reassignment
  uint max_sendrecv;
  process_reassignment(partitions, &max_sendrecv);

  // Determine if the computation gains from repartitioning/reassignment
  if (tf > 0.0 && tb > 0.0 && ts > 0.0)
  {
    if (!computational_gain(mesh, weight, partitions, max_sendrecv, tf, tb, ts))
    {
      return;
    }
  }

  // Distribute mesh according to new partition function
  if (dolfin_get<bool>("Load balancer redistribute"))
  {
    MeshData D(mesh);
    D.add(cell_marker);
    for (Array<VertexFunctionPPair>::iterator it = vertex_functions.begin();
         it != vertex_functions.end(); ++it)
    {
      *(it->second) = *(it->first); D.add(*(it->second));
    }
    mesh.distribute(partitions, D);
  }
  else
  {
    swap( LoadBalancer::partitions(mesh), partitions );
  }

  if (dolfin_get<bool>("Load balancer report"))
  {
    weight_function(mesh, cell_marker, weight, &w_local, type);

    MPI::all_reduce<MPI::max>( w_local, w_max );
    MPI::all_reduce<MPI::sum>( w_local, w_sum );

    w_avg = (real) w_sum / (real) MPI::size();
    real new_imbalance = (real) w_max / (real) w_avg;
    message("%0.2f percent load imbalance after repartitioning.",
            (new_imbalance - 1.0) * 100);
  }
}

//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<uint, Cell>& weight,
                           Array<CellFunctionPPair>& cell_functions,
                           Array<VertexFunctionPPair>& vertex_functions)
{
  begin("Load balancing with cell weights and cell/vertex functions");

  // Repartition mesh
  MeshValues<uint, Cell> partitions(mesh);
  mesh.partition(partitions, weight);

  // Calculate process reassignment
  uint max_sendrecv;
  process_reassignment(partitions, &max_sendrecv);

  // Distribute mesh according to new partition function
  MeshData D(mesh);
  for (Array<CellFunctionPPair>::iterator it = cell_functions.begin();
       it != cell_functions.end(); ++it)
  {
    *(it->second) = *(it->first); D.add(*(it->second));
  }
  for (Array<VertexFunctionPPair>::iterator it = vertex_functions.begin();
       it != vertex_functions.end(); ++it)
  {
    *(it->second) = *(it->first); D.add(*(it->second));
  }
  mesh.distribute(partitions, D);

  end();
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                           Array<CellFunctionPPair>& cell_functions,
                           Array<VertexFunctionPPair>& vertex_functions,
                           Type type)
{
  balance(mesh, cell_marker, cell_functions, vertex_functions, 0.0, 0.0, 0.0,
          type);
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                           Array<CellFunctionPPair>& cell_functions,
                           Array<VertexFunctionPPair>& vertex_functions,
                           real tf, real tb, real ts, Type type)
{
  // Construct weight function
  MeshValues<uint, Cell> weight(mesh);
  uint w_local, w_sum, w_max;
  real w_avg;
  weight_function(mesh, cell_marker, weight, &w_local, type);

  // Preliminary evalution of load imbalance
  MPI::all_reduce<MPI::max>( w_local, w_max );
  MPI::all_reduce<MPI::sum>( w_local, w_sum );

  w_avg = (real) w_sum / (real) MPI::size();

  const real threshold = 1.04;
  real imbalance = (real) w_max / (real) w_avg;
  if (threshold > imbalance)
  {
    message("Load imbalance %0.2f percent, below threshold.",
            (imbalance - 1.0) * 100);
    return;
  }
  else
  {
    message("Repartitioning %0.2f percent load imbalance.",
            (imbalance - 1.0) * 100);
  }

  // Repartition mesh
  MeshValues<uint, Cell> partitions(mesh);
  mesh.partition(partitions, weight);

  // Calculate process reassignment
  uint max_sendrecv;
  process_reassignment(partitions, &max_sendrecv);

  // Determine if the computation gains from repartitioning/reassignment
  if (tf > 0.0 && tb > 0.0 && ts > 0.0)
  {
    if (!computational_gain(mesh, weight, partitions, max_sendrecv, tf, tb, ts))
    {
      return;
    }
  }

  // Distribute mesh according to new partition function
  if (dolfin_get<bool>("Load balancer redistribute"))
  {
    MeshData D(mesh);
    D.add(cell_marker);
    for (Array<CellFunctionPPair>::iterator it = cell_functions.begin();
         it != cell_functions.end(); ++it)
    {
      *(it->second) = *(it->first); D.add(*(it->second));
    }
    for (Array<VertexFunctionPPair>::iterator it = vertex_functions.begin();
         it != vertex_functions.end(); ++it)
    {
      *(it->second) = *(it->first); D.add(*(it->second));
    }
    mesh.distribute(partitions, D);
  }
  else
  {
    swap( LoadBalancer::partitions(mesh), partitions );
  }

  if (dolfin_get<bool>("Load balancer report"))
  {
    weight_function(mesh, cell_marker, weight, &w_local, type);

    MPI::all_reduce<MPI::max>( w_local, w_max );
    MPI::all_reduce<MPI::sum>( w_local, w_sum );

    w_avg = (real) w_sum / (real) MPI::size();
    real new_imbalance = (real) w_max / (real) w_avg;
    message("%0.2f percent load imbalance after repartitioning.",
            (new_imbalance - 1.0) * 100);
  }
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<uint, Cell>& weight,
                           Array<CellFunctionPPair>& cell_functions)
{
  begin("Load balancing with cell weights and cell functions");

  // Repartition mesh
  MeshValues<uint, Cell> partitions(mesh);
  mesh.partition(partitions, weight);

  // Calculate process reassignment
  uint max_sendrecv;
  process_reassignment(partitions, &max_sendrecv);

  // Distribute mesh according to new partition function
  MeshData D(mesh);
  for (Array<CellFunctionPPair>::iterator it = cell_functions.begin();
       it != cell_functions.end(); ++it)
  {
    *(it->second) = *(it->first); D.add(*(it->second));
  }
  mesh.distribute(partitions, D);

  end();
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                           Array<CellFunctionPPair>& cell_functions, Type type)
{
  balance(mesh, cell_marker, cell_functions, 0.0, 0.0, 0.0, type);
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh& mesh, MeshValues<bool, Cell>& cell_marker,
                           Array<CellFunctionPPair>& cell_functions, real tf,
                           real tb, real ts, Type type)
{
  // Construct weight function
  MeshValues<uint, Cell> weight(mesh);
  uint w_local, w_sum, w_max;
  real w_avg;
  weight_function(mesh, cell_marker, weight, &w_local, type);

  // Preliminary evalution of load imbalance
  MPI::all_reduce<MPI::max>( w_local, w_max );
  MPI::all_reduce<MPI::sum>( w_local, w_sum );

  w_avg = (real) w_sum / (real) MPI::size();

  const real threshold = 1.04;
  real imbalance = (real) w_max / (real) w_avg;
  if (threshold > imbalance)
  {
    message("Load imbalance %0.2f percent, below threshold.",
            (imbalance - 1.0) * 100);
    return;
  }
  else
  {
    message("Repartitioning %0.2f percent load imbalance.",
            (imbalance - 1.0) * 100);
  }

  // Repartition mesh
  MeshValues<uint, Cell> partitions(mesh);
  mesh.partition(partitions, weight);

  // Calculate process reassignment
  uint max_sendrecv;
  process_reassignment(partitions, &max_sendrecv);

  // Determine if the computation gains from repartitioning/reassignment
  if (tf > 0.0 && tb > 0.0 && ts > 0.0)
  {
    if (!computational_gain(mesh, weight, partitions, max_sendrecv, tf, tb, ts))
    {
      return;
    }
  }

  // Distribute mesh according to new partition function
  if (dolfin_get<bool>("Load balancer redistribute"))
  {
    MeshData D(mesh);
    D.add(cell_marker);
    for (Array<CellFunctionPPair>::iterator it = cell_functions.begin();
         it != cell_functions.end(); ++it)
    {
      it->second = it->first; D.add(*(it->second));
    }
    mesh.distribute(partitions, D);
  }
  else
  {
    swap( LoadBalancer::partitions(mesh), partitions );
  }

  if (dolfin_get<bool>("Load balancer report"))
  {
    weight_function(mesh, cell_marker, weight, &w_local, type);

    MPI::all_reduce<MPI::max>( w_local, w_max );
    MPI::all_reduce<MPI::sum>( w_local, w_sum );

    w_avg = (real) w_sum / (real) MPI::size();
    real new_imbalance = (real) w_max / (real) w_avg;
    message("%0.2f percent load imbalance after repartitioning.",
            (new_imbalance - 1.0) * 100);
  }
}
//---------------------------------------------------------------------------
void LoadBalancer::weight_function(Mesh& mesh,
                                   MeshValues<bool, Cell>& cell_marker,
                                   MeshValues<uint, Cell>& weight,
                                   uint* w_sum, Type type)
{
  // Calculated weights for the dual graph with a simulated mesh refinement
  real max, l;
  uint index = 0;
  weight = 1;
  //  *w_sum = mesh.numCells();

  // Make sure cell - edge connectivity is created, note we have to
  // explicitly create the connectivity here since not all ranks will
  // call the EdgeIterator, thus causing a deadlock in MeshRenumber
  mesh.init(mesh.type().dim(), 1);

  if (type == Default)
  {
    MeshValues<bool, Cell> used_cell(mesh);
    MeshValues<bool, Edge> used_edge(mesh);
    used_cell = false;
    used_edge = false;

    for (CellIterator c(mesh); !c.end(); ++c)
    {
      if (cell_marker(*c) && !used_cell(*c))
      {
        max = 0.0;
        for (EdgeIterator e(*c); !e.end(); ++e)
        {
          if (!used_edge(*e))
          {
            l = e->length();
            if (max < l)
            {
              max = l;
              index = e->index();
            }
          }
        }
        if (max == 0.0) continue;
        Edge le(mesh, index);
        for (CellIterator nc(le); !nc.end(); ++nc)
        {
          if (!used_cell(*nc))
          {
            //	  *w_sum++;
            weight(*nc) = weight(*nc) + 1;
            used_cell(*nc) = true;
            for(EdgeIterator e(*nc); !e.end(); ++e)
              used_edge(*e) = true;
          }
        }
      }
    }
  }
  else if (type == LEPP)
  {
    for (CellIterator c(mesh); !c.end(); ++c)
    {
      if (cell_marker(*c))
      {
        max = 0.0;
        for (EdgeIterator e(*c); !e.end(); ++e)
        {
          l = e->length();
          if (max < l)
          {
            max = l;
            index = e->index();
          }
        }
        Edge le(mesh, index);
        for (CellIterator oc(le); !oc.end(); ++oc)
          weight_lepp(mesh, *oc, le, weight, 0);
      }
    }
  }
  else if (type == EdgeCollapse)
  {
    for (CellIterator c_it(mesh); !c_it.end(); ++c_it)
    {
      if (cell_marker(*c_it))
      {
        // cell marked for coarsening gets increased weight
        weight(*c_it) = 2u;

        // all neighboring cells also get increased weight
        for ( CellIterator nc_it(*c_it) ; !nc_it.end() ; ++nc_it )
          weight(*nc_it) = 2u;
      }
    }
  }
  else
  {
    error("Unknown Type for LoadBalancer.");
  }

  *w_sum = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
    *w_sum += weight(*c);

}
//-----------------------------------------------------------------------------
void LoadBalancer::weight_lepp(Mesh& mesh, Cell& c, Edge& ce,
                               MeshValues<uint, Cell>& weight, uint depth)
{
  // weight(c.index()) = weight(c.index()) + 1;
  real l;
  real max = 0.0;
  uint index = 0;
  for (EdgeIterator e(c); !e.end(); ++e)
  {
    l = e->length();
    if (max < l)
    {
      max = l;
      index = e->index();
    }
  }

  Edge le(mesh, index);

  if (le.index() == ce.index() || depth > 1) return;

  weight(c.index()) = weight(c.index()) + 1;

  depth++;

  for (CellIterator nc(le); !nc.end(); ++nc)
    weight_lepp(mesh, *nc, le, weight, depth);
}
//-----------------------------------------------------------------------------
void LoadBalancer::process_reassignment(MeshFunction<uint>& partitions,
                                        uint* max_sendrecv)
{
  uint pe_size = MPI::size();
  uint rank = MPI::rank();

  // Construct my row in the similarity matrix
  uint *sim_row = new uint[pe_size];
  uint *sorted_indices = new uint[pe_size];
  memset(sim_row, 0, pe_size * sizeof(uint));
  for (uint i = 0; i < partitions.size(); i++)
    sim_row[partitions(i)]++;

  // Gather similarity matrix on root node
  uint m = pe_size * pe_size;
  uint *SiM = new uint[m]; // FIXME, this should only be needon on rank == 0
  //  MPI_Gather(sim_row, pe_size, MPI_UNSIGNED,
  //	     SiM, pe_size, MPI_UNSIGNED, 0, MPI::DOLFIN_COMM);

  MPI::all_gather( sim_row, pe_size, SiM, pe_size );

  uint *sorted = new uint[m];
  pradixsort_matrix(&sorted_indices[0], &SiM[0], pe_size);
  MPI::check_error( MPI_Gather(sorted_indices, pe_size, MPI_UNSIGNED, sorted,
                               pe_size, MPI_UNSIGNED, 0, MPI::DOLFIN_COMM) );

  uint *map = new uint[pe_size];

  if (rank == 0)
  {

    uint num_assigned = 0;
    uint idx, idy;
    bool *unassigned_x = new bool[pe_size];
    bool *unassigned_y = new bool[pe_size];
    for (uint i = 0; i < pe_size; i++)
    {
      unassigned_y[i] = true;
      unassigned_x[i] = true;
    }

    uint tail = m;
    // Greedy approximation for the MWBG problem
    while (num_assigned < pe_size)
    {
      for (uint i = tail; i > 0; i--)
      {
        idx = (sorted[i - 1]) / pe_size;
        idy = (sorted[i - 1]) % pe_size;
        if (unassigned_x[idx] && unassigned_y[idy])
        {
          unassigned_y[idy] = false;
          unassigned_x[idx] = false;
          map[idx] = idy;
          tail = i;
          num_assigned++;
          break;
        }
      }
    }

    delete[] unassigned_x;
    delete[] unassigned_y;
  }

  MPI::check_error( MPI_Bcast(map, pe_size, MPI_UNSIGNED, 0, MPI::DOLFIN_COMM) );

  // Reassign processors
  for (uint i = 0; i < partitions.size(); i++)
    partitions(i) = map[partitions(i)];

  // Calculate maximum number to send from processor
  *max_sendrecv = 0;
  for (uint i = 0; i < pe_size; i++)
    if (map[i] != rank) *max_sendrecv = std::max(*max_sendrecv, sim_row[i]);

  delete[] map;
  delete[] sim_row;
  delete[] sorted_indices;
  delete[] sorted;
  delete[] SiM; // FIXME, this should only be needon on rank == 0
}
//-----------------------------------------------------------------------------
bool LoadBalancer::computational_gain(Mesh& mesh,
                                      MeshValues<uint, Cell>& weight,
                                      MeshValues<uint, Cell>& partitions,
                                      uint max_sendrecv, real tf, real tb,
                                      real ts)
{
  uint pe_size = MPI::size();

  uint w_old = 0;
  uint *tmp_w = new uint[pe_size];
  memset(tmp_w, 0, pe_size * sizeof(uint));
  for (CellIterator c(mesh); !c.end(); ++c)
  {
    w_old += weight(*c);
    tmp_w[partitions(*c)] += weight(*c);
  }

  uint w_new;
  for (uint i = 0; i < pe_size; i++)
    MPI::check_error( MPI_Reduce(&tmp_w[i], &w_new, 1, MPI_UNSIGNED, MPI_SUM, i,
                                 MPI::DOLFIN_COMM) );

  uint w_oldmax, w_newmax;
  MPI::all_reduce<MPI::max>( w_old, w_oldmax );
  MPI::all_reduce<MPI::max>( w_new, w_newmax );

  uint tmp = max_sendrecv; // No MPI aliasing on BG/L
  MPI::all_reduce<MPI::max>( tmp, max_sendrecv );
  message("**** %d - %d = %d  w_nmax / w_omax = %f maxsr = %d ****", w_oldmax,
          w_newmax, w_oldmax - w_newmax, (real) w_newmax / (real) w_oldmax,
          max_sendrecv);

  real ndims = (real) mesh.type().num_vertices(mesh.topology_dimension());
  real comp_cost = tf * ((ndims * ndims) * ((real) (w_oldmax - w_newmax)));
  //  comp_cost +=  ts + tb *((real) (w_oldmax - w_newmax));
  real comm_cost = ts + tb * (ndims * (real) max_sendrecv);

  message("comp_cost %g comm_cost %g", comp_cost, comm_cost);

  delete[] tmp_w;
  return (comp_cost > comm_cost);
}
//-----------------------------------------------------------------------------
void LoadBalancer::radixsort_matrix(uint* res, uint* Matrix, uint m, bool desc)
{

  uint *index = new uint[m];
  uint *tmp = new uint[m];

  for (uint i = 0; i < m; i++)
    index[i] = i;

  uint *count = new uint[256];
  uint *map = new uint[256];
  for (uint i = 0; i < 4; i++)
  {
    memcpy(tmp, index, m * sizeof(uint));
    memset(count, 0, 256 * sizeof(uint));
    for (uint j = 0; j < m; j++)
      count[((Matrix[tmp[j]]) >> (8 * i)) & 0xff]++;

    map[0] = 0;
    for (uint j = 1; j < 256; j++)
      map[j] = map[j - 1] + count[j - 1];

    for (uint j = 0; j < m; j++)
      index[map[((Matrix[tmp[j]]) >> (8 * i)) & 0xff]++] = tmp[j];
  }

  if (desc)
  {
    uint j = m - 1;
    for (uint i = 0; i < m; i++)
      res[i] = index[j--];
  }
  else
  {
    memcpy(res, index, m * sizeof(uint));
  }

  delete[] map;
  delete[] count;
  delete[] tmp;
  delete[] index;
}
//-----------------------------------------------------------------------------
void LoadBalancer::pradixsort_matrix(uint* res, uint* Matrix, uint m)
{

  uint *index = new uint[m];
  uint *tmp = new uint[m];

  uint rank = MPI::rank();
  uint pe_size = MPI::size();

  uint *recvbuff = new uint[2 * pe_size];
  Array<uint> *sendbuff = new Array<uint> [pe_size];

  for (uint i = 0; i < m; i++)
    index[i] = i + rank * pe_size;

  uint *count = new uint[256];
  uint *glb_count = new uint[256];
  uint *map = new uint[256];
  uint offset, glb_scan;
  for (uint i = 0; i < 4; i++)
  {
    memcpy(tmp, index, m * sizeof(uint));
    memset(count, 0, 256 * sizeof(uint));
    for (uint j = 0; j < m; j++)
      count[((Matrix[tmp[j]]) >> (8 * i)) & 0xff]++;

    MPI::all_reduce<MPI::sum>( count, glb_count, 256 );

    offset = 0;
    for (uint j = 0; j < 256; j++)
    {
      MPI::check_error( MPI_Scan(&count[j], &glb_scan, 1, MPI_UNSIGNED, MPI_SUM,
                                 MPI::DOLFIN_COMM) );
      map[j] = glb_scan + offset;
      offset += glb_count[j];
    }

    for (uint j = 0; j < 256; j++)
      map[j]--;

    for (uint j = m; j > 0; j--)
    {
      uint glb_index = map[((Matrix[tmp[j - 1]]) >> (8 * i)) & 0xff]--;
      uint target = floor((double) glb_index / (double) pe_size);
      if (target != rank)
      {
        sendbuff[target].push_back(glb_index);
        sendbuff[target].push_back(tmp[j - 1]);
      }
      else
      {
        index[glb_index % pe_size] = tmp[j - 1];
      }
    }
    uint src, dest;
    int recv_size;

    for (uint j = 1; j < pe_size; j++)
    {

      src = (rank - j + pe_size) % pe_size;
      dest = (rank + j) % pe_size;

      recv_size = MPI::sendrecv( &sendbuff[dest][0], sendbuff[dest].size(), dest,
                                 recvbuff, 2 * pe_size, src, 0 );

      for (int k = 0; k < recv_size; k += 2)
      {
        index[recvbuff[k] % pe_size] = recvbuff[k + 1];
      }
      sendbuff[dest].clear();
    }
  }

  memcpy(res, index, m * sizeof(uint));

  for (uint i = 0; i < pe_size; i++)
    sendbuff[i].clear();
  delete[] sendbuff;

  delete[] map;
  delete[] glb_count;
  delete[] count;
  delete[] recvbuff;
  delete[] tmp;
  delete[] index;

}
//-----------------------------------------------------------------------------
#else
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh&, MeshValues<uint, Cell>&)
{
  warning("Load balancing only implemented for MPI");
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh&, MeshValues<bool, Cell>&, Type)
{
  warning("Load balancing only implemented for MPI");
}
//-----------------------------------------------------------------------------
void LoadBalancer::balance(Mesh&, MeshValues<bool, Cell>&,
                           real, real, real, Type)
{
  warning("Load balancing only implemented for MPI");
}
//-----------------------------------------------------------------------------
#endif
//-----------------------------------------------------------------------------
MeshValues<uint, Cell>& LoadBalancer::partitions(Mesh& mesh)
{
  _ordered_map<Mesh *, MeshValues<uint, Cell> *>::iterator it = s_.find(&mesh);
  if (it == s_.end())
  {
    MeshValues<uint, Cell> * v = new MeshValues<uint, Cell>(mesh, PE::rank());
    s_[&mesh] = v;
    return *v;
  }
  return *(it->second);
}
//-----------------------------------------------------------------------------
bool LoadBalancer::clear(Mesh& mesh)
{
  _ordered_map<Mesh *, MeshValues<uint, Cell> *>::iterator it = s_.find(&mesh);
  if (it != s_.end())
  {
    delete it->second;
    s_.erase(it);
    return true;
  }
  return false;
}
//-----------------------------------------------------------------------------

